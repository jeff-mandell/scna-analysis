library(data.table)
library(pbapply)
library(ggplot2)
library(ggrepel)
library(devtools)

# Should be running in RStudio with the dev version cancereffectsizeR set as the active project.
usethis::proj_set('~/cancereffectsizeR/') # change path as necessary
setwd(usethis::proj_get())
load_all()
load_all('~/ces.refset.hg38') 
stopifnot(packageVersion('cancereffectsizeR') > as.package_version('3.0.0.9000'))
stopifnot(packageVersion('ces.refset.hg38') > as.package_version('1.3.0'))

# Load output from prep_SCNA_analysis.R; adjust path as needed
output_dir = '~/Desktop/repos/scna-analysis/untracked/outputs/'
prepped_data = readRDS(file.path(output_dir, 'prepped_cna.rds'))
cna_calls = copy(prepped_data$prepped_calls$calls)

event_type = 'increase' # increases for now
seg_rates = get_seg_increase_rates(cna_calls = cna_calls,
                                               increase_burden = prepped_data$cna_burdens$increase_burden)

# We take our gene coordinates and get rid of some genes so that no remaining genes overlap.
# See the function itself for the approach used.
gene_coord = get_disjoint_gene_coord(ces.refset.hg38$cancer_gene_coord)

## run_round will run a collection of inferences across various loci in a chromosome
## Supply fixed_genes, which must be all in the same chromosome, or a chromosome (chr).
# The genes in fixed_genes get used in every inference in the round. The remaining same-chromosome
# genes that present in gene_coord will be cycled through. That is, each inference will use all the 
# fixed genes and one additional gene.
# To speed runtime, set prop_noncancer to something <1, such as .1. By default, all TSG/oncogenes
# in gene_coord will be used, and the given proportion of additional genes will also be used.
run_round = function(fixed_genes = character(), chr = character(), gene_coord = gene_coord, 
                     prop_noncancer = 1) {
  if(! rlang::is_scalar_double(prop_noncancer) || prop_noncancer < 0 || prop_noncancer > 1) {
    stop('prop_noncancer must be scalar numeric on [0, 1].')
  }
  
  if(! is.data.table(gene_coord) || ! all(c('gene', 'cancer_anno', 'chr', 'start', 'end') %in% names(gene_coord))) {
    stop('gene_coord must be data.table with fields gene, cancer_anno, chr, start, end.')
  }
  
  if(! xor(identical(fixed_genes, character()), identical(chr, character()))) {
    stop('Must specify exactly one of chr and fixed_genes.')
  }
  
  if(! is.character(fixed_genes)) {
    stop('fixed_genes must be type character.')
  }
  if(any(! fixed_genes %in% gene_coord$gene)) {
    stop('All fixed_genes must be present in gene_coord.')
  }
  
  
  if(length(fixed_genes) > 0) {
    curr_chr = gene_coord[gene %in% fixed_genes, unique(chr)]
  } else {
    curr_chr = chr
  }
  stopifnot(length(curr_chr) == 1)
  valid_chr = unique(gene_coord$chr)
  if(! curr_chr %in% valid_chr) {
    stop('Specified chr not present in gene_coord.')
  }
  gene_coord = gene_coord[chr == curr_chr]

  noncancer = gene_coord[! cancer_anno %in% c('TSG', 'oncogene') & ! gene %in% fixed_genes][, .SD[.I %% floor(1/prop_noncancer) == 0]]
  gene_coord = rbind(gene_coord[cancer_anno %in% c('TSG', 'oncogene') | gene %in% fixed_genes], noncancer)
  
  other_genes = setdiff(gene_coord$gene, fixed_genes)
  output = rbindlist(pblapply(other_genes, 
                              function(x) {
                                run_seg_multi(
                                  selection_loci = gene_coord[gene %in% c(fixed_genes, x)],
                                  cna_calls = cna_calls,
                                  seg_rates = seg_rates
                                )
                              }, cl = 1), idcol = 'run')
  output[, is_fixed := gene %in% fixed_genes]
  
  fixed_output = unique(output[is_fixed == T, .(cancer_anno, start, end, si = median(si)), by = 'gene'])
  new_output = output[is_fixed == FALSE, .(gene, cancer_anno, start, end, si)]
  
  all_output = rbind(fixed_output, new_output)
  all_output[, chr_label := paste0('chr', curr_chr)]
  all_output[cancer_anno %in% c('TSG', 'oncogene') | gene %in% fixed_genes, gene_label := gene]
  all_output[, cancer_anno := factor(cancer_anno, levels = c('oncogene', 'other', 'TSG', 'noncancer'))]
  all_output[, midpoint := floor((start + end)/2)]
  return(all_output[])
}

# A plotting function
chr_coordinates = copy(prepped_data$prepped_calls$effective_coordinates)
chr_coordinates[, chr_label := paste0('chr', chr)]
plot_round = function(round_output) {
  chr_info = chr_coordinates[round_output, on = 'chr_label'][1]
  palette =  c(RColorBrewer::brewer.pal(3, "Spectral"), 'grey')
  gg = ggplot(round_output, aes(x = midpoint, y = si)) + 
    geom_rect(data = chr_info, aes(xmin = cen_start, xmax = cen_end, 
                  ymin = -Inf, ymax = Inf), fill = 'gray95',
              inherit.aes = FALSE) + 
    geom_point(aes(color = cancer_anno), size = 2.5) +
    xlab('Position') + ylab('Scaled selection for copy increase') + 
    scale_color_manual(name = 'Cancer annotation', values = palette) + 
    facet_wrap(~chr_label) +
    geom_text_repel(na.rm = T, aes(label = gene_label), size = 2, 
                    min.segment.length = 0, point.padding = .5, nudge_y = .07, segment.colour = 'gray20') +
    scale_x_continuous(labels = scales::label_comma()) +
    theme_light()
  return(gg)
}


r1 = run_round(chr = '7', gene_coord = gene_coord, prop_noncancer = .2)
ggr1 = plot_round(r1)

# After visual inspection, deciding to fix EGFR. We'll use the same genes as round 1.
r2 = run_round(fixed_genes = 'EGFR', gene_coord = gene_coord[gene %in% r1$gene])
ggr2 = plot_round(r2)

# IKZF1, a TSG, has one of the lowest effects. (CCT6A is even lower, but we'll ignore that for now.)
r3 = run_round(fixed_genes = c('EGFR', 'IKZF1'), gene_coord = gene_coord[gene %in% r1$gene])
ggr3 = plot_round(r3)

# CCT6A is still sticking out, so we'll use it next.
for_r4_choice = r3[! gene %in% c('EGFR', 'IKZF1')][, z := scale(si)][order(-abs(z))]

r4 = run_round(fixed_genes = c('EGFR', 'IKZF1', 'CCT6A'),
               gene_coord = gene_coord[gene %in% r1$gene])
r4[gene == 'CCT6A' | si > 5, gene_label := gene]
ggr4 = plot_round(r4)

## Combine all plots
all_plots = list(ggr1, ggr2, ggr3, ggr4)
all_plots = lapply(all_plots,
                   function(x) x + xlab('') + ylab('Effect') +
                     theme(strip.background = element_blank(),
                           strip.text = element_blank()))
grid = cowplot::plot_grid(plotlist = all_plots, ncol = 1,
                   labels = paste0('Round ', 1:4), label_size = 10,
                   label_x = 0, vjust = 0.3)
# Add empty plot for sufficient whitespace at the top
grid = cowplot::plot_grid(ggplot() + theme_void(), grid, rel_heights = c(.04, .96), nrow = 2)



