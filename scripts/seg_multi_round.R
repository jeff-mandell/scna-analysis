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

seg_rates = get_segmentation_rates(cna_calls = cna_calls, cna_burdens = prepped_data$cna_burdens)

### CNA rates by genomic interval (one chromosome at a time)
## Example 1: Tile the chromosome.
## gw = get_genomic_windows(refset = ces.refset.hg38, window_size = 1e5, annotate_intervals = TRUE)
## chr7_rates = get_cna_rates(chr = '7', cna_calls = cna_calls, seg_rates = seg_rates,
##                            rate_intervals = gw)


## Example 2: Use coordinates taken from gene transcripts
# We take our gene coordinates and get rid of some genes so that no remaining genes overlap.
# See the function itself for the approach used.
gc = get_disjoint_gene_coord(gene_coord = ces.refset.hg38$cancer_gene_coord)
gc = gc[, .(chr, start, end, gene1 = gene, gene1_anno = cancer_anno)]
gc[, rn := 1:.N]
subset = gc[rn %% 5 == 0 | gene1_anno %in% c('oncogene', 'TSG')]
gc[, rn := NULL]
chr7_rates = get_cna_rates(chr = '7', cna_calls = cna_calls, seg_rates = seg_rates, rate_intervals = subset)

# REQUIRED for current model: Get rid of sample-segments that cause problems for model.
chr7_rates = clean_rates(chr7_rates)


# Now we have all segments that can be used in every inference.

## run_round will run a collection of inferences across the input chromosomal intervals
# always_selected: data.table(); see example below.
# rates: from get_cna_rates(); for now, must also be processed as shown above.
run_round = function(always_selected = data.table(), rates = NULL, cores = 1, debug = FALSE) {
  if(! is.data.table(always_selected)) {
    stop('always_selected should be data.table.')
  }
  
  for_events = unique(rates$rates[, .(event_type, range_id)])
  if(always_selected[, .N] > 0) {
    stopifnot(all(c('range_id', 'event_type') %in% names(always_selected)))
    for_events = for_events[! always_selected, on = c('range_id', 'event_type')]
    events_to_test = split(always_selected$range_id, always_selected$event_type)
  } else {
    events_to_test = list()
  }
  events_to_cycle = split(for_events$range_id, for_events$event_type)
  
  output = lapply(names(events_to_cycle),
         function(curr_event_type) {
           message("Testing models with one more ", curr_event_type, ' effect term...')
           pblapply(events_to_cycle[[curr_event_type]],
                     function(curr_event) {
                       curr_to_test = copy(events_to_test)
                       curr_to_test[[curr_event_type]] = c(curr_to_test[[curr_event_type]], curr_event)
                       return(run_seg_multi(events_to_test = curr_to_test, rates = rates, debug = debug))
                     }, cl = cores)
         })
  output = unlist(output, recursive = FALSE, use.names = FALSE)
  
  si_out = rbindlist(lapply(output,  '[[', 'si'), idcol = 'run')
  ll = sapply(output, function(x) x$info$loglik)
  num_samples = sapply(output, function(x) length(x$info$included_samples))
  si_out[, loglik := ll[run]] # same order
  si_out[, num_samples := num_samples[run]]
  if(always_selected[, .N] > 0) {
    si_out[always_selected, is_fixed := TRUE, on = c('range_id', 'event_type')]
    si_out[is.na(is_fixed), is_fixed := FALSE]
  } else {
    si_out[, is_fixed := FALSE]
  }

  
  anno_cols = intersect(c('gene1', 'all_genes'), names(si_out))
  setcolorder(si_out, c('run', 'range_id', 'event_type', 'si', anno_cols))
  
  return(si_out)
}

# For visualizations
chr_coordinates = copy(prepped_data$prepped_calls$effective_coordinates)
chr_coordinates[, chr_label := paste0('chr', chr)]
plot_round = function(round_output) {
  curr_chr = unique(round_output$chr)
  stopifnot(length(curr_chr) == 1)
  if(! curr_chr %like% '^chr') {
    curr_chr = paste0('chr', curr_chr)
  }
  chr_info = chr_coordinates[curr_chr, on = 'chr_label'][1]
  fp = round_output[, .(mp = floor((start + end)/2)[1], gene = gene1[1], si = median(si), 
                        is_fixed = is_fixed[1], cancer_anno = gene1_anno[1]),
                    by = c('range_id', 'event_type')]
  fp[cancer_anno %in% c('oncogene', 'TSG'), gene_label := c(gene, rep(NA, .N - 1)), 
     by = c('gene', 'event_type')]
  palette =  c(RColorBrewer::brewer.pal(3, "Spectral"))
  
  fp[, inner_fill := fcase(cancer_anno == 'oncogene', palette[1],
                           cancer_anno == 'other', palette[2],
                           cancer_anno == 'TSG', palette[3],
                           default = 'gray')]
  
  fp[is.na(cancer_anno) | cancer_anno == 'noncancer', cancer_anno := 'noncancer or\nintergenic']
  palette[4] = 'gray'
  fill_legend = unique(fp[, .(cancer_anno, inner_fill)])
  fill_legend[, legend_order := match(inner_fill, palette)]
  fill_legend = fill_legend[order(legend_order)]
  
  for_fixed_labels = fp[is_fixed == T]
  fp[is_fixed == T, gene_label := NA]
  for_fixed_labels[is.na(gene_label), gene_label := gene] # even plotting unimportant genes

  gg = ggplot(fp, aes(x = mp, y = si)) + 
    geom_rect(data = chr_info, aes(xmin = cen_start, xmax = cen_end, 
                  ymin = -Inf, ymax = Inf), fill = 'gray95',
              inherit.aes = FALSE) + 
    geom_point(aes(color = inner_fill), size = 2.5) +
    xlab(paste0(curr_chr, ' position')) + ylab('Scaled selection for copy change') + 
    scale_color_identity(name = 'Cancer annotation', breaks = fill_legend$inner_fill,
                        labels = fill_legend$cancer_anno,
                        guide = guide_legend()) + 
    facet_wrap(~chr_label) +
    geom_text_repel(na.rm = T, aes(label = gene_label), size = 2, 
                    min.segment.length = 0, point.padding = .5, nudge_y = .07, segment.colour = 'gray20') +
    geom_label_repel(data = for_fixed_labels, na.rm = T, aes(label = gene_label), size = 2.5,
                     min.segment.length = 0, point.padding = .5, nudge_y = .07, segment.colour = 'gray20',
                     color = 'white',
                     fill = 'darkmagenta', box.padding = .15) +
    scale_x_continuous(labels = scales::label_comma()) +
    theme_light() + facet_wrap(~event_type, ncol = 1, scales = 'free')
  return(gg)
}


# Run inference over all intervals (decreases and increases).
# Adjust cores as necessary for your system.
r1 = run_round(rates = chr7_rates, cores = 2, debug = T)
p1 = plot_round(r1)
ggsave(file.path(output_dir, 'chr7_r1.png'), p1,
       width = 8, height = 6)

# EGFR decrease (negatively selected) is the third-best by likelihood, and we'll choose it.
stopifnot(r1[order(-loglik), gene1[3] == 'EGFR'])
to_fix_in_r2 = r1[order(-loglik)][gene1 == 'EGFR', .(range_id, event_type)][1]

r2 = run_round(always_selected = to_fix_in_r2, rates = chr7_rates, cores = 2)

p2 = plot_round(r2)
ggsave(file.path(output_dir, 'chr7_two_rounds.png'), p2,
       width = 8, height = 6)

to_fix_in_r3 = rbind(to_fix_in_r2,
                     r2[order(-loglik)][is_fixed == FALSE, .(range_id, event_type)][1])

r3 = run_round(always_selected = to_fix_in_r3, rates = chr7_rates, cores = 2)

# Examine other genes in top runs
# Issue: By likelihood, the best models use sites right next to EGFR. My guess is that this way,
# the model doesn't have to explain as much.
r2[order(-loglik)][is_fixed == F]

saveRDS(list(r1 = r1, r2 = r2, chr7_rates = chr7_rates),
        file.path(output_dir, 'curr_chr7_output.rds'))

## Combine all plots (to be continued...)
all_plots = list(p1, p2)
all_plots = lapply(all_plots,
                   function(x) x + xlab('') + ylab('Effect') +
                     theme(strip.background = element_blank(),
                           strip.text = element_blank()))
grid = cowplot::plot_grid(plotlist = all_plots, ncol = 1,
                   labels = paste0('Round ', 1:4), label_size = 10,
                   label_x = 0, vjust = 0.3)
# # Add empty plot for sufficient whitespace at the top
# grid = cowplot::plot_grid(ggplot() + theme_void(), grid, rel_heights = c(.04, .96), nrow = 2)



# tmp = run_seg_multi(events_to_test = list(increase = c(to_test$range_id, 'chr7.10'),
#                                     decrease = to_test$range_id), rates = chr7_rates)
#               
# chr13_rates = get_cna_rates(chr = '13', cna_calls = cna_calls, seg_rates = seg_rates, rate_intervals = subset)
# chr13_r1 = run_round(rates = chr13_rates)
# 
# # Top result is selection against RB1 increase
# chr13_fix_r2 = data.table(range_id = 'chr13.33', event_type = 'increase')
# chr13_r2 = run_round(chr13_fix_r2, rates = chr13_rates, cores = 2)

