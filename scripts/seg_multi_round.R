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
chr7_rates = get_cna_rates(chr = '7', cna_calls = cna_calls, seg_rates = seg_rates, rate_intervals = gc)



## run_round will run a collection of inferences across the input chromosomal intervals
# fixed_events: data.table(); see example below.
# rates: from get_cna_rates()
run_round = function(fixed_events = data.table(), rates = NULL, cores = 1) {
  if(! is.data.table(fixed_events)) {
    stop('fixed_events should be data.table.')
  }
  
  for_events = unique(rates$rates[, .(event_type, range_id)])
  if(fixed_events[, .N] > 0) {
    stopifnot(all(c('range_id', 'event_type') %in% names(fixed_events)))
    for_events = for_events[! fixed_events$range_id, on = 'range_id']
    events_to_test = split(fixed_events$range_id, fixed_events$event_type)
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
                       return(run_seg_multi(events_to_test = curr_to_test, rates = rates))
                     }, cl = cores)
         })
  output = unlist(output, recursive = FALSE, use.names = FALSE)
  
  si_out = rbindlist(lapply(output,  '[[', 'si'), idcol = 'run')
  ll = sapply(output, function(x) x$info$loglik)
  num_samples = sapply(output, function(x) length(x$info$included_samples))
  si_out[, loglik := ll[run]] # same order
  si_out[, num_samples := num_samples[run]]
  si_out[, is_fixed := range_id %in% fixed_events$range_id]
  
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
  

  # fp[, border_color := fcase(is_fixed == TRUE, 'darkmagenta',
  #    default = inner_fill)]
  
  for_fixed_labels = fp[is_fixed == T]
  fp[is_fixed == T, gene_label := NA]
  
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
r1 = run_round(rates = chr7_rates, cores = 4)
p1 = plot_round(r1)

# EGFR decrease (negatively selected) is the third-best by likelihood, and we'll choose it.
to_fix_in_r2 = r1[order(-loglik)][gene1 == 'EGFR', .(range_id, event_type)][1]
r2 = run_round(fixed_events = to_fix_in_r2, rates = chr7_rates, cores = 4)

# Keep in mind, inferences with both increase and decrease effects are currently wrong.
p2 = plot_round(r2)

## Combine all plots (to be continued...)
# all_plots = list(p1, p2)
# all_plots = lapply(all_plots,
#                    function(x) x + xlab('') + ylab('Effect') +
#                      theme(strip.background = element_blank(),
#                            strip.text = element_blank()))
# grid = cowplot::plot_grid(plotlist = all_plots, ncol = 1,
#                    labels = paste0('Round ', 1:4), label_size = 10,
#                    label_x = 0, vjust = 0.3)
# # Add empty plot for sufficient whitespace at the top
# grid = cowplot::plot_grid(ggplot() + theme_void(), grid, rel_heights = c(.04, .96), nrow = 2)



