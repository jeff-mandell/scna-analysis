## Run from cancereffectsizeR dev directory
## You also need the latest development version of ces.refset.hg38 and the SCNA data repository  
library(devtools)
library(data.table)
library(pbapply)
library(ggplot2)
library(ggrepel)
setwd(usethis::proj_get())
load_all()
load_all('~/ces.refset.hg38') # UPDATE as necessary

# Load output from prep_SCNA_analysis.R (adjust path as needed)
output_dir = '~/Desktop/repos/scna-analysis/untracked/outputs/'
prepped_data = readRDS(file.path(output_dir, 'prepped_cna.rds'))
calls = copy(prepped_data$prepped_calls$calls)

# Data directory from SCNA repo (adjust as needed)
data_dir = '~/Desktop/repos/scna-analysis/prepped_data' 
patient_info = fread(file.path(data_dir, 'TCGA_patient_info.txt'))
setnames(patient_info, c('patient_id', 'project_id'),  c('sample', 'project'))
patient_info = patient_info[sample %in% calls$sample]

# We will calculate local mutation rates over 100Mb-sized windows
hectobase_ranges = get_genomic_windows(window_size = 1e5, refset = 'ces.refset.hg38',
                                      ranges = prepped_data$prepped_calls$effective_coordinates)



# This local rates function currently uses particular information about TCGA projects.
# In the future, the function will likely use different method.
# The CNA calls must also be annotated with project.
calls[patient_info, project := project, on = 'sample']
patient_info[prepped_data$prepped_calls$ploidy_calls, is_nondiploid := simple_ploidy > 2, on = 'sample']

project_sample_counts = patient_info[, .(.N, N_nondiploid = sum(is_nondiploid)), by = 'project']

r_chr10 = get_local_rates(calls = calls, chr_ranges = hectobase_ranges, chr = '10',
                              project_sample_counts = project_sample_counts,
                              gene_coord = ces.refset.hg38$cancer_gene_coord)

# Temporary issue: get some decrease rates of 0. Setting them to lowest nonzero rates and updating r0.
lowest_nonzero = r_chr10[r > 0 & rate_type == 'decrease', min(r)]
r_chr10[rate_type == 'decrease' & r == 0, let(r = lowest_nonzero, r0 = 1 - lowest_nonzero)]

project_by_sample = patient_info[, .(sample, project)]

# Take cancer-annotated genes and every fifth gene regardless.
chr10_test_genes =  ces.refset.hg38$cancer_gene_coord[chr == '10'][, .SD[! cancer_anno %in% c('TSG', 'oncogene') | .I %% 10 == 0]][, gene]
chr10_test_genes = chr10_test_genes[chr10_test_genes %in% r_chr10$gene]
gene_coord = ces.refset.hg38$cancer_gene_coord[chr10_test_genes, on = 'gene']

# This would be significantly faster if we went by window instead of gene and then looked up gene by window.
chr10_prob_model = pblapply(chr10_test_genes, prep_gene_probability_model, 
                           event_type = 'decrease', prepped_data = prepped_data, 
                           rates = r_chr10[rate_type == 'decrease'], gene_coord = gene_coord,
                           project_by_sample = project_by_sample, cl = 1)
names(chr10_prob_model) = chr10_test_genes


# Adjust cores as necessary
si_fits = pblapply(chr10_prob_model, 
                 function(x) {
                   output = 'begin'
                   output = tryCatch(
                     expr = {
                       return(coef(bbmle::mle2(x$lik, start = c(gamma = 1), 
                                               method = 'L-BFGS-B', lower = 1e-6, upper = 1e6)))
                     }, error = function(e) {
                       return('failed')
                     })
                   return(output)
                 }, cl = 4
)
names(si_fits) = chr10_test_genes
chr10_si = data.table(gene = names(si_fits), si  = unlist(si_fits))

chr10_si[gene_coord, let(chr = chr, midpoint = floor((start + end)/2), gene_type = cancer_anno), on = 'gene']
chr10_si[, chr_label := paste0('chr', chr)]
chr10_si[gene  == 'PTEN', gene_label := gene]
chr10_si[, gene_type := factor(gene_type, levels = c('oncogene', 'other', 'TSG', 'noncancer'))]
palette =  c(RColorBrewer::brewer.pal(3, "Spectral"), 'grey')
chr10_si = chr10_si[order(-gene_type)] # so that the important genes don't get overplotted by noncancer

set.seed(910)
gg_prob_round1 = ggplot(chr10_si, aes(x = midpoint, y = si)) + geom_point(aes(color = gene_type), size = 1) +
  xlab('Position') + ylab('Scaled selection for copy decrease') + scale_color_manual(name = 'Cancer annotation', values = palette) +
  geom_text_repel(na.rm = T, aes(label = gene_label), size = 2, nudge_y = .15, segment.colour = 'gray20') +
  theme_light() + 
  facet_wrap(~chr_label)
ggsave(file.path(output_dir, 'chr10_prob_r1.png'), gg_prob_round1, width = 8, height = 4)

## Demo of running round 2 of inference. It needs some fixes, so this is just illustrative.
# Run round 2, using our wisdom to select PTEN as the fixed effect.
r1_info = c(chr10_prob_model$PTEN, list(si = si_fits$PTEN))
other_chr10_genes = setdiff(chr10_test_genes, 'PTEN')
chr10_run2 = pblapply(other_chr10_genes, prep_gene_probability_model, 
                            event_type = 'decrease', prepped_data = prepped_data, 
                            rates = r_chr10[rate_type == 'decrease'], gene_coord = gene_coord,
                            project_by_sample = project_by_sample, r1_info = r1_info, cl = 1)
names(chr10_run2) = other_chr10_genes

r2_fits = pblapply(chr10_run2, 
                   function(x) {
                     output = 'begin'
                     output = tryCatch(
                       expr = {
                         return(coef(bbmle::mle2(x$lik, start = c(gamma = 1), 
                                                 method = 'L-BFGS-B', lower = 1e-6, upper = 1e6)))
                       }, error = function(e) {
                         return('failed')
                       })
                     return(output)
                   }, cl = 4
)

# Some inferences failed. For now, just exclude.
r2_si = data.table(gene = names(r2_fits), si_tmp  = unlist(r2_fits))
r2_si = r2_si[si_tmp != 'failed']
r2_si[, si := as.numeric(si_tmp)]
r2_si[gene_coord, let(chr = chr, midpoint = floor((start + end)/2), gene_type = cancer_anno), on = 'gene']
r2_si[, chr_label := paste0('chr', chr)]
r2_si[, gene_type := factor(gene_type, levels = c('oncogene', 'other', 'TSG', 'noncancer'))]
palette =  c(RColorBrewer::brewer.pal(3, "Spectral"), 'grey')
r2_si = r2_si[order(-gene_type)] # so that the important genes don't get overplotted by noncancer


r2_plot = ggplot(r2_si, aes(x = midpoint, y = si)) + geom_point(aes(color = gene_type), size = 1) +
  xlab('Position') + ylab('Scaled selection for copy decrease') + scale_color_manual(name = 'Cancer annotation', values = palette) +
  theme_light() + 
  facet_wrap(~chr_label) + geom_hline(yintercept = 1)




# Alternative model, based on segmentation rates.
# Currently, this model makes no use of regional mutation rates.
seg_rates = get_seg_increase_rates(cna_calls = prepped_data$prepped_calls$calls,
                                   increase_burden = prepped_data$cna_burdens$increase_burden)
# chr7 has EGFR, BRAF
chr7_test_genes = ces.refset.hg38$cancer_gene_coord[chr == '7', gene]
seg_prep = pblapply(chr7_test_genes, prep_seg_lik, event_type = 'increase', cna_calls = prepped_data$prepped_calls$calls,
                    seg_rates = seg_rates)
names(seg_prep) = chr7_test_genes
seg_fit = lapply(seg_prep,
                 function(x) {
                   output = 'begin'
                   output = tryCatch(
                     expr = {
                       return(coef(bbmle::mle2(x$lik, start = c(gamma = 1), vecpar = T, 
                                               lower = 1e-6, upper = 1e9, method = 'L-BFGS-B')))
                     }, error = function(e) {
                       return('failed')
                     })
                   return(output)
                 }
)
seg_si = data.table(gene = names(seg_fit), si = unlist(seg_fit))
seg_si[ces.refset.hg38$cancer_gene_coord, let(chr = chr, midpoint = floor((start + end)/2), gene_type = cancer_anno), on = 'gene']

seg_si[, chr_label := paste0('chr', chr)]
seg_si[gene %in% c('BRAF', 'EGFR', 'KRAS', 'MDM2', 'CDK6'), gene_label := gene]
seg_si[, gene_type := factor(gene_type, levels = c('oncogene', 'other', 'TSG', 'noncancer'))]
palette =  c(RColorBrewer::brewer.pal(3, "Spectral"), 'grey')

# so that the important genes don't get overplotted by noncancer
# Plot of segment rate model effects (chr7)
seg_gg = ggplot(seg_si, aes(x = midpoint, y = si)) + geom_point(aes(color = gene_type), size = 1) +
  xlab('Position') + ylab('Scaled selection for copy increase') + scale_color_manual(name = 'Cancer annotation', values = palette) + 
  facet_wrap(~chr_label) +
  geom_text_repel(na.rm = T, aes(label = gene_label), size = 2, nudge_y = .15, segment.colour = 'gray20') +
  theme_light()

ggsave(file.path(output_dir, 'chr7_cna_increase_segmentation_method.png'), seg_gg, width = 8, height = 4)



# Let's try IKZF1, a TSG at ~7:5e7
gene_name = 'IKZF1'
cna_calls = copy(prepped_data$prepped_calls$calls)
event_type = 'increase'

seg_rates = seg_rates # already calculated above
gene_coord = copy(ces.refset.hg38$cancer_gene_coord)
r1_info = list(interval = seg_prep$EGFR$interval, si = seg_fit$EGFR[[1]])

prep_seg_lik_r2 = function(gene_name, cna_calls = NULL, event_type = 'increase', seg_rates,
                        gene_coord = ces.refset.hg38$cancer_gene_coord, r1_info = NULL) {
  total_rates = copy(seg_rates$total_rates)
  seg_rates = copy(seg_rates$seg_rates)
  setkey(cna_calls, chr, start, end)
  
  # This will be moved into CES internal data later.
  sr_info = data.table(pretty_label = c('L < 100 kb', '100 kb < L < 1 Mb', '1 Mb < L < 10 Mb', '10 Mb < L < 40 Mb', '40 Mb < L'),
                       cosmic_label = c('0-100kb', '100kb-1Mb', '1Mb-10Mb', '10Mb-40Mb', '>40Mb'),
                       bin_id = c('r1', 'r2', 'r3', 'r4', 'r5'),
                       bin_order = 1:5)
  
  groups = sr_info$cosmic_label # size ranges
  
  stopifnot(all(c('sample', 'size_range', 'seg_rate') %in% names(seg_rates)))
  curr_rates = seg_rates[size_range %in% groups]
  
  if(curr_rates[, .N] == 0) {
    warning('No rates found; returning empty table')
    return(NULL)
  }
  
  gr = gene_coord[gene == gene_name, .(chr, start, end)]
  setkey(gr, chr, start, end)
  curr_ol = foverlaps(cna_calls, gr, nomatch = NULL)
  
  curr_samples = unique(curr_ol$sample)
  curr_n = length(curr_samples)
  
  curr_ol[curr_rates, neutral_rate := seg_rate, on = c('sample', 'size_range')]
  
  #neutral_prob = curr_rates[, .(exp = (rate/n_nondiploid) * curr_n), by = 'size_range']
  
  
  # If more than one event in a sample, we will count the largest event
  if(event_type == 'decrease') {
    curr_ol = curr_ol[total_copy < exp_total_copy]
  } else if(event_type == 'loh') {
    curr_ol = curr_ol[nMinor == 0] # Yeah, total_copy might increase! Yet something was lost.
  } else if(event_type == 'increase') {
    curr_ol = curr_ol[total_copy > exp_total_copy]
  }
  
  outer_ol = copy(curr_ol)
  
  if(! is.null(r1_info)) {
    
    curr_interval = gene_coord[gene == gene_name, .(chr, int_start = start, int_end = end)]
    setkey(curr_interval, chr, int_start, int_end)
    other_interval = r1_info$interval
    si_other = r1_info$si
    
    if(curr_interval$int_start > other_interval$end) {
      dist_to_other = curr_interval$int_start - other_interval$end
      rel_pos = 'right'
    } else if(curr_interval$int_end < other_interval$start) {
      dist_to_other = other_interval$start - curr_interval$int_end
      rel_pos = 'left'
    } else {
      stop('Gene intervals overlap.')
    }
    
    p_by_size = list()
    # order doesn't matter here
    for(curr_bin_id in sr_info$bin_id) {
      curr_size_range = sr_info[bin_id == curr_bin_id, cosmic_label]
      curr_rel_neut = cna_calls[size_range == curr_size_range]
      setkey(curr_rel_neut, 'chr', 'start', 'end')
      curr_ol = foverlaps(curr_rel_neut, curr_interval, nomatch = NULL)
      ol_both = foverlaps(curr_ol, other_interval, nomatch = NULL)
      if(ol_both[, .N] == 0) {
        p_by_size[[curr_bin_id]] = 0
      } else {
        p_by_size[[curr_bin_id]] = ol_both[, .N]/curr_ol[, .N]
      }
    }
    
    p_by_size = setNames(as.numeric(p_by_size), names(p_by_size))
    curr_rates[sr_info, bin_id := bin_id, on = c(size_range = 'cosmic_label')]
    curr_rates[, curr_p := p_by_size[bin_id]]
    curr_rates[, seg_rate := seg_rate * (1 - seg_rate + si_other * seg_rate)]
    
    min_by_type = curr_rates[seg_rate != 0, .(min(seg_rate)), by = 'bin_id']
    curr_rates[seg_rate == 0, seg_rate := min_by_type[bin_id, V1, on = 'bin_id']]
    rm(curr_rel_neut)
  }
  curr_ol = copy(outer_ol)
  
  neutral_prob = unique(curr_ol[, .(sample, size_range, neutral_rate)]) # sample might have two same-size events
  neutral_prob[sr_info, bin_order := bin_order, on = c(size_range = 'cosmic_label')]
  
  # Get rid of smaller size entries where samples have multiple entries 
  neutral_prob = unique(neutral_prob[order(-bin_order)], by = 'sample')
  
  samples_with = list()
  
  # We go from largest to smallest bin
  for(i in sort(sr_info$bin_order, decreasing = TRUE)) {
    curr_bin = sr_info[bin_order == i, bin_id]
    curr_size_range = sr_info[bin_order == i, cosmic_label]
    
    curr_samples_with = curr_ol[size_range == curr_size_range, unique(sample)]
    
    # Larger events in the same samples take precedence (already handled above, really)
    samples_with[[curr_bin]] = setdiff(curr_samples_with, unlist(samples_with))
  }
  samples_with$with_any = unlist(samples_with)
  samples_with$with_none = setdiff(curr_samples, samples_with$with_any)
  
  
  # for lik
  combined_rate = setNames(total_rates$seg_rate, total_rates$sample)[curr_samples]
  with_none = samples_with$with_none
  neutral_prob[, total_rate := combined_rate[sample]]
  cond_prob = setNames(neutral_prob[, neutral_rate/total_rate], neutral_prob$sample)
  
  lik = function(gamma) {
    sum_log_lik = 0
    if (length(with_none) > 0) {
      sum_log_lik = -1 * sum(gamma * combined_rate[with_none])
    }
    
    for(bin_id in sr_info$bin_id) {
      curr_samples = samples_with[[bin_id]]
      curr_combined_rate = combined_rate[curr_samples] * gamma
      prob_any = 1 - exp(-curr_combined_rate)
      prob_curr = prob_any * cond_prob[curr_samples]
      sum_log_lik = sum_log_lik + sum(log(prob_curr))
    }
    return(-1 * sum_log_lik)
  }
  formals(lik)[["gamma"]] = 1
  bbmle::parnames(lik) = "gamma"
  gene_interval = gr[, .(chr = chr[1], start = min(start), end = max(end))]
  output = list(lik = lik, 
                groups = samples_with,
                num_covering = length(curr_samples), n_neutral_dec = sum(curr_rates$rate),
                interval = gene_interval,
                gene_name = gene_name)
  return(output)
}

chr7_test_genes = setdiff(ces.refset.hg38$cancer_gene_coord[chr == '7', gene], 'EGFR')
seg_prep = pblapply(chr7_test_genes, prep_seg_lik_r2, event_type = 'increase', 
                    cna_calls = prepped_data$prepped_calls$calls,
                    seg_rates = seg_rates, r1_info = r1_info)
names(seg_prep) = chr7_test_genes
seg_fit = lapply(seg_prep,
                 function(x) {
                   output = 'begin'
                   output = tryCatch(
                     expr = {
                       return(coef(bbmle::mle2(x$lik, start = c(gamma = 1), vecpar = T, 
                                               lower = 1e-6, upper = 1e9, method = 'L-BFGS-B')))
                     }, error = function(e) {
                       return('failed')
                     })
                   return(output)
                 }
)


seg_si = data.table(gene = names(seg_fit), si = unlist(seg_fit))
seg_si[ces.refset.hg38$cancer_gene_coord, let(chr = chr, midpoint = floor((start + end)/2), gene_type = cancer_anno), on = 'gene']

seg_si[, chr_label := paste0('chr', chr)]
seg_si[gene %in% c('BRAF', 'EGFR', 'KRAS', 'MDM2', 'CDK6'), gene_label := gene]
seg_si[, gene_type := factor(gene_type, levels = c('oncogene', 'other', 'TSG', 'noncancer'))]
palette =  c(RColorBrewer::brewer.pal(3, "Spectral"), 'grey')


gg2 = ggplot(seg_si, aes(x = midpoint, y = si)) + geom_point(aes(color = gene_type), size = 1) +
  xlab('Position') + ylab('Scaled selection for copy increase') + scale_color_manual(name = 'Cancer annotation', values = palette) + 
  facet_wrap(~chr_label) +
  geom_text_repel(na.rm = T, aes(label = gene_label), size = 2, nudge_y = .15, segment.colour = 'gray20') +
  theme_light()
