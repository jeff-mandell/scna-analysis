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

# We will calculate local mutation rates over megabase-sized windows
megabase_ranges = get_genomic_windows(window_size = 1e6, refset = 'ces.refset.hg38',
                                      ranges = prepped_data$prepped_calls$effective_coordinates)



# This local rates function currently uses particular information about TCGA projects.
# In the future, the function will likely use different method.
# The CNA calls must also be annotated with project.
calls[patient_info, project := project, on = 'sample']
patient_info[prepped_data$prepped_calls$ploidy_calls, is_nondiploid := simple_ploidy > 2, on = 'sample']

project_sample_counts = patient_info[, .(.N, N_nondiploid = sum(is_nondiploid)), by = 'project']

r_chr10 = get_local_rates(calls = calls, chr_ranges = megabase_ranges, chr = '10',
                              project_sample_counts = project_sample_counts,
                              gene_coord = ces.refset.hg38$cancer_gene_coord)

# Temporary issue: get some decrease rates of 0. Setting them to lowest nonzero rates and updating r0.
lowest_nonzero = r_chr10[r > 0 & rate_type == 'decrease', min(r)]
r_chr10[rate_type == 'decrease' & r == 0, let(r = lowest_nonzero, r0 = 1 - lowest_nonzero)]

project_by_sample = patient_info[, .(sample, project)]

# Take cancer-annotated genes and every fifth gene regardless.
chr10_test_genes =  ces.refset.hg38$cancer_gene_coord[chr == '10'][, .SD[cancer_anno != 'noncancer' | .I %% 5 == 0]][, gene]
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

