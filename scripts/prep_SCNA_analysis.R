## Run from cancereffectsizeR dev directory
## You also need the latest development version of ces.refset.hg38 and the SCNA data repository  
library(devtools)
library(data.table)
setwd(usethis::proj_get())
load_all()
load_all('~/ces.refset.hg38') # Change path as necessary

# The version of ces.refset.hg38 that you need resides on the dev branch.
# Install with remotes::install_github("Townsend-Lab-Yale/ces.refset.hg38@dev")
stopifnot(packageVersion('ces.refset.hg38') >= as.package_version('1.3.0.9000'))

# cancereffectsizeR v3 also from dev branch. Rather than installing with install_github,
# clone the repo and use devtools as shown above.
stopifnot(packageVersion('cancereffectsizeR') >= as.package_version('3.0.0.9000'))

# You also need a custom installation of BISCUT.
if(! require('BISCUT') || packageVersion('BISCUT') < as.package_version('1.2')) {
  remotes::install_github('jeff-mandell/BISCUT-py3')
}

# Data folder within scna-analysis (available as another Github repo).
# And a directory for output files (can be any directory).
data_dir = '~/Desktop/repos/scna-analysis/prepped_data' 
output_dir = '~/Desktop/repos/scna-analysis/untracked/outputs'
dir.create(output_dir, showWarnings = FALSE)


cna_calls = fread(file.path(data_dir, 'TCGA_ASCAT3_hg38.txt.gz'))

# Load in sex (needed for expected chrX copy)
patient_info = fread(file.path(data_dir, 'TCGA_patient_info.txt'))

cna_calls[patient_info, sex := sex, on = 'patient_id']

# 27 samples lack sex information. For simplicity, drop them.
stopifnot(cna_calls[is.na(sex), uniqueN(patient_id) == 27])
cna_calls = cna_calls[! is.na(sex)]

# Filter to just patients in use (>99%)
patient_info = patient_info[patient_id %in% cna_calls$patient_id]
prepped_calls = prep_ASCAT3_segments(segments = cna_calls, refset = 'ces.refset.hg38')

# Adjust cores as necessary for your system.
prepped_calls = call_large_events(prepped_calls = prepped_calls, arm_chr_threshold = .99, 
                                 cores = 4, account_biscut_regions = TRUE, 
                                 biscut_dir = file.path(output_dir, 'biscut-out'),
                                 run_biscut = FALSE)



# Run MutationalPatterns (takes a bit).
# Note that CN signatures were derived from hg38 data. (Probably not a problem.)
cn_sig_def_file = system.file('extdata/COSMIC_v3.4_CN_GRCh37.txt', package ='ces.refset.hg38')

# To-do: should set seed here for reproducibility.
signature_output = cn_signature_extraction(sig_def = cn_sig_def_file, 
                                           cna_segments = prepped_calls$calls)


# For every sample and segment size class, get expected burden (under CN signature reconstruction)
# and rel_burden (proportion of overall burden among all samples within the size class).
# Output is 2-length list: separate data.table of burdens for increase_burden and decrease_burden.
cna_burdens = cna_class_relative_rates(cna_calls = prepped_calls$calls,
                                       cna_recon = signature_output$reconstruction,
                                       ploidy_calls = prepped_calls$ploidy_calls)

# Save stuff for next time.
# The unneeded plots that get included in MP output are cumulatively huge, so we'll remove.
signature_output$mp_out$sim_decay_fig = 'deleted'
to_save = list(prepped_calls = prepped_calls, signature_output = signature_output,
               cna_burdens = cna_burdens)
saveRDS(to_save, file.path(output_dir, 'prepped_cna.rds'))









