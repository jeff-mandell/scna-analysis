library(data.table)
library(cancereffectsizeR)

# Run from directory where you wish to save the files.
projects = fread('reference/all_TCGA_projects.txt', header = F)[[1]]
projects = setdiff(projects, 'TCGA-SKCM')
for (proj in projects) {
  output_name = paste0(proj, '.maf.gz')
  if(! file.exists(output_name)) {
    get_TCGA_project_MAF(project = proj, filename = output_name)
  }
}

# Special handling for TCGA-SKCM
skcm_maf = 'TCGA-SKCM_with_met_exclude_multisample.maf.gz'
if (! file.exists(skcm_maf)) {
  tmp_maf = paste0(tempfile(), '.maf')
  get_TCGA_project_MAF(project = 'TCGA-SKCM', filename = tmp_maf, exclude_TCGA_nonprimary = FALSE)
  
  # 2 patients have 2 samples each. For simplicity, we'll remove these patients.
  maf_to_edit = fread(tmp_maf)
  stopifnot(maf_to_edit[multisample_patient == T, 
                        uniqueN(Tumor_Sample_Barcode) == 4 && uniqueN(patient_id) == 2])
  maf_to_edit = maf_to_edit[multisample_patient == FALSE]
  fwrite(maf_to_edit, skcm_maf, sep = "\t")
  unlink(tmp_maf)
}

