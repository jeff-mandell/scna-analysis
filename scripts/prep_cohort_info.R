library(data.table)

clinical = fread('data/clinical.cohort.2024-09-11/clinical.tsv.gz', na.strings = "'--")
cna_calls = fread('prepped_data/TCGA_ASCAT3_hg38.txt.gz')

# Load in sex (needed for expected chrX copy)
clinical = unique(clinical[case_submitter_id %in% cna_calls$patient_id, 
                           .(patient_id = case_submitter_id, project_id,
                             sex = fcase(gender == 'female', 'F', gender == 'male', 'M', default = NA))])
clinical = clinical[patient_id %in% cna_calls$patient_id]
fwrite(clinical, 'prepped_data/TCGA_patient_info.txt', sep = "\t", na = 'NA')

