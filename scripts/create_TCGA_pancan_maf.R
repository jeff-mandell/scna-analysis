library(data.table)
library(cancereffectsizeR)

# Untracked symbolic link to directory of TCGA MAF files.
# See get_all_TCGA_maf.R to create the MAF file directory.
tcga_files = list.files('data/TCGA_link/', pattern = 'maf\\.gz$', full.names = T, recursive = F)
stopifnot(length(tcga_files) == 33) # Should be 33 projects
names(tcga_files) = sub('TCGA-', '', sub('\\.maf.*', '', basename(tcga_files)))
big_maf = rbindlist(lapply(tcga_files, fread), fill = T, idcol = 'project')
big_maf[project == 'SKCM_with_met_exclude_multisample', project := 'SKCM']
big_maf = preload_maf(maf = big_maf, refset = 'ces.refset.hg38', keep_extra_columns = T)

# Let's remove problems and apply usual filtering
big_maf = big_maf[is.na(problem)]
big_maf = big_maf[germline_variant_site == F][repetitive_region == FALSE | cosmic_site_tier %in% 1:3]
big_maf[HGVSp_Short == '', HGVSp_Short := NA]
smaller_maf = big_maf[, .(patient_id = Unique_Patient_Identifier, Chromosome, Start_Position, Reference_Allele, Tumor_Allele,
     variant_id, variant_type, project, gene = Hugo_Symbol, t_ref_count, t_alt_count, aachange = sub('^p\\.', '', HGVSp_Short), pid = ENSP, 
     SIFT, PolyPhen, One_Consequence)]
smaller_maf[! is.na(One_Consequence), is_synonymous := One_Consequence == 'synonymous_variant']
smaller_maf[! is.na(One_Consequence), is_lof := One_Consequence %in% c('start_loss', 'stop_loss', 'stop_gained', 'frameshift_variant')]
smaller_maf[! is.na(One_Consequence), is_splice_variant := One_Consequence %like% 'splice']

# Fix project annotation (missing on DBS/other variants)
project_to_sample = unique(smaller_maf[! is.na(project), .(patient_id, project)])
smaller_maf[project_to_sample, project := i.project, on = 'patient_id']
smaller_maf[, is_missense := One_Consequence == 'missense_variant']

# Note that annotations are missing for DBS/other variants. Adding canonical gene annotations for these.
cds_intervals = ces.refset.hg38$transcripts[is_mane == T & type == 'CDS', .(chr, start, end, gene = gene_name)]
setkey(cds_intervals, 'chr', 'start', 'end')

# Will keep whatever gene annotation came with SNVs.
to_update = smaller_maf[is.na(gene) & variant_type != 'snv']
to_update[, let(start = Start_Position, end = Start_Position)]
setkey(to_update, 'Chromosome', 'start', 'end')
gene_overlaps = foverlaps(to_update[, .(Chromosome, Start_Position, Reference_Allele, Tumor_Allele, patient_id,
                                        start, end)], cds_intervals, type = 'within', nomatch = NULL)
gene_overlaps[, c('i.start', 'i.end', 'start', 'end') := NULL]
multi_gene_overlaps = gene_overlaps[, .SD[.N > 1], by = setdiff(names(gene_overlaps), 'gene')]

# Just a few multiple-overlap cases (and it's 3 loci that together have 13 gene annotations)
stopifnot(multi_gene_overlaps[, .N] == 13)

# We'll arbitrarily take the first of these.
gene_overlaps = gene_overlaps[, .SD[1], by = setdiff(names(gene_overlaps), 'gene')]

prev = copy(smaller_maf)
smaller_maf[gene_overlaps, new_gene := i.gene, on = setdiff(names(gene_overlaps), 'gene')]
smaller_maf[is.na(gene), gene := new_gene]
smaller_maf[, new_gene := NULL]

fwrite(smaller_maf, 'prepped_data/panTCGA.maf.gz', sep = "\t")

