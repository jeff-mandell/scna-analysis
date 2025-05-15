library(data.table)
filters = sprintf(
  '{
    "op": "and",
    "content": [
      {
        "op": "in",
        "content":{
          "field": "access",
          "value": ["open"]
        }
      },
      {
        "op":"=",
        "content":{
          "field":"files.data_type",
          "value": ["Allele-specific Copy Number Segment"]
        }
      },
      {
        "op": "in",
        "content":{
          "field": "analysis.workflow_type",
          "value": ["ASCAT3"]
        }
      }
    ]
  }')
num_files = format(1e5, scientific = F) # file limit passed to query must be given as character
files_endpt = 'https://api.gdc.cancer.gov/files'
response = httr::GET(files_endpt, 
                     query = list(filters = filters, 
                                  fields = "file_name,md5sum,cases.samples.submitter_id", size = num_files, format = 'JSON'))

stopifnot(response$status_code == 200)
content = rjson::fromJSON(rawToChar(response$content))
files = rbindlist(lapply(content$data$hits, '[', c("id", "file_name", "md5sum")))
files[, patient_id := sapply(content$data$hits, function(x) substr(x$cases[[1]]$samples[[1]]$submitter_id, 1, 12))]
stopifnot(uniqueN(files$patient_id) == files[, .N])


stopifnot(all(files$file_name %like% 'ascat3'))

dest_dir = 'data/TCGA_ASCAT'
data_endpt = 'https://api.gdc.cancer.gov/data'

files[, path := paste0(dest_dir, '/', file_name)]
files[, url := paste0(data_endpt, '/', id)]

pbapply::pbmapply(
  function(url, path, filename) {
    # Try up to 5 times to download a file
    num_tries = 5
    
    for(i in 1:(num_tries - 1)) {
      unlink(path) # in case file exists from previous attempt
      code = 1L
      tryCatch(
        { code = utils::download.file(url = url, destfile = path, quiet = TRUE, mode = 'wb')},
        error = function(e) NULL, warning = function(w) NULL
      )
      if(is.integer(code) && code == 0) return()
      Sys.sleep(1) # wait a second
    }
    
    # On last try, we'll let warnings/errors go through.
    message("\nHaving some trouble with a download. Waiting 10 seconds and giving it one last try....\n")
    Sys.sleep(10)
    unlink(path)
    withCallingHandlers(
      {
        code = utils::download.file(url = url, destfile = path, quiet = TRUE, mode = 'wb')
      }, error = function(e) {
        msg = paste0(strwrap(paste0("File ", filename, " failed to download from ",
                                    url, " (tried ", num_tries, " times).")), collapse = "\n")
        warning("\n", msg, "\nLast errors/warnings:\n", call. = FALSE, immediate. = TRUE)
        e
      }, warning = function(w) w
    )
  }, files$url, files$path, files$file_name)


stopifnot(all(file.exists(files$path)))

actual_sums = tools::md5sum(files$path)
files[names(actual_sums), obs_sum := actual_sums, on = 'path']
files[, failed := obs_sum != md5sum]
stopifnot(files[, sum(failed) == 0])

to_read = setNames(files$path, files$patient_id)
combined = rbindlist(lapply(to_read, fread), idcol = 'patient_id')

stopifnot(uniqueN(combined$GDC_Aliquot) == uniqueN(combined[, .(GDC_Aliquot, patient_id)]))

# Keep track of what GDC file each record comes from
combined[files, file_id := id, on = 'patient_id']
fwrite(combined, 'prepped_data/TCGA_ASCAT3_hg38.txt.gz', sep = "\t")

