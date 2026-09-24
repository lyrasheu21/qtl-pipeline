# Nominal associations of the FDR-significant phenotypes of one set, written
# with named columns. The chunk files are filtered with awk on column 1 (the
# phenotype, or the group for --grp-best runs) before R reads them.
source(file.path(snakemake@scriptdir, "load_qtlpost.R"))

keys <- snakemake@input[["fdr_keys"]]
raw <- snakemake@output[["raw"]]
if (file.size(keys) > 0) {
  chunks <- paste(shQuote(unlist(snakemake@input[["chunks"]])), collapse = " ")
  cmd <- sprintf("for f in %s; do gzip -dc \"$f\"; done | awk 'NR == FNR { ids[$1]; next } $1 in ids' %s - > %s",
                 chunks, shQuote(keys), shQuote(raw))
  status <- system(cmd)
  if (status != 0) stop("Filtering nominal chunks failed")
} else {
  file.create(raw)
}
nom <- read_qtltools(raw, "nominal")
message(sprintf("%d nominal rows for %d FDR phenotypes", nrow(nom), uniqueN(nom$phe_id)))
write_tsv(nom, snakemake@output[["table"]])
