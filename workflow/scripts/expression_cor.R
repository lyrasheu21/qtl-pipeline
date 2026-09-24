# Correlation across donors between each phenotype of one set and the
# expression of its gene in its group (qtlpost::phenotype_expression_cor).
source(file.path(snakemake@scriptdir, "load_qtlpost.R"))

pattern <- snakemake@params[["expression_pattern"]]
bed <- read_bed_matrix(snakemake@input[["bed"]])
annotation <- fread(snakemake@input[["annotation"]])
get_expression <- function(g) {
  path <- fill_pattern(pattern, group = g, subclass = g)
  if (!file.exists(path)) { message("No expression file for ", g, ": ", path); return(NULL) }
  read_expression(path)
}
res <- phenotype_expression_cor(bed$mat, annotation, get_expression, snakemake@params[["min_donors"]])
res[, set := snakemake@wildcards[["set"]]]
message(sprintf("%d phenotype-gene correlations", nrow(res)))
write_tsv(res, snakemake@output[[1]])
