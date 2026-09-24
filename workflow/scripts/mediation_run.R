# Mediation (genotype -> anchor expression -> partner expression) and
# interaction models for the selected hits (qtlpost::run_mediation).
source(file.path(snakemake@scriptdir, "load_qtlpost.R"))

set.seed(snakemake@params[["seed"]])
pattern <- snakemake@params[["expression_pattern"]]
hits <- read_tsv(snakemake@input[["hits"]])
genotypes <- read_tsv(snakemake@input[["genotypes"]], colClasses = "character")

cache <- new.env()
get_expression <- function(g) {
  if (!exists(g, envir = cache, inherits = FALSE)) {
    path <- fill_pattern(pattern, group = g, subclass = g)
    if (!file.exists(path)) message("Expression file missing: ", path)
    assign(g, if (file.exists(path)) read_expression(path) else NULL, envir = cache)
  }
  get(g, envir = cache)
}

res <- run_mediation(hits, genotypes, get_expression, sims = snakemake@params[["sims"]])
message(sprintf("Mediation results: %d (ACME p < 0.05: %d)", nrow(res$mediation), sum(res$mediation$acme_p < 0.05)))
write_tsv(res$mediation, snakemake@output[["mediation"]])
write_tsv(res$interaction, snakemake@output[["interaction"]])
