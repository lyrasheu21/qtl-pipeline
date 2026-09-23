# pi1 replication: look up each FDR-significant cciQTL (lead SNP, anchor gene,
# anchor subclass) in the eQTL resource for that subclass, then estimate
# pi1 = 1 - pi0 from the matched eQTL p-values.

for (f in list.files(snakemake@params[["rdir"]], pattern = "\\.R$", full.names = TRUE)) source(f)
start_logging(snakemake@log[[1]])

s <- snakemake@params[["settings"]]
qt <- rbind(fread(snakemake@input[["sending"]]), fread(snakemake@input[["receiving"]]))
hits <- qt[fdr_pass == TRUE, .(phe_id, direction, subclass, gene, snp = var_id, chr = var_chr,
                               cciqtl_beta = slope, cciqtl_q = qval)]
message(sprintf("FDR cciQTLs to look up: %d", nrow(hits)))

has_header <- isTRUE(s$has_header)
col_names <- unlist(s$columns)
read_eqtl <- function(path, pairs) {
  # Filter the (possibly large) eQTL file to the SNP-gene pairs of interest with awk.
  header <- if (has_header) strsplit(readLines(gzcon(file(path, "rb")), n = 1), "\t")[[1]] else col_names
  snp_i <- match(s$snp_column, header)
  gene_i <- match(s$gene_column, header)
  if (is.na(snp_i) || is.na(gene_i)) stop("SNP/gene column not found in ", path)
  pair_file <- tempfile()
  writeLines(pairs, pair_file)
  cmd <- sprintf("gzip -dcf %s | awk -F'\\t' 'NR == FNR { t[$1]; next } ($%d \"_\" $%d) in t' %s -",
                 shQuote(path), snp_i, gene_i, pair_file)
  lines <- system(cmd, intern = TRUE)
  unlink(pair_file)
  if (!length(lines)) return(NULL)
  dt <- fread(text = paste(lines, collapse = "\n"), header = FALSE, sep = "\t")
  setnames(dt, header[seq_len(ncol(dt))])
  dt
}

groups <- unique(hits[, .(subclass, chr)])
eqtl <- rbindlist(lapply(seq_len(nrow(groups)), function(i) {
  sc <- groups$subclass[i]
  chr <- groups$chr[i]
  path <- fill_pattern(s$eqtl_pattern, subclass = sc, chr = chr)
  if (!file.exists(path)) {
    message("eQTL file missing, skipped: ", path)
    return(NULL)
  }
  pairs <- unique(hits[subclass == sc & chr == groups$chr[i], paste0(snp, "_", gene)])
  dt <- read_eqtl(path, pairs)
  if (is.null(dt)) return(NULL)
  p <- if (nzchar(s$p_column %||% "")) as.numeric(dt[[s$p_column]]) else 2 * pnorm(-abs(as.numeric(dt[[s$z_column]])))
  data.table(subclass = sc, snp = as.character(dt[[s$snp_column]]), gene = as.character(dt[[s$gene_column]]),
             eqtl_beta = as.numeric(dt[[s$beta_column]]), eqtl_p = p)
}))
if (nrow(eqtl) == 0) {
  eqtl <- data.table(subclass = character(), snp = character(), gene = character(),
                     eqtl_beta = numeric(), eqtl_p = numeric())
}
eqtl <- unique(eqtl, by = c("subclass", "snp", "gene"))

matched <- merge(hits, eqtl, by = c("subclass", "snp", "gene"), all.x = TRUE)
write_tsv(matched, snakemake@output[["matched"]])

pi1 <- rbind(
  matched[, .(n_cciqtl = .N, n_matched = sum(!is.na(eqtl_p)), pi1 = pi1_estimate(eqtl_p)), by = direction],
  matched[, .(direction = "both", n_cciqtl = .N, n_matched = sum(!is.na(eqtl_p)), pi1 = pi1_estimate(eqtl_p))]
)
print(pi1)
write_tsv(pi1, snakemake@output[["pi1"]])

plot_data <- matched[!is.na(eqtl_beta)]
if (nrow(plot_data) < 3) {
  empty_plot(snakemake@output[["scatter"]], "Fewer than 3 cciQTLs matched in the eQTL resource")
} else {
  rho <- plot_data[, .(rho = suppressWarnings(cor(eqtl_beta, cciqtl_beta, method = "spearman"))), by = direction]
  labels <- setNames(sprintf("%s (Spearman rho = %.2f)", rho$direction, rho$rho), rho$direction)
  p <- ggplot(plot_data, aes(eqtl_beta, cciqtl_beta, colour = direction)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey40") +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, formula = y ~ x) +
    scale_colour_manual(values = DIRECTION_COLORS, labels = labels) +
    labs(x = "eQTL beta", y = "cciQTL beta") +
    theme_pipeline()
  save_plot(p, snakemake@output[["scatter"]], width = 7, height = 5)
}
