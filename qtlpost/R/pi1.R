# pi1 replication of QTL hits in an external eQTL resource.

#' Look up QTL hits in an eQTL resource and estimate pi1 per set.
#'
#' @param hits data.table with columns set, group, gene, snp, chr, qtl_beta.
#' @param pattern eQTL file pattern with `{group}` (or `{subclass}`) and `{chr}` placeholders.
#' @param columns column names of the eQTL files when they have no header.
#' @param snp_column,gene_column,beta_column,z_column,p_column eQTL column names;
#'   give `p_column` or `z_column` (p is then computed from z).
#' @return list(matched = hits with eqtl_beta/eqtl_p, pi1 = per-set and overall pi1)
pi1_replication <- function(hits, pattern, has_header = FALSE, columns = NULL,
                            snp_column = "snp", gene_column = "gene", beta_column = "beta",
                            z_column = "z", p_column = "") {
  read_eqtl <- function(path, pairs) {
    header <- if (has_header) strsplit(readLines(gzcon(file(path, "rb")), n = 1), "\t")[[1]] else columns
    snp_i <- match(snp_column, header)
    gene_i <- match(gene_column, header)
    if (is.na(snp_i) || is.na(gene_i)) stop("SNP/gene column not found in ", path)
    pair_file <- tempfile()
    writeLines(pairs, pair_file)
    on.exit(unlink(pair_file))
    # filter the (possibly large) file to the SNP-gene pairs of interest with awk
    cmd <- sprintf("gzip -dcf %s | awk -F'\\t' 'NR == FNR { t[$1]; next } ($%d \"_\" $%d) in t' %s -",
                   shQuote(path), snp_i, gene_i, pair_file)
    lines <- system(cmd, intern = TRUE)
    if (!length(lines)) return(NULL)
    dt <- fread(text = paste(lines, collapse = "\n"), header = FALSE, sep = "\t")
    setnames(dt, header[seq_len(ncol(dt))])
    dt
  }

  lookups <- unique(hits[, .(group, chr)])
  eqtl <- rbindlist(lapply(seq_len(nrow(lookups)), function(i) {
    g <- lookups$group[i]
    chr <- lookups$chr[i]
    path <- fill_pattern(pattern, group = g, subclass = g, chr = chr)
    if (!file.exists(path)) {
      message("eQTL file missing, skipped: ", path)
      return(NULL)
    }
    pairs <- unique(hits[group == g & chr == lookups$chr[i], paste0(snp, "_", gene)])
    dt <- read_eqtl(path, pairs)
    if (is.null(dt)) return(NULL)
    p <- if (nzchar(p_column %||% "")) as.numeric(dt[[p_column]]) else 2 * pnorm(-abs(as.numeric(dt[[z_column]])))
    data.table(group = g, snp = as.character(dt[[snp_column]]), gene = as.character(dt[[gene_column]]),
               eqtl_beta = as.numeric(dt[[beta_column]]), eqtl_p = p)
  }))
  if (nrow(eqtl) == 0) {
    eqtl <- data.table(group = character(), snp = character(), gene = character(),
                       eqtl_beta = numeric(), eqtl_p = numeric())
  }
  eqtl <- unique(eqtl, by = c("group", "snp", "gene"))
  matched <- merge(hits, eqtl, by = c("group", "snp", "gene"), all.x = TRUE)
  pi1 <- rbind(
    matched[, .(n_hits = .N, n_matched = sum(!is.na(eqtl_p)), pi1 = pi1_estimate(eqtl_p)), by = set],
    matched[, .(set = "all", n_hits = .N, n_matched = sum(!is.na(eqtl_p)), pi1 = pi1_estimate(eqtl_p))]
  )
  list(matched = matched, pi1 = pi1)
}
