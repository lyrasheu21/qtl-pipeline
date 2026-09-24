test_that("sample maps rename columns and average replicates", {
  mat <- matrix(c(1, 10, 2, 20, 3, 30), nrow = 2, dimnames = list(c("p1", "p2"), c("s1", "s2", "s3")))
  map <- data.frame(sample = c("s1", "s2", "s3"), donor = c("D1", "D1", "D2"))
  out <- apply_sample_map(mat, map, "sample", "donor")
  expect_equal(colnames(out), c("D1", "D2"))
  expect_equal(out["p1", "D1"], mean(c(1, 2)))
  expect_equal(out["p2", "D2"], 30)
})

test_that("unmapped samples are dropped", {
  mat <- matrix(1:4, nrow = 1, dimnames = list("p", c("s1", "s2", "s3", "s4")))
  out <- apply_sample_map(mat, data.frame(a = c("s1", "s3"), b = c("D1", "D3")), "a", "b")
  expect_equal(colnames(out), c("D1", "D3"))
})

test_that("factor files are read in either orientation", {
  samples <- c("D1", "D2", "D3")
  rows_are_samples <- tempfile(fileext = ".tsv")
  fwrite(data.table(SubID = samples, PC1 = 1:3, PC2 = 4:6), rows_are_samples, sep = "\t")
  a <- read_sample_factors(rows_are_samples, samples)
  rows_are_pcs <- tempfile(fileext = ".pca")
  fwrite(data.table(SampleID = c("x_PC1", "x_PC2"), D1 = c(1, 4), D2 = c(2, 5), D3 = c(3, 6)), rows_are_pcs, sep = "\t")
  b <- read_sample_factors(rows_are_pcs, samples)
  expect_equal(rownames(b), samples)
  expect_equal(unname(a[, 1]), unname(b[, 1]))
  expect_equal(unname(a[, 2]), unname(b[, 2]))
})

test_that("factors can be read from an RData object", {
  hidden_factors <- matrix(1:6, nrow = 3, dimnames = list(c("S1", "S2", "S3"), NULL))
  path <- tempfile(fileext = ".RData")
  save(hidden_factors, file = path)
  f <- read_sample_factors(path, c("S1", "S2", "S3"), object = "hidden_factors")
  expect_equal(dim(f), c(3L, 2L))
  expect_equal(colnames(f), c("F1", "F2"))
})

test_that("fill_pattern fills named placeholders", {
  expect_equal(fill_pattern("e_{group}_chr{chr}.gz", group = "Astro", chr = 7), "e_Astro_chr7.gz")
})

test_that("BED files with ';' in phenotype IDs are read as tab-separated", {
  path <- tempfile(fileext = ".bed")
  writeLines(c("#chr\tstart\tend\tpid\tgid\tstrand\tD1\tD2",
               "1\t100\t101\tchr1:100-101;+\tGENE\t+\t1.5\t2.5"), path)
  b <- read_bed_matrix(path)
  expect_equal(rownames(b$mat), "chr1:100-101;+")
  expect_equal(unname(b$mat[1, ]), c(1.5, 2.5))
})
