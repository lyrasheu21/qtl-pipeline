test_that("gene track finds genes in the window (argument names shadow annotation columns)", {
  annot <- data.table(gene_name = c("A", "B", "C"), gene_id = c("g1", "g2", "g3"),
                      chrom = c("chr1", "chr1", "chr2"), start = c(100, 5000, 100),
                      end = c(200, 6000, 200), strand = "+")
  annot[, chrom_norm := sub("^chr", "", chrom)]
  p <- gene_track(annot, 1L, 0, 1000, highlight = "A")
  expect_equal(p$data$gene_name, "A")
  expect_equal(gene_track(annot, "chr1", 0, 10000)$data$gene_name, c("A", "B"))
  expect_length(gene_track(annot, 3, 0, 1000)$layers, 0)  # no genes: empty track
})
