# test script for quietTCRgenes.R - testcases are NOT comprehensive!

test_that("quietTCRgenes works", {
  
  data("trex_example")
  
  features <- rownames(trex_example@assays$RNA$counts)
  
  expect_equal(
    quietTCRgenes(features),
    getdata("quietTCRgenes", "quietTCRgenes_feature.vector")
  )
  
  trex_example@assays$RNA@var.features <- features
  Seurat::DefaultAssay(trex_example) <- "RNA"
  
  trex_example <- quietTCRgenes(trex_example)
  
  
  expect_equal(
    quietTCRgenes(features),
    trex_example@assays$RNA@var.features
  )
})