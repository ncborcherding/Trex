# test script for quietTCRgenes.R - testcases are NOT comprehensive!

test_that("quietTCRgenes works", {
  
  data("trex_example")
  SeuratObject::DefaultAssay(trex_example) <- "RNA"
  
  features <- rownames(trex_example@assays$RNA$counts)
  
  expect_equal(
    quietTCRgenes(features),
    getdata("quietTCRgenes", "quietTCRgenes_feature.vector")
  )
  
  SeuratObject::VariableFeatures(trex_example) <- features
  
  trex_example <- quietTCRgenes(trex_example)
  
  expect_equal(
    quietTCRgenes(features),
    SeuratObject::VariableFeatures(trex_example)
  )
})