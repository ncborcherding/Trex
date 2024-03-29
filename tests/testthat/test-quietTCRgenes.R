# test script for quietTCRgenes.R - testcases are NOT comprehensive!

test_that("quietTCRgenes works", {
  
  data("trex_example")
  
  features <- rownames(trex_example@assays$RNA$counts)
  
  expect_equal(
    quietTCRgenes(features),
    getdata("quietTCRgenes", "quietTCRgenes_feature.vector")
  )
  
})