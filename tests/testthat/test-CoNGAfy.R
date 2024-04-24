# test script for CoNGAfy.R - testcases are NOT comprehensive!

test_that("CoNGAfy works", {
  
  data("trex_example")
  
  
  conga_reduction <- CoNGAfy(trex_example)
  
  expect_equal(
    conga_reduction@meta.data,
    getdata("CoNGAfy", "CoNGAfy_meta.data")
  )
  
  expect_equal(
    conga_reduction@assays$RNA@layers$counts,
    getdata("CoNGAfy", "CoNGAfy_counts"),
    tolerance=1e-2
  )
  
  conga_mean_reduction <- CoNGAfy(trex_example,
                                  method = "mean")
  
  expect_equal(
    conga_mean_reduction@meta.data,
    getdata("CoNGAfy", "CoNGAfy_mean_meta.data")
  )
  
  expect_equal(
    conga_mean_reduction@assays$RNA@layers$counts,
    getdata("CoNGAfy", "CoNGAfy_mean_counts"),
    tolerance=1e-2
  )
  
  conga_ADT_reduction <- CoNGAfy(trex_example,
                                 assay = "ADT")
  
  expect_equal(
    conga_ADT_reduction@assays$ADT@layers$counts,
    getdata("CoNGAfy", "CoNGAfy_ADT_counts"),
    tolerance=1e-2
  )
  
  sce <- getdata("runTrex", "SCE.object")
  conga_SCE_reduction <- CoNGAfy(sce,
                                 assay = "RNA")
})