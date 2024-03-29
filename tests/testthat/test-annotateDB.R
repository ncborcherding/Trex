# test script for annotateDB.R - testcases are NOT comprehensive!

test_that("annotateDB works", {
  
  data("trex_example")
  
  trex_example <- annotateDB(trex_example)
  
  expect_equal(
    trex_example@meta.data,
    getdata("annotateDB", "annotateDB_TRB_meta.data")
  )
  
})