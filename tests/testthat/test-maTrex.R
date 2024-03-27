# test script for maTrex.R - testcases are NOT comprehensive!

test_that("maTrex works", {
  data("trex_example")
  
  set.seed(42)
  
  maTrex_result1 <- maTrex(trex_example, 
                           chains = "TRB",
                           method = "encoder",
                           encoder.model = "VAE",
                           encoder.input = "AF")
  
  expect_equal(
    maTrex_result1,
    getdata("runTrex", "maTrex_TRB_VAE_AF", tolerance=1e-3)
  )
  
  maTrex_result2 <- maTrex(trex_example, 
                           chains = "TRA",
                           method = "encoder",
                           encoder.model = "AE",
                           encoder.input = "OHE")
  
  expect_equal(
    maTrex_result2,
    getdata("runTrex", "maTrex_TRA_AE_OHE", tolerance=1e-3)
  )
})