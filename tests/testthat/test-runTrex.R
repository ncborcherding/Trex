# test script for runTrex.R - testcases are NOT comprehensive!

test_that("runTrex works with seurat objects", {
	data("trex_example")
  set.seed(42)
  
  trex_example <- runTrex(trex_example, 
                          chains = "TRA",
                          method = "encoder",
                          encoder.model = "VAE",
                          encoder.input = "AF", 
                          reduction.name = "TRA_VAE_AF")
                       
	expect_equal(
		trex_example@reductions$TRA_VAE_AF@cell.embeddings,
		getdata("runTrex", "runTrex_TRA_VAE_AF_reduction"),
		tolerance=1e-2
	)
	
	trex_example <- runTrex(trex_example, 
	                        chains = "TRB",
	                        method = "encoder",
	                        encoder.model = "AE",
	                        encoder.input = "KF", 
	                        reduction.name = "TRB_AE_KF")
	
	expect_equal(
	  trex_example@reductions$TRB_AE_KF@cell.embeddings,
	  getdata("runTrex", "runTrex_TRB_AE_kF_reduction"),
	  tolerance=1e-2
	)
	
	trex_example <- runTrex(trex_example, 
	                        chains = "TRB",
	                        method = "encoder",
	                        encoder.model = "VAE",
	                        encoder.input = "OHE", 
	                        reduction.name = "TRB_VAE_OHE")
	
	expect_equal(
	  trex_example@reductions$TRB_VAE_OHE@cell.embeddings,
	  getdata("runTrex", "runTrex_TRB_VAE_OHE_reduction"),
	  tolerance=1e-2
	)
	
	trex_example <- runTrex(trex_example, 
	                        chains = "TRB",
	                        method = "geometric",
	                        reduction.name = "TRB_Geometric")
	
	expect_equal(
	  trex_example@reductions$TRB_Geometric@cell.embeddings,
	  getdata("runTrex", "runTrex_TRB_geometric_reduction"),
	  tolerance=1e-2
	)
	

	sce <- getdata("runTrex", "SCE.object")
	sce <- runTrex(sce, 
	               chains = "TRA",
	               method = "geometric",
	               reduction.name = "TRA_Geometric")
	
	expect_equal(
	  sce@int_colData$reducedDims$TRA_Geometric,
	  getdata("runTrex", "runTrex_SCE_TRA_geometric_reduction"),
	  tolerance=1e-2
	)
	
})
