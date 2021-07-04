#AF = Atchley factors, KF = Kidera factors


runTrex <- function(sc, 
                   chains = "both", 
                   edit.method = "lv",
                   AA.propeties = c("AF", "KF", "other"),
                   c.trim = 0,
                   n.trim = 0,
                   threshold = 0.85,
                   add.INKT = TRUE,
                   add.MAIT = TRUE, 
                   species = "human") {
    TCR <- getTCR(sc, chains)
    print("Calculating the Edit Distance for CDR3 AA sequence...")
    multi.network <- distanceMatrix(TCR, c.trim, n.trim, threshold)
    
    
    if (add.INKT) {
        print("Calculating the INKT gene usage...")
        tmpscore <- scoreINKT(TCR, species = species)
        if (length(which(tmpscore$score > 0)) != 0) {
            tmp.knn <- gene.to.knn(tmpscore)
            multi.network <- add.to.network(multi.network, tmp.knn, "INKT") 
        }
        barcodes <- tmpscore$barcode
        tmpscore <- as.data.frame(tmpscore[,2])
        colnames(tmpscore) <- "IKNT.score"
        rownames(tmpscore) <- barcodes
        sc <- add.meta.data(sc, tmpscore)
    }
    if (add.MAIT) {
        print("Calculating the MAIT gene usage...")
        tmpscore <- scoreMAIT(TCR, species = species)
        if (length(which(tmpscore$score > 0)) != 0) {
            tmp.knn <- gene.to.knn(tmpscore)
            multi.network <- add.to.network(multi.network, tmp.knn, "INKT") 
        }
        barcodes <- tmpscore$barcode
        tmpscore <- as.data.frame(tmpscore[,2])
        colnames(tmpscore) <- "MAIT.score"
        rownames(tmpscore) <- barcodes
        sc <- add.meta.data(sc, tmpscore)
    }

}