#!/usr/bin/env Rscript


#cellProp=0.3
#FC=1.4
#shape="hotspot"
for (cellProp in c(0.1, 0.2, 0.3)){
  for (FC in c(1.2, 1.3, 1.4, 1.5, 1.6)){
    for (shape in c("hotspot", "streak")){
    print(paste0(shape, ": cellProp=", cellProp, " FC=", FC))
    
    # the first 50/200 genes are marker genes
    readcounts <- as.matrix(read.table(file=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/simulation/", shape, "/Cell_prop_",cellProp,"/count_", shape, "_200_FC", FC, ".txt")))
    info <- read.table(file=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/simulation/", shape, "/Cell_prop_",cellProp,"/info_", shape, "_200_FC", FC, ".txt"))
    rownames(info) <- colnames(readcounts)
    rownames(info) <- colnames(readcounts) <- paste0("C",1:ncol(readcounts))
    rownames(readcounts) <-  paste0("G",1:nrow(readcounts))
    ED <- as.matrix(dist(info))
    lrang <- SPARK::ComputeGaussianPL(ED, compute_distance=FALSE)
    colnames(info) <- c("x", "y")
    info_back <- info

    DESpace_stlearn_list <- list()
 

    nPerm <- 1000
    for (Perm_i in c(0:nPerm)){
     
    print(paste0("Perm_", Perm_i))
    #if (Perm_i %% 100==0) cat(Perm_i,"..")
    set.seed(Perm_i)
    sam_indx <- sample(1:200)

    if (Perm_i==0) {info <- info_back}
    if (Perm_i>0) {info$x <- info_back$x[sam_indx]; info$y <- info_back$y[sam_indx]}

    ## DESpace_BayesSpace
    # BiocManager::install("DESpace")
    library(SingleCellExperiment)
    library(BayesSpace)
    colData <- data.frame(row=info$x, col=info$y)
    rownames(colData) <- colnames(readcounts)
    rowData <- data.frame(GeneID=rownames(readcounts))
    # rowData <- 
    sce <- SingleCellExperiment(assays=list(counts=as.matrix(readcounts)),
                                rowData=rowData,
                                colData=colData)

    ##DESpace
    library(DESpace)
  
    stLearn_results <- read.csv(paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/simulation/", shape, "/Cell_prop_",cellProp,"/FC", FC, "/stlearn/stlearn_clusters_", shape, "_FC", FC, "_Perm_", Perm_i, ".txt"), row.names=1)
    rownames(stLearn_results) <- gsub("cell", "C", rownames(stLearn_results))
    # Match colData(spe) and stLearn results
    stLearn_results1 <- stLearn_results[match(rownames(colData(sce)), 
                                      rownames(stLearn_results)), ]
    colData(sce)$stLearn_clusters <- stLearn_results1$X_pca_kmeans



	library(DESpace)
	hotspot_DESpace_stlearn <- DESpace_test(spe = sce,
                                        spatial_cluster = "stLearn_clusters")
	head(hotspot_DESpace_stlearn$gene_results, 3)
	DESpace_stlearn.res <-  hotspot_DESpace_stlearn$gene_results
	idx <- match(rownames(readcounts), rownames(DESpace_stlearn.res))
	DESpace_stlearn.res <- DESpace_stlearn.res[idx,]

      DESpace_stlearn_list[[Perm_i+1]] <- DESpace_stlearn.res

      names(DESpace_stlearn_list)[Perm_i+1] <- paste0("Perm_i_", Perm_i)
  }
 save(DESpace_stlearn_list, file=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/simulation/", shape, "/Cell_prop_",cellProp,"/FC",FC,"/Perm_DESpace_stlearn_", shape, ".RData"))

}
}
}

