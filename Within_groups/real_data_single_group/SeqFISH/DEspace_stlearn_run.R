#!/usr/bin/env Rscript

#args = commandArgs(trailingOnly=TRUE)
#print(args)

Perm_i_min=0
Perm_i_max=1000

# the first 50/200 genes are marker genes
readcounts <- as.matrix(read.table(file=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/real_data_single_group/SeqFISH/count_SeqFISH.spa.txt"), header=T))
info <- read.table(file=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/real_data_single_group/SeqFISH/info_SeqFISH.spa.txt"), row.names=1)
rownames(info) <- colnames(readcounts)
rownames(info) <- colnames(readcounts) 
ED <- as.matrix(dist(info))
lrang <- SPARK::ComputeGaussianPL(ED, compute_distance=FALSE)
colnames(info) <- c("x", "y")
info_back <- info

DESpace_stlearn_list <- list()
 
nPerm <- 1000
for (Perm_i in Perm_i_min:Perm_i_max){
     
    print(paste0("Perm_", Perm_i))
    #if (Perm_i %% 100==0) cat(Perm_i,"..")
    set.seed(Perm_i)
    sam_indx <- sample(1:nrow(info))

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
    library(edgeR)
    files <- list.files("/home/qinf2/GuoshuaiCai/SPADE/DESpace")
    for (j in files){
      source(paste0("/home/qinf2/GuoshuaiCai/SPADE/DESpace/", j))
    }
  
    stLearn_results <- read.csv(paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/real_data_single_group/SeqFISH/stlearn/stlearn_clusters_SeqFISH_Perm_", Perm_i, ".txt"), row.names=1)
    #rownames(stLearn_results) <- gsub("cell", "C", rownames(stLearn_results))
    # Match colData(spe) and stLearn results
    stLearn_results1 <- stLearn_results[match(rownames(colData(sce)), 
                                      rownames(stLearn_results)), ]
    colData(sce)$stLearn_clusters <- stLearn_results1$X_pca_kmeans


    hotspot_DESpace_stlearn <- DESpace_test(spe = sce,
                                        spatial_cluster = "stLearn_clusters")
    head(hotspot_DESpace_stlearn$gene_results, 3)
    DESpace_stlearn.res <-  hotspot_DESpace_stlearn$gene_results
    idx <- match(rownames(readcounts), rownames(DESpace_stlearn.res))
    DESpace_stlearn.res <- DESpace_stlearn.res[idx,]

    DESpace_stlearn_list[[Perm_i+1]] <- DESpace_stlearn.res

    save(DESpace_stlearn.res, file=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/real_data_single_group/SeqFISH/DESpace_stlearn/Perm_", Perm_i,"_DESpace_stlearn_SeqFISH.RData"))

}


