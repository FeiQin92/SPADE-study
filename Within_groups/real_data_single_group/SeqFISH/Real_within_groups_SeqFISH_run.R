#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
print(args)

Perm_i_min=as.numeric(args[1])*10
Perm_i_max=(as.numeric(args[1])+1)*10

    
# the first 50/200 genes are marker genes
readcounts <- as.matrix(read.table(file=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/real_data_single_group/SeqFISH/count_SeqFISH.spa.txt"), header=T))
info <- read.table(file=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/real_data_single_group/SeqFISH/info_SeqFISH.spa.txt"), row.names=1)
rownames(info) <- colnames(readcounts)
rownames(info) <- colnames(readcounts) 
ED <- as.matrix(dist(info))
lrang <- SPARK::ComputeGaussianPL(ED, compute_distance=FALSE)
colnames(info) <- c("x", "y")
info_back <- info

nPerm <- 1000
for (Perm_i in Perm_i_min:Perm_i_max){

    print(paste0("Perm_", Perm_i))
    #if (Perm_i %% 100==0) cat(Perm_i,"..")
    set.seed(Perm_i)
    sam_indx <- sample(1:nrow(info_back))

    if (Perm_i==0) {info <- info_back}
    if (Perm_i>0) {info$x <- info_back$x[sam_indx]; info$y <- info_back$y[sam_indx]}

    # SPARK method
    library(SPARK)
    spark <- CreateSPARKObject(counts=readcounts,
                                location=info,
                                percentage = 0.1,
                                min_total_counts = 10)
     
    spark@lib_size <- apply(readcounts, 2, sum)
     
    spark <- spark.vc(spark,
                      covariates = NULL,
                      lib_size = spark@lib_size,
                      num_core = 1,
                      verbose = F)
    
    spark_res <- spark.test(object=spark,
                        check_positive = T,
                        verbose = F)
     
    pval <- spark_res@res_mtest

      
    SPARK_res <- spark_res@res_mtest





    ### SPADE estimation 
    library(SPADE)
    regdata <- SPADE_norm(readcounts=as.matrix(readcounts), info=info)
    Est <- SPADE_estimate(expr_data=regdata, info=info)
    res_all_SPADE <- SPADE_test(object=regdata, location=info, para=Est)

    SPADE_list <- list(Est, res_all_SPADE)





    ## MERINGUE
    library(MERINGUE)
    w <- getSpatialNeighbors(info, filterDist = lrang[6])
    #apply(w, 2, sum)
    I <- getSpatialPatterns(readcounts, w)
    results.filter <- filterSpatialPatterns(mat = readcounts,
                                            I = I,
                                            w = w,
                                            adjustPv = TRUE,
                                            alpha = 0.05,
                                            minPercentCells = 0.05,
                                            verbose = TRUE)
       
    MERINGUE_res <- I




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
    
    
    set.seed(2024)
    sce_pre <- spatialPreprocess(sce, platform="ST", n.PCs=7, n.HVGs=50,
                                 log.normalize=T)


    #sce_pre<- qTune(sce_pre, qs=seq(2, 15), platform="ST", d=7)
    
    set.seed(2024)
    sce_cluster <- spatialCluster(sce_pre, q=2, platform="ST", d=7,
                               init.method="mclust", model="t", gamma=2,
                               nrep=1000, burn.in=100,
                               save.chain=TRUE)

    
    ##DESpace
    #library(DESpace)
    library(edgeR)
    files <- list.files("/home/qinf2/GuoshuaiCai/SPADE/DESpace")
    for (j in files){
      source(paste0("/home/qinf2/GuoshuaiCai/SPADE/DESpace/", j))
    }
  
    hotspot_DESpace_BayesSpace <- DESpace_test(spe = sce_cluster,
                                               spatial_cluster = "spatial.cluster")
    head(hotspot_DESpace_BayesSpace$gene_results, 3)
    
    DESpace_BayesSpace.res <-  hotspot_DESpace_BayesSpace$gene_results
    idx <- match(rownames(readcounts), rownames(DESpace_BayesSpace.res))
    DESpace_BayesSpace.res <- DESpace_BayesSpace.res[idx,]
    

    save(SPARK_res, file=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/real_data_single_group/SeqFISH/SPARK/Perm_", Perm_i,"_SPARK_SeqFISH.RData"))
    save(SPADE_list, file=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/real_data_single_group/SeqFISH/SPADE/Perm_", Perm_i,"_SPADE_SeqFISH.RData"))
    save(MERINGUE_res, file=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/real_data_single_group/SeqFISH/MERINGUE/Perm_", Perm_i,"_MERINGUE_SeqFISH.RData"))
    save(DESpace_BayesSpace.res, file=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/real_data_single_group/SeqFISH/DESpace_BayesSpace/Perm_", Perm_i,"_DESpace_BayesSpace_SeqFISH.RData"))

}
