#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
print(args)

FC=args[2]
Perm_i_min=as.numeric(args[1])*10
Perm_i_max=(as.numeric(args[1])+1)*10

print(paste0("FC=", FC))

# the first 50/200 genes are marker genes
readcounts <- as.matrix(read.table(file=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/Simu_real_data/MOB/count_MOB_FC", FC, ".txt")))
info <- read.table(file=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/Simu_real_data/MOB/info_MOB_FC", FC, ".txt"))
rownames(info) <- colnames(readcounts)
rownames(info) <- colnames(readcounts) <- paste0("C",1:ncol(readcounts))
rownames(readcounts) <-  paste0("G",1:nrow(readcounts))
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

  actual <- factor(c(rep(1, 50), rep(0, 150)))
  
  pred <- factor(as.numeric(pval$adjusted_pvalue < 0.05))
  library(caret)
  hotspot_SPARK <- confusionMatrix(pred, actual, mode = "everything", positive="1")
  hotspot_SPARK1 <- data.frame(FC=FC, method="SPARK",
                               precision=hotspot_SPARK$byClass["Precision"],
                               Recall=hotspot_SPARK$byClass["Recall"],
                               F1=hotspot_SPARK$byClass["F1"])
  hotspot_SPARK1


  
  
  
  ## SPADE estimation
  library(SPADE)
  regdata <- SPADE_norm(readcounts=as.matrix(readcounts), info=info)
  Est <- SPADE_estimate(expr_data=regdata, info=info)
  res_all_SPADE <- SPADE_test(object=regdata, location=info, para=Est)

  SPADE_list <- list(Est, res_all_SPADE)

  pred <- factor(as.numeric(res_all_SPADE$Adjust.Pvalue < 0.05))

  library(caret)
  test_SPADE <- confusionMatrix(pred, actual, mode = "everything", positive="1")
  hotspot_SPADE1 <- data.frame(FC=FC, method="SPADE",
                              precision=test_SPADE$byClass["Precision"],
                              Recall=test_SPADE$byClass["Recall"],
                              F1=test_SPADE$byClass["F1"])
  hotspot_SPADE1
  
  
  

  
  
  ## MERINGUE
  library(MERINGUE)
  w <- getSpatialNeighbors(info, filterDist = lrang[6])
  apply(w, 2, sum)
  I <- getSpatialPatterns(readcounts, w)
  results.filter <- filterSpatialPatterns(mat = readcounts,
                                          I = I,
                                          w = w,
                                          adjustPv = TRUE,
                                          alpha = 0.05,
                                          minPercentCells = 0.05,
                                          verbose = TRUE)
  
  MERINGUE_res <- list(I, results.filter)
  
  pred <- rep(0, 200)
  pred[as.numeric(substr(results.filter,2,nchar(results.filter)))] <- 1
  pred <- factor(pred)
  
  library(caret)
  test_MERINGUE <- confusionMatrix(pred, actual, mode = "everything", positive="1")
  hotspot_MERINGUE1 <- data.frame(FC=FC, method="MERINGUE", 
                                  precision=test_MERINGUE$byClass["Precision"], 
                                  Recall=test_MERINGUE$byClass["Recall"], 
                                  F1=test_MERINGUE$byClass["F1"])
  
  hotspot_MERINGUE1
  



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
  ##DESpace
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

  pred <- factor(as.numeric(DESpace_BayesSpace.res$FDR < 0.05))

  library(caret)
  test_DESpace_BayesSpace <- confusionMatrix(pred, actual, mode = "everything", positive="1")
  hotspot_DESpace_BayesSpace <- data.frame(FC=FC, method="DESpace_BayesSpace",
                                           precision=test_DESpace_BayesSpace$byClass["Precision"],
                                           Recall=test_DESpace_BayesSpace$byClass["Recall"],
                                           F1=test_DESpace_BayesSpace$byClass["F1"])
  hotspot_DESpace_BayesSpace
  hotspot_MERINGUE1
  hotspot_SPARK1
  hotspot_SPADE1

  
  
  DESpace_BayesSpace_res <- DESpace_BayesSpace.res
  
  save(SPARK_res, file=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/Simu_real_data/MOB/FC",FC,"/SPARK/Perm_", Perm_i,"_SPARK_MOB.RData"))
  save(SPADE_list, file=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/Simu_real_data/MOB/FC",FC,"/SPADE/Perm_", Perm_i,"_SPADE_MOB.RData"))
  save(MERINGUE_res, file=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/Simu_real_data/MOB/FC",FC,"/MERINGUE/Perm_", Perm_i,"_MERINGUE_MOB.RData"))
  save(DESpace_BayesSpace_res, file=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/Simu_real_data/MOB/FC",FC,"/DESpace_BayesSpace/Perm_", Perm_i,"_DESpace_BayesSpace_MOB.RData"))
  
}











