#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
print(args)

cellProp=args[1]
FC=args[2]

#cellProp=0.3
#FC=1.4
#Eval.final <- NULL
#for (cellProp in c(0.1, 0.2, 0.3)){
#  for (FC in c(1.2, 1.3, 1.4, 1.5, 1.6)){
    print(paste0("cellProp=", cellProp, " FC=", FC))
    
    # the first 50/200 genes are marker genes
    # Each row is a gene, and each column is a cell
    readcounts <- as.matrix(read.table(file=paste0("../simulation/streak/Cell_prop_",cellProp,"/count_streak_200_FC", FC, ".txt")))
    info <- read.table(file=paste0("../simulation/streak/Cell_prop_",cellProp,"/info_streak_200_FC", FC, ".txt"))
    rownames(info) <- colnames(readcounts)
    rownames(info) <- colnames(readcounts) <- paste0("C",1:ncol(readcounts))
    rownames(readcounts) <-  paste0("G",1:nrow(readcounts))
    ED <- as.matrix(dist(info))
    lrang <- SPARK::ComputeGaussianPL(ED, compute_distance=FALSE)
    colnames(info) <- c("x", "y")
    info_back <- info

    SPARK_list <- list()
    MERINGUE_list <- list()
    SPADE_list <- list()
    DESpace_BayesSpace_list <- list()

    nPerm <- 1000
    for (Perm_i in c(0:1000)){
     
    	print(paste0("Perm_", Perm_i))
    	#if (Perm_i %% 100==0) cat(Perm_i,"..")
    	set.seed(Perm_i)
    	sam_indx <- sample(1:200)

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

    	actual <- factor(c(rep(1, 50), rep(0, 150)))
    	## F1 score to evaluate results
    	pred <- factor(as.numeric(pval$adjusted_pvalue < 0.05))
    	library(caret)
    	streak_SPARK <- confusionMatrix(pred, actual, mode = "everything", positive="1")
    	streak_SPARK1 <- data.frame(cellProp=cellProp, FC=FC, method="SPARK",
                             precision=streak_SPARK$byClass["Precision"],
                             Recall=streak_SPARK$byClass["Recall"],
                             F1=streak_SPARK$byClass["F1"])
    	streak_SPARK1

    	SPARK_list[[Perm_i+1]] <- spark_res@res_mtest


    	## SPADE method
    	library(SPADE)
    	regdata <- SPADE_norm(readcounts=as.matrix(readcounts), info=info)

    	Est <- SPADE_estimate(expr_data=regdata, info=info)
    	res_all_SPADE <- SPADE_test(object=regdata, location=info, para=Est)

    	actual <- factor(c(rep(1, 50), rep(0, 150)))
   	pred <- factor(as.numeric(res_all_SPADE$Adjust.Pvalue < 0.05))

    	library(caret)
    	test_SPADE <- confusionMatrix(pred, actual, mode = "everything", positive="1")
    	streak_SPADE1 <- data.frame(cellProp=cellProp, FC=FC, method="SPADE",
                             precision=test_SPADE$byClass["Precision"],
                             Recall=test_SPADE$byClass["Recall"],
                             F1=test_SPADE$byClass["F1"])
    	streak_SPADE1
    	streak_SPARK1
    	SPADE_list[[Perm_i+1]] <- list(Est, res_all_SPADE)



    	## MERINGUE method
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
    	actual <- factor(c(rep(1, 50), rep(0, 150)))
    	pred <- rep(0, 200)
    	pred[as.numeric(substr(results.filter,2,nchar(results.filter)))] <- 1
    	pred <- factor(pred)

    	library(caret)
    	test_MERINGUE <- confusionMatrix(pred, actual, mode = "everything", positive="1")
   	streak_MERINGUE1 <- data.frame(cellProp=cellProp, FC=FC, method="MERINGUE", 
                                precision=test_MERINGUE$byClass["Precision"], 
                                Recall=test_MERINGUE$byClass["Recall"], 
                                F1=test_MERINGUE$byClass["F1"])

    	streak_MERINGUE1

    	MERINGUE_list[[Perm_i+1]] <- list(I, results.filter)



    	# DESpace_BayesSpace
    	#BiocManager::install("DESpace")
    	library(SingleCellExperiment)
    	library(BayesSpace)
    	colData <- data.frame(row=info$x, col=info$y)
    	rownames(colData) <- colnames(readcounts)
    	rowData <- data.frame(GeneID=rownames(readcounts))

    	sce <- SingleCellExperiment(assays=list(counts=as.matrix(readcounts)),
                            rowData=rowData,
                            colData=colData)


    	set.seed(2024)
    	sce_pre <- spatialPreprocess(sce, platform="ST", n.PCs=7, n.HVGs=50,
                             log.normalize=T)

    	set.seed(2024)
    	sce_cluster <- spatialCluster(sce_pre, q=2, platform="ST", d=7,
                              init.method="mclust", model="t", gamma=2,
                              nrep=1000, burn.in=100,
                              save.chain=TRUE)


    	##DESpace
    	library(DESpace)

    	streak_DESpace_BayesSpace <- DESpace_test(spe = sce_cluster,
                                           spatial_cluster = "spatial.cluster")
    	head(streak_DESpace_BayesSpace$gene_results, 3)

    	DESpace_BayesSpace.res <-  streak_DESpace_BayesSpace$gene_results
    	idx <- match(rownames(readcounts), rownames(DESpace_BayesSpace.res))
    	DESpace_BayesSpace.res <- DESpace_BayesSpace.res[idx,]

    	actual <- factor(c(rep(1, 50), rep(0, 150)))
    	pred <- factor(as.numeric(DESpace_BayesSpace.res$FDR < 0.05))

    	library(caret)
    	test_DESpace_BayesSpace <- confusionMatrix(pred, actual, mode = "everything", positive="1")
    	streak_DESpace_BayesSpace <- data.frame(cellProp=cellProp, FC=FC, method="DESpace_BayesSpace",
                                         precision=test_DESpace_BayesSpace$byClass["Precision"],
                                         Recall=test_DESpace_BayesSpace$byClass["Recall"],
                                         F1=test_DESpace_BayesSpace$byClass["F1"])
    	streak_DESpace_BayesSpace

    	DESpace_BayesSpace_list[[Perm_i+1]] <- DESpace_BayesSpace.res


    	names(SPARK_list)[Perm_i+1] <-  paste0("Perm_i_", Perm_i)
    	names(MERINGUE_list)[Perm_i+1] <-  paste0("Perm_i_", Perm_i)
    	names(SPADE_list)[Perm_i+1] <-  paste0("Perm_i_", Perm_i)
    	names(DESpace_BayesSpace_list)[Perm_i+1] <-  paste0("Perm_i_", Perm_i)
    }


    save(SPARK_list, file=paste0("../simulation/streak/Cell_prop_",cellProp,"/FC",FC,"/Perm_SPARK_streak.RData"))
    save(SPADE_list, file=paste0("../simulation/streak/Cell_prop_",cellProp,"/FC",FC,"/Perm_SPADE_streak.RData"))
    save(MERINGUE_list, file=paste0("../simulation/streak/Cell_prop_",cellProp,"/FC",FC,"/Perm_MERINGUE_streak.RData"))
    save(DESpace_BayesSpace_list, file=paste0("../simulation/streak/Cell_prop_",cellProp,"/FC",FC,"/Perm_DESpace_BayesSpace_streak.RData"))

    #Eval.final <- rbind(Eval.final, streak_SPADE1, streak_SPARK1, streak_MERINGUE1, streak_DESpace_BayesSpace)

#  }
#}

