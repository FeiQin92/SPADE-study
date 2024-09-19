cellProp=0.3
FC=1.4
## Gen position file for permutation test
for (cellProp in c(0.1, 0.2, 0.3)){
  for (FC in c(1.2, 1.3, 1.4, 1.5, 1.6)){
    for (shape in c("hotspot", "streak")){
    	print(paste0("cellProp=", cellProp, " FC=", FC))
    
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

    	info_all <- info_back
        colnames(info_all) <- paste0("Perm_", 0, c("_x", "_y")) 
    	for (Perm_i in 1:1000){

    		set.seed(Perm_i)
    		sam_idx <- sample(1:200)

       		info1 <- info_back[sam_idx,]

    		colnames(info1) <- paste0("Perm_", Perm_i, c("_x", "_y")) 
    		info_all <- cbind(info_all, info1)
    	}
        rownames(info_all) <- paste0("cell", 1:nrow(info_all))

    	write.table(info_all, file=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/simulation/", shape, "/Cell_prop_",cellProp,"/info_", shape, "_200_FC", FC, "_Perm1000.txt"), quote=F, sep="\t", col.names=F)
      } 
   }
}








