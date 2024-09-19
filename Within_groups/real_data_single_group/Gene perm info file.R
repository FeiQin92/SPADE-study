## Gen position file for permutation test
for (dataset in c("MOB", "SeqFISH")){
    
        readcounts <- as.matrix(read.table(file=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/real_data_single_group/", dataset, "/count_", dataset, ".spa.txt"), header=T))
        info <- read.table(file=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/real_data_single_group/", dataset, "/info_", dataset, ".spa.txt"), row.names=1)
        all(colnames(readcounts)==rownames(info))
   		colnames(info) <- c("x", "y")
   		info_back <- info

    	info_all <- info_back
        colnames(info_all) <- paste0("Perm_", 0, c("_x", "_y")) 
    	for (Perm_i in 1:1000){

    		set.seed(Perm_i)
    		sam_idx <- sample(1:nrow(info))

       		info1 <- info_back[sam_idx,]

    		colnames(info1) <- paste0("Perm_", Perm_i, c("_x", "_y")) 
    		info_all <- cbind(info_all, info1)
    	}
        rownames(info_all) <- rownames(info)

    	write.table(info_all, file=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/real_data_single_group/", dataset, "/info_", dataset, "_Perm1000.txt"), quote=F, sep="\t", col.names=F)
} 










