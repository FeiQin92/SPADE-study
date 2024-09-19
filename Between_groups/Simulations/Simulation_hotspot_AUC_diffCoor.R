library("dplyr")
library("gstat")
library("sp")
library("spacetime")
library("STRbook")
library("tidyr")

library(spatstat)
library(Matrix)
library(MASS)
library(FRK)
source("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\Code\\utilities.R")


setwd("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\data\\SPARK\\raw_data")

MOB <- t(read.table("Rep11_MOB_count_matrix-1.tsv"))
MOB1 <- MOB[,apply(MOB, 2, sum) >= 10]
MOB2 <- MOB1[(apply(MOB1>0, 1, mean) >= 0.3) ,]
set.seed(2022)
MOB3 <- MOB2[, sample(1:ncol(MOB2), 200)]

library(splatter)
# library(SCRIP)
params <- splatEstimate(MOB3)

sim_trend <-  splatSimulate(params=params)
sim_trend_RC <- counts(sim_trend)
save(sim_trend_RC, file="E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\hotspot\\simu.data.400Cells.RData")
sim_trend_RC <- get(load("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\hotspot\\simu.data.400Cells.RData"))
sim_trend_RC1 <- sim_trend_RC[apply(sim_trend_RC>1, 1, mean) >= 0.3,]
dim(sim_trend_RC1)

numGenes = 200
lowSignal = 25
highSignal = 25

numCell = 200

for (cellProp in c(0.1, 0.2, 0.3)){
  for (ieffect in c(1.5, 2.0, 2.5, 3.0, 3.5)){
    # itheta = 10
    set.seed(31)
    sim_trend_RC400 <- sim_trend_RC1[sample(1:nrow(sim_trend_RC1), 400),]
    
    low_expr    <- c(1)
    high_expr   <- c(2)

    set.seed(31)
    pp_T <- rpoispp(0.00020,  win=owin(c(1,1000),c(1,1000)))
    pp_T = add_markdist_hotspot(pp_T,low_marks=low_expr, high_marks=high_expr, nMarkers=(lowSignal+highSignal), 
                                ncell=numGenes, cell_proportion=cellProp)
    
    set.seed(69)
    pp_C <- rpoispp(0.00020,  win=owin(c(1,1000),c(1,1000)))
    pp_C = add_markdist_hotspot(pp_C,low_marks=low_expr, high_marks=high_expr, nMarkers=(lowSignal+highSignal), 
                                ncell=numGenes, cell_proportion=cellProp)
    
    expr_mat_C <- matrix(NA,nrow=pp_C$n,ncol=numGenes)
    expr_mat_T <- matrix(NA,nrow=pp_T$n,ncol=numGenes)
    for(igene in 1:numGenes){
      
      spiked_in_cells_T   <- which(pp_T$markformat[,igene]==2)
      spiked_in_cells_C   <- which(pp_C$markformat[,igene]==2)
      
      a <- sim_trend_RC400[igene,]
      if(igene <= highSignal){
        a[spiked_in_cells_T]  <- sim_trend_RC400[igene,spiked_in_cells_T]*ieffect
      }else if(igene <=lowSignal+highSignal) {
        a[spiked_in_cells_T]  <- sim_trend_RC400[igene,spiked_in_cells_T]/ieffect
      }
      
      b <- sim_trend_RC400[igene+200, 1:pp_C$n]
      
      expr_mat_T[,igene] <- a
      expr_mat_C[,igene] <- b
      
   
    }
    
    colnames(expr_mat_T) <- colnames(expr_mat_C)  <- paste0("gene",1:numGenes)
    info_T                <- cbind.data.frame(x=pp_T$x, y=pp_T$y)
    info_T <- round(info_T)

    info_C                <- cbind.data.frame(x=pp_C$x, y=pp_C$y)
    info_C <- round(info_C)

    
    rownames(info_T) <-  rownames(expr_mat_T) <- paste0("C",1:nrow(expr_mat_T))
    rownames(info_C) <- rownames(expr_mat_C) <- paste0("C",1:nrow(expr_mat_C))
    count_T=t(expr_mat_T)
    count_C=t(expr_mat_C)
    
    write.table(pp_T$markformat,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\count_hotspotEX_200_FC_hr_indexT", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")
    write.table(pp_C$markformat,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\count_hotspotEX_200_FC_hr_indexC", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")
    
    write.table(count_T,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\count_hotspotEX_200_FC_hr_T", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")
    write.table(count_C,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\count_hotspotEX_200_FC_hr_C", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")
    
    write.table(info_T,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\info_hotspotEX_200_FC_hr_T", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")
    write.table(info_C,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\info_hotspotEX_200_FC_hr_C", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")
      }
}







## SPADE method with LLR test based on exchange both length scale theta and Tao 1 (both shape and strength)
############################################################################
###########################################################################

cellProp=0.2
FC=2.0
for (cellProp in c(0.1, 0.2, 0.3)){
  for (FC in c(1.5, 2.0, 2.5, 3.0, 3.5)){
    print(paste0("cellProp=", cellProp, "; FC=", FC))
    # the first 100/300 genes are marker genes
    readcounts_T <- as.matrix(read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\count_hotspotEX_200_FC_hr_T", FC, ".txt")))
    readcounts_C <- as.matrix(read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\count_hotspotEX_200_FC_hr_C", FC, ".txt")))
    
    info_T <- read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\info_hotspotEX_200_FC_hr_T", FC, ".txt"))
    info_C <- read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\info_hotspotEX_200_FC_hr_C", FC, ".txt"))
    
    rownames(info_T) <- colnames(readcounts_T) <-  paste0("C",1:ncol(readcounts_T))
    rownames(info_C) <-  colnames(readcounts_C) <- paste0("C",1:ncol(readcounts_C))
    
    rownames(readcounts_T) <- rownames(readcounts_C) <-  paste0("G",1:nrow(readcounts_T))
    
    library(spatialDE)
    regdata_T <- SPADE_norm(readcounts=as.matrix(readcounts_T), info=info_T)
    regdata_C <- SPADE_norm(readcounts=as.matrix(readcounts_C), info=info_C)
    res <- SPADE_DE(regdata_T, regdata_C, info_T, info_C, mode="Shape&Strength")
    
    save(res, file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\SPDE_estimate_hr_",cellProp,"_", FC, "EX_ShapeStrength.RData" ))
    
  }
}




cellProp=0.2
FC=2.0
for (cellProp in c(0.1, 0.2, 0.3)){
  for (FC in c(1.5, 2.0, 2.5, 3.0, 3.5)){
    print(paste0("cellProp=", cellProp, "; FC=", FC))
    # the first 100/300 genes are marker genes
    readcounts_T <- as.matrix(read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\count_hotspotEX_200_FC_hr_T", FC, ".txt")))
    readcounts_C <- as.matrix(read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\count_hotspotEX_200_FC_hr_C", FC, ".txt")))
    
    info_T <- read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\info_hotspotEX_200_FC_hr_T", FC, ".txt"))
    info_C <- read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\info_hotspotEX_200_FC_hr_C", FC, ".txt"))
    
    rownames(info_T) <- colnames(readcounts_T) <-  paste0("C",1:ncol(readcounts_T))
    rownames(info_C) <- colnames(readcounts_C) <- paste0("C",1:ncol(readcounts_C))
    
    rownames(readcounts_T) <- rownames(readcounts_C) <-  paste0("G",1:nrow(readcounts_T))
    
    library(spatialDE)
    regdata_T <- SPADE_norm(readcounts=as.matrix(readcounts_T), info=info_T)
    regdata_C <- SPADE_norm(readcounts=as.matrix(readcounts_C), info=info_C)
    res <- SPADE_DE(regdata_T, regdata_C, info_T, info_T, mode="Shape&Strength")

    
    # SPDE_estimate_hr_EX_samecoor <- SPDE_estimate_DE1(readcounts1=regdata_T, readcounts2=regdata_C, location1=info_T,  location2=info_T)
    save(res, file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\SPDE_estimate_hr_",cellProp,"_", FC, "EX.samecoor_ShapeStrength.RData" ))
    
  }
}





library(ROCR)
AUC.all <- NULL
for (cellProp in c(0.1, 0.2, 0.3)){
  for (FC in c(1.5, 2.0, 2.5, 3.0, 3.5)){
    print(paste0("cellProp=", cellProp, "; FC=", FC))
    truth=c(rep(0, 50), rep(1, 150))
    
    SPDE_estimate_hr <- get(load(paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\SPDE_estimate_hr_",cellProp,"_", FC, "EX_ShapeStrength.RData" )))

    
    pred <- prediction(SPDE_estimate_hr$Adjust.Pvalue, truth)
    auc_ROCR <- performance(pred, measure = "auc")
    auc_ROCR_SPDE_ex <- auc_ROCR@y.values[[1]]
    auc_ROCR_SPDE_ex
    
    
    SPDE_estimate_hr <- get(load(paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\SPDE_estimate_hr_",cellProp,"_", FC, "EX.samecoor_ShapeStrength.RData" )))

    pred <- prediction(SPDE_estimate_hr$Adjust.Pvalue, truth)
    auc_ROCR <- performance(pred, measure = "auc")
    auc_ROCR_SPDE_ex_samecoor <- auc_ROCR@y.values[[1]]
    auc_ROCR_SPDE_ex_samecoor
    
    
    AUC.all1 <- data.frame(AUC=c(auc_ROCR_SPDE_ex, auc_ROCR_SPDE_ex_samecoor) , Methods=c("Diff_coor", "Same_coor"))
    AUC.all1$CellProp <- cellProp
    AUC.all1$FC <- FC
    AUC.all <- rbind(AUC.all, AUC.all1)
  }
}    
    




