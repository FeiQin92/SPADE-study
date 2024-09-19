library("dplyr")
library("gstat")
library("sp")
library("spacetime")
# library("STRbook")
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
    pp <- rpoispp(0.00020,  win=owin(c(1,1000),c(1,1000)))
    pp_T = add_markdist_hotspot(pp,low_marks=low_expr, high_marks=high_expr, nMarkers=(lowSignal+highSignal), ncell=numGenes, 
                                cell_proportion=cellProp)
    

    ## Specify hotspot in different locations
    set.seed(69)
    pp_C = add_markdist_hotspot2(pp,low_marks=low_expr, high_marks=high_expr, nMarkers=(lowSignal+highSignal), ncell=numGenes, 
                                 cell_proportion=cellProp, seed=69)

    
    expr_mat_C <- matrix(NA,nrow=pp_C$n,ncol=numGenes)
    expr_mat_T <- matrix(NA,nrow=pp_T$n,ncol=numGenes)
    expr_mat_T2 <- matrix(NA,nrow=pp_T$n,ncol=numGenes)
    for(igene in 1:numGenes){
      
      spiked_in_cells_T   <- which(pp_T$markformat[,igene]==2)
      spiked_in_cells_C   <- which(pp_C$markformat[,igene]==2)
      
      a <- sim_trend_RC400[igene,]
      if(igene <= highSignal){
        a[spiked_in_cells_T]  <- sim_trend_RC400[igene,spiked_in_cells_T]*ieffect
      }else if(igene <=lowSignal+highSignal) {
        a[spiked_in_cells_T]  <- sim_trend_RC400[igene,spiked_in_cells_T]/ieffect
      }
      
      #Add hotspots with the different hotspots locations#
      b <- sim_trend_RC400[igene,]
      if(igene <= highSignal){
        b[spiked_in_cells_C]  <- sim_trend_RC400[igene,spiked_in_cells_C]*ieffect
      }else if(igene <=lowSignal+highSignal) {
        b[spiked_in_cells_C]  <- sim_trend_RC400[igene,spiked_in_cells_C]/ieffect
      }
      
      #Add hotspots with the same hotspots location#
      b2 <- sim_trend_RC400[igene,]
      if(igene <= highSignal){
        b2[spiked_in_cells_T]  <- sim_trend_RC400[igene,spiked_in_cells_T]*ieffect
      }else if(igene <=lowSignal+highSignal) {
        b2[spiked_in_cells_T]  <- sim_trend_RC400[igene,spiked_in_cells_T]/ieffect
      }
      
      expr_mat_T[,igene] <- a
      expr_mat_C[,igene] <- b
      expr_mat_T2[,igene] <- b2
      
    }
    
    colnames(expr_mat_T) <- colnames(expr_mat_C)  <- paste0("gene",1:numGenes)
    info_T                <- cbind.data.frame(x=pp_T$x, y=pp_T$y)
    info_T <- round(info_T)
    
    # expr_mat_C <- expr_mat_C[sample(ncol(expr_mat_C), nrow(info_C)),]
    
    rownames(info_T) <-  rownames(expr_mat_T) <- paste0("C",1:nrow(expr_mat_T))
    rownames(info_C) <- rownames(expr_mat_C) <- paste0("C",1:nrow(expr_mat_C))
    count_T=t(expr_mat_T)
    count_T2=t(expr_mat_T2)
    count_C=t(expr_mat_C)
    
    write.table(pp_T$markformat,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\count_hotspotEX_200_FC_hr_indexT", as.character(ieffect), "_T2.txt"),row.names = F,col.names = F,quote = F,sep="\t")
    write.table(pp_C$markformat,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\count_hotspotEX_200_FC_hr_indexC", as.character(ieffect), "_T2.txt"),row.names = F,col.names = F,quote = F,sep="\t")
    
    write.table(count_T,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\count_hotspotEX_200_FC_hr_T", as.character(ieffect), "_T2.txt"),row.names = F,col.names = F,quote = F,sep="\t")

    # hotspot with the same locations
    write.table(count_T2,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\count_hotspotEX_200_FC_hr_T2", as.character(ieffect), "_T2.txt"),row.names = F,col.names = F,quote = F,sep="\t")

    # hotspot with different locations
    write.table(count_C,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\count_hotspotEX_200_FC_hr_C", as.character(ieffect), "_T2.txt"),row.names = F,col.names = F,quote = F,sep="\t")
    
    write.table(info_T,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\info_hotspotEX_200_FC_hr_T", as.character(ieffect), "_T2.txt"),row.names = F,col.names = F,quote = F,sep="\t")
  }
}






### exchange both length scale theta and gamma values
#########################################################
#########################################################
#########################################################

## hotspot with different directions and same coordinates
cellProp=0.2
FC=2.0
for (cellProp in c(0.1, 0.2, 0.3)){
  for (FC in c(1.5, 2.0, 2.5, 3.0, 3.5)){
    print(paste0("cellProp=", cellProp, "; FC=", FC))
    # all the 200 genes in both groups are genes with markers but with different locations
    readcounts_T <- as.matrix(read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\count_hotspotEX_200_FC_hr_T", FC, "_T2.txt")))
    readcounts_C <- as.matrix(read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\count_hotspotEX_200_FC_hr_C", FC, "_T2.txt")))
    
    info_T <- read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\info_hotspotEX_200_FC_hr_T", FC, "_T2.txt"))
    # info_C <- read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\info_hotspotEX_200_FC_hr_C", FC, "_T2.txt"))
    
    rownames(info) <- colnames(readcounts_T) <-  paste0("C",1:ncol(readcounts_T))
    colnames(readcounts_C) <- paste0("C",1:ncol(readcounts_C))
    
    rownames(readcounts_T) <- rownames(readcounts_C) <-  paste0("G",1:nrow(readcounts_T))
    
    library(spatialDE)
    regdata_T <- SPADE_norm(readcounts=as.matrix(readcounts_T), info=info)
    regdata_C <- SPADE_norm(readcounts=as.matrix(readcounts_C), info=info)
    res <- SPADE_DE(regdata_T, regdata_C, info, info, mode="Shape&Strength")
  
    save(res, file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\SPDE_estimate_hr_",cellProp,"_", FC, "EX.diffDir_ShapeStrength.RData" ))
    
  }
}



cellProp=0.2
FC=2.0
for (cellProp in c(0.1, 0.2, 0.3)){
  for (FC in c(1.5, 2.0, 2.5, 3.0, 3.5)){
    print(paste0("cellProp=", cellProp, "; FC=", FC))
    # all the 200 genes in both groups are genes with markers but with different locations
    readcounts_T <- as.matrix(read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\count_hotspotEX_200_FC_hr_T", FC, "_T2.txt")))
    readcounts_C <- as.matrix(read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\count_hotspotEX_200_FC_hr_T2", FC, "_T2.txt")))
    
    info <- read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\info_hotspotEX_200_FC_hr_T", FC, "_T2.txt"))
    # info_C <- read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\info_hotspotEX_200_FC_hr_C", FC, "_T2.txt"))
    
    rownames(info) <- colnames(readcounts_T) <-  paste0("C",1:ncol(readcounts_T))
                      colnames(readcounts_C) <- paste0("C",1:ncol(readcounts_C))
    
    rownames(readcounts_T) <- rownames(readcounts_C) <-  paste0("G",1:nrow(readcounts_T))
    
    library(spatialDE)
    regdata_T <- SPADE_norm(readcounts=as.matrix(readcounts_T), info=info)
    regdata_C <- SPADE_norm(readcounts=as.matrix(readcounts_C), info=info)
    res <- SPADE_DE(regdata_T, regdata_C, info, info, mode="Shape&Strength")
  
    save(res, file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\SPDE_estimate_hr_",cellProp,"_", FC, "EX.sameDir_ShapeStrength.RData" ))
    
  }
}




library(ROCR)
FPR.all <- NULL
for (cellProp in c(0.1, 0.2, 0.3)){
  for (FC in c(1.5, 2.0, 2.5, 3.0, 3.5)){
    print(paste0("cellProp=", cellProp, "; FC=", FC))
    truth=c(rep(0, 50), rep(1, 150))
    
    SPDE_estimate_hr <- get(load(paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\SPDE_estimate_hr_",cellProp,"_", FC, "EX.diffDir_ShapeStrength.RData" )))  
    FPR_diffDir <- sum(SPDE_estimate_hr$Adjust.Pvalue < 0.05)/200
    FPR_diffDir

    
    SPDE_estimate_hr <- get(load(paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspotEX\\Cell_prop_",cellProp,"\\SPDE_estimate_hr_",cellProp,"_", FC, "EX.sameDir_ShapeStrength.RData" )))
    FPR_sameDir <- sum(SPDE_estimate_hr$Adjust.Pvalue < 0.05)/200
    FPR_sameDir
    
    FPR.all1 <- data.frame(FPR=c(FPR_diffDir, FPR_sameDir) , Methods=c("Diff_dir","Same_dir"))
    FPR.all1$CellProp <- cellProp
    FPR.all1$FC <- FC
    FPR.all <- rbind(FPR.all, FPR.all1)
  }
}    









