library(spatstat)

source("utilities.R")

MOB <- t(read.table("Rep11_MOB_count_matrix-1.tsv"))
MOB1 <- MOB[,apply(MOB, 2, sum) >= 10]
MOB2 <- MOB1[(apply(MOB1>0, 1, mean) >= 0.5) ,]
set.seed(2022)
MOB3 <- MOB2[, sample(1:ncol(MOB2), 200)]


library(splatter)
## Simulate pattern-free data
# params <- splatEstimate(MOB3)
# sim_trend <-  splatSimulate(params=params, nGenes=200)
# sim_trend_RC <- counts(sim_trend)
# save(sim_trend_RC, file="simu.data.200Cells.RData")
sim_trend_RC <- get(load("simu.data.200Cells.RData"))
#sim_trend_RC1 <- sim_trend_RC[apply(sim_trend_RC>1, 1, mean) >= 0.3,]
dim(sim_trend_RC1)

numGenes = 200
lowSignal = 25
highSignal = 25
numCell = 200

#### hotspot pattern
for (cellProp in c(0.1, 0.2, 0.3)){
for (ieffect in c(1.2, 1.3, 1.4, 1.5, 1.6)){
  cat("\n ## of Cells",numCell,"\n")
  set.seed(31)
  sim_trend_RC200 <- sim_trend_RC1[sample(1:nrow(sim_trend_RC1), 200),]
  dim(sim_trend_RC200)
 
  set.seed(31)
  pp <- rpoispp(0.00020,  win=owin(c(1,1000),c(1,1000)))
  low_expr    <- c(1)
  high_expr   <- c(2)
  
  pp = add_markdist_hotspot(pp,low_marks=low_expr, high_marks=high_expr, nMarkers=(lowSignal+highSignal), ncell=numGenes, cell_proportion=cellProp)
  
  expr_mat <- matrix(NA,nrow=pp$n,ncol=numGenes)
  for(igene in 1:numGenes){
    back_ground_cells   <- which(pp$markformat[,igene]==1)
    spiked_in_cells     <- which(pp$markformat[,igene]==2)
    
    set.seed(igene)
    a <- sim_trend_RC200[igene,]
    if(igene <= highSignal){
      print("1")
      set.seed(igene)
      a[spiked_in_cells]  <- sim_trend_RC200[igene,spiked_in_cells]*ieffect
    }else if(igene <=lowSignal+highSignal) {
      print("2")
      set.seed(igene)
      a[spiked_in_cells]  <- sim_trend_RC200[igene,spiked_in_cells]/ieffect
    }else if(igene %% 2 == 0) {
      print("3")
      set.seed(igene)
      ram_idx <- sample(1:numCell, round(numCell*cellProp))
      a[ram_idx]  <- sim_trend_RC200[igene,ram_idx]*ieffect
    }else if(igene %% 2 == 1) {
      print("4")
      set.seed(igene)
      ram_idx <- sample(1:numCell, round(numCell*cellProp))
      a[ram_idx]  <- sim_trend_RC200[igene,ram_idx]/ieffect
    } 
    expr_mat[,igene] <- a
  }

  colnames(expr_mat)  <- paste0("gene",1:numGenes)
  info                <- cbind.data.frame(x=pp$x,y=pp$y)
  rownames(info)      <- rownames(expr_mat) <- paste0("C",1:nrow(expr_mat))
  info <- round(info)
  count=t(expr_mat)
  
  write.table(count,file=paste0("../simulation/hotspot/Cell_prop_",cellProp,"/count_hotspot_200_FC", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")
  write.table(info,file=paste0("../simulation/hotspot/Cell_prop_",cellProp,"/info_hotspot_200_FC", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")

  rownames(count) <- paste0("gene", 1:nrow(info))
  colnames(count) <- paste0("cell", 1:ncol(count))
  rownames(info) <- paste0("cell", 1:nrow(info))
  
  write.table(count,file=paste0("../simulation/hotspot/Cell_prop_",cellProp,"/count_hotspot_200_FC", as.character(ieffect), ".spa.txt"), quote = F)
  write.table(info,file=paste0("../simulation/hotspot/Cell_prop_",cellProp,"/info_hotspot_200_FC", as.character(ieffect), ".spa.txt"),col.names = F, quote = F)
  
  count1 <- t(count)
  write.csv(count1, file=paste0("../simulation/hotspot/Cell_prop_",cellProp,"/count_hotspot_200_FC", as.character(ieffect), ".csv"))
  write.csv(info, file=paste0("../simulation/hotspot/Cell_prop_",cellProp,"/info_hotspot_200_FC", as.character(ieffect), ".csv"))
  
}
}









#### Streak pattern
for (cellProp in c(0.1, 0.2, 0.3)){
for (ieffect in c(1.2, 1.3, 1.4, 1.5, 1.6)){
  cat("\n ## of Cells",numCell,"\n")
  set.seed(31)
  sim_trend_RC200 <- sim_trend_RC1[sample(1:nrow(sim_trend_RC1), 200),]
  dim(sim_trend_RC200)
 
  set.seed(31)
  pp <- rpoispp(0.00020,  win=owin(c(1,1000),c(1,1000)))
  low_expr    <- c(1)
  high_expr   <- c(2)
  
  pp = add_markdist_streak(pp,low_marks=low_expr, high_marks=high_expr, nMarkers=(lowSignal+highSignal), ncell=numGenes, cell_proportion=cellProp)
  
  expr_mat <- matrix(NA,nrow=pp$n,ncol=numGenes)
  for(igene in 1:numGenes){
    back_ground_cells   <- which(pp$markformat[,igene]==1)
    spiked_in_cells     <- which(pp$markformat[,igene]==2)
    
    set.seed(igene)
    a <- sim_trend_RC200[igene,]
    if(igene <= highSignal){
      print("1")
      set.seed(igene)
      a[spiked_in_cells]  <- sim_trend_RC200[igene,spiked_in_cells]*ieffect
    }else if(igene <=lowSignal+highSignal) {
      print("2")
      set.seed(igene)
      a[spiked_in_cells]  <- sim_trend_RC200[igene,spiked_in_cells]/ieffect
    }else if(igene %% 2 == 0) {
      print("3")
      set.seed(igene)
      ram_idx <- sample(1:numCell, round(numCell*cellProp))
      a[ram_idx]  <- sim_trend_RC200[igene,ram_idx]*ieffect
    }else if(igene %% 2 == 1) {
      print("4")
      set.seed(igene)
      ram_idx <- sample(1:numCell, round(numCell*cellProp))
      a[ram_idx]  <- sim_trend_RC200[igene,ram_idx]/ieffect
    } 
    expr_mat[,igene] <- a
  }

  colnames(expr_mat)  <- paste0("gene",1:numGenes)
  info                <- cbind.data.frame(x=pp$x,y=pp$y)
  rownames(info)      <- rownames(expr_mat) <- paste0("C",1:nrow(expr_mat))
  info <- round(info)
  count=t(expr_mat)
  
  write.table(count,file=paste0("../simulation/streak/Cell_prop_",cellProp,"/count_streak_200_FC", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")
  write.table(info,file=paste0("../simulation/streak/Cell_prop_",cellProp,"/info_streak_200_FC", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")

  rownames(count) <- paste0("gene", 1:nrow(info))
  colnames(count) <- paste0("cell", 1:ncol(count))
  rownames(info) <- paste0("cell", 1:nrow(info))
  
  write.table(count,file=paste0("../simulation/streak/Cell_prop_",cellProp,"/count_streak_200_FC", as.character(ieffect), ".spa.txt"), quote = F)
  write.table(info,file=paste0("../simulation/streak/Cell_prop_",cellProp,"/info_streak_200_FC", as.character(ieffect), ".spa.txt"),col.names = F, quote = F)
  

  count1 <- t(count)
  write.csv(count1, file=paste0("../simulation/streak/Cell_prop_",cellProp,"/count_streak_200_FC", as.character(ieffect), ".csv"))
  write.csv(info, file=paste0("../simulation/streak/Cell_prop_",cellProp,"/info_streak_200_FC", as.character(ieffect), ".csv"))
  
}
}






