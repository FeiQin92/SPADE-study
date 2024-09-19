## QQ Plots
##==============
library(tidyr)
library(ggplot2)
library(reshape2)
library(scales) 
qplot_gg <- function(pvalues,legend.position="top",ylab=NULL,xlab=NULL,ax.txt.size=20,ax.title.size=20,len.title.size=1.5,len.txt.size=1,len.title=NULL,col.alpha=1,pt.size=2,
    col.base=NULL,cl=1,single_col=NULL,keeporder=T,factor_level=NULL,self_label=NULL){

    if(is.null(col.base)){
        if(cl==1){
            # col.base <- c("hotpink","mediumorchid2","#F4A582","#B8E186","#8dd3c7","wheat")
            col.base <- c("hotpink","mediumorchid2","lightskyblue")
        }else{
            if(is.list(pvalues)){
                col.base <- hue_pal()(length(pvalues))
            }else{
                col.base <- hue_pal()(10)
            } 
        }
    }

    if(is.null(ylab)){ylab=expression(paste("Observed ",-log[10],"(",italic(p),"-value)","\n"))}
    if(is.null(xlab)){xlab=expression(paste("Expected ",-log[10],"(",italic(p),"-value)"))}
    if(is.list(pvalues)){
        nn          <- length(pvalues) 
        exp.vec     <- c()
        for(iset in 1:nn){
            n           <- length(pvalues[[iset]])
            exp.x       <- -log10((rank(pvalues[[iset]], ties.method="first")-.5)/n)
            exp.vec     <- c(exp.vec,exp.x)
        }

        pd              <- melt(pvalues)
        pd$pve          <- exp.vec
        pd$pvo          <- -log10(pd$value)

        if(keeporder){
            if(is.null(factor_level)){
                factor_level <- unique(pd$L1)
            }
            pd$L1       <- factor(pd$L1,levels=factor_level)
        }

        n = length(unique(exp.vec))
        df1 <- data.frame(
                expected = -log10(ppoints(n)),
                clower   = -log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:n, shape2 = n:1)),
                cupper   = -log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:n, shape2 = n:1))
                )
        p1 <- ggplot(df1) +geom_ribbon(aes(x=expected, ymax=cupper, ymin=clower), fill="darkgrey", alpha=.5)+
                geom_abline(intercept = 0, slope = 1, size = 1.5,col="white")
        p2 <- p1 + geom_point(data=pd, aes(pve, pvo,color = L1), size = pt.size,alpha = col.alpha) +
                scale_color_manual(values=col.base,labels = self_label)+
                scale_x_continuous(xlab) +
                scale_y_continuous(ylab) +
                labs(color=len.title)+
                theme_bw()+
                theme(axis.text.x = element_text(size=ax.txt.size),axis.text.y=element_text(size=ax.txt.size),axis.title.y = element_text(size=ax.title.size), axis.title.x = element_text(size=ax.title.size),
                        legend.text= element_text(size=rel(len.txt.size)),legend.title = element_text(size=rel(len.title.size)),
                        panel.border = element_blank(), panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position=legend.position,legend.text.align = 0)
    
    }else{
        n   <- length(pvalues)
        df1 <- data.frame(
            expected = -log10(ppoints(n)),
            observed = -log10(sort(pvalues)),
            clower   = -log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:n, shape2 = n:1)),
            cupper   = -log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:n, shape2 = n:1))
        )

        if(is.null(single_col)){single_col = "lightpink"}

        p1  <- ggplot(df1) +geom_ribbon(aes(x=expected, ymax=cupper, ymin=clower), fill="darkgrey", alpha=.5)+
                geom_abline(intercept = 0, slope = 1, size = 1.5,col="white")+
                geom_point(aes(expected, observed), size = 3,col=single_col,alpha = col.alpha) 

        p2  <- p1 + scale_x_continuous(xlab) +
                scale_y_continuous(ylab) +
                theme_bw()+
                theme(axis.text.x = element_text(size=ax.txt.size),axis.text.y=element_text(size=ax.txt.size),
                        axis.title.y = element_text(size=ax.title.size), axis.title.x = element_text(size=ax.title.size),
                        legend.text= element_text(size=rel(len.txt.size)),legend.title = element_text(size=rel(len.title.size)),
                        panel.border = element_blank(), panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position=legend.position)
    }
    return(p2)
}




## Calculate Empirical FDR
##============================

fdr.cal <- function(altpval,nullpval,rn=NULL){ 
    if(!is.null(rn)){
            order.rn <- rn[order(altpval)]
        }else{
            order.rn <- NULL
        }
    altpval.order   <- altpval[order(altpval)]
    numperm         <- ncol(nullpval)
    nullpval.order  <- sapply(1:numperm,function(x){nullpval[,x][order(nullpval[,x])]})
    out             <- c()
    for(i in 1:nrow(nullpval)){
        s = sum(nullpval.order[1:i,]<=altpval.order[i])
        e   <- s/(i*numperm)
        out <- c(out,e)
    }
    if(is.null(order.rn)){
        return(out)
    }else{
        res <- cbind.data.frame(rn=order.rn,fdr=out)
        return(res)
    }
    
}



## Simulate data with hotspot
##============================
add_markdist_hotspot <- function(pp,low_marks,high_marks,nMarkers,ncell,cell_proportion = 0.2){
    x       <- pp[["x"]]
    y       <- pp[["y"]]
    x_max   <- pp[["window"]][["xrange"]][2]
    y_max   <- pp[["window"]][["yrange"]][2]

    nmarks  <- length(low_marks)
    npoints <- pp[['n']]
    
    # x_half = x_max * 0.5
    # y_half = y_max * 0.5
    set.seed(2021)
    x_middle = matrix(x_max * runif(ncell, 0.1, 0.9), ncell, 1)
    set.seed(2022)
    y_middle = matrix(y_max * runif(ncell, 0.1, 0.9), ncell, 1)
    
    marx = as.data.frame(matrix(NA, nrow = npoints, ncol = nmarks))
    for (j_gene in 1:ncell) {
      marx[, j_gene] <- low_marks
      if (j_gene <= nMarkers){
        # for (i in 1:50){
          set.seed(j_gene)
        #   message(runif(1, cell_proportion-0.15, cell_proportion+0.15))
        # }
  
        center_point <- which.min(sapply(1:npoints,function(p){(x[p]- x_middle[j_gene, 1])^2+(y[p]-y_middle[j_gene, 1])^2}))
        
        cell_proportion1  <- runif(1, cell_proportion/2, cell_proportion*3/2)
        all_dist <- sapply(1:npoints,function(p){(x[p]-x[center_point])^2+(y[p]-y[center_point])^2})
        all_dist_order <- all_dist[order(all_dist)]
        
        high_ind <- which(all_dist %in% all_dist_order[1:round(npoints*cell_proportion1)])
        low_ind = setdiff(1:npoints, high_ind)
        
        marx[low_ind, j_gene]   = low_marks
        marx[high_ind, j_gene]  = high_marks
      }  
    }
    pp[["markformat"]] = marx
    
    return(pp)
}




## Specify hotspot in different locations
add_markdist_hotspot2 <- function(pp,low_marks,high_marks,nMarkers,ncell,cell_proportion = 0.2, seed){
  x       <- pp[["x"]]
  y       <- pp[["y"]]
  x_max   <- pp[["window"]][["xrange"]][2]
  y_max   <- pp[["window"]][["yrange"]][2]
  
  nmarks  <- length(low_marks)
  npoints <- pp[['n']]
  
  # x_half = x_max * 0.5
  # y_half = y_max * 0.5
  set.seed(seed)
  x_middle = matrix(x_max * runif(ncell, 0.1, 0.9), ncell, 1)
  set.seed(seed+1)
  y_middle = matrix(y_max * runif(ncell, 0.1, 0.9), ncell, 1)
  
  marx = as.data.frame(matrix(NA, nrow = npoints, ncol = nmarks))
  for (j_gene in 1:ncell) {
    marx[, j_gene] <- low_marks
    if (j_gene <= nMarkers){
      # for (i in 1:50){
      set.seed(j_gene)
      #   message(runif(1, cell_proportion-0.15, cell_proportion+0.15))
      # }
      
      center_point <- which.min(sapply(1:npoints,function(p){(x[p]- x_middle[j_gene, 1])^2+(y[p]-y_middle[j_gene, 1])^2}))
      
      cell_proportion1  <- runif(1, cell_proportion/2, cell_proportion*3/2)
      all_dist <- sapply(1:npoints,function(p){(x[p]-x[center_point])^2+(y[p]-y[center_point])^2})
      all_dist_order <- all_dist[order(all_dist)]
      
      high_ind <- which(all_dist %in% all_dist_order[1:round(npoints*cell_proportion1)])
      low_ind = setdiff(1:npoints, high_ind)
      
      marx[low_ind, j_gene]   = low_marks
      marx[high_ind, j_gene]  = high_marks
    }  
  }
  pp[["markformat"]] = marx
  
  return(pp)
}




## Simulate data with two hotspots
##============================
add_markdist_bi_hotspot <- function(pp,low_marks,high_marks,nMarkers,ncell,cell_proportion = 0.2){
  x       <- pp[["x"]]
  y       <- pp[["y"]]
  x_max   <- pp[["window"]][["xrange"]][2]
  y_max   <- pp[["window"]][["yrange"]][2]
  
  nmarks  <- length(low_marks)
  npoints <- pp[['n']]
  
  set.seed(2021)
  x_middle = matrix(x_max * runif(ncell*2, 0.1, 0.9), ncell, 2)
  set.seed(2022)
  y_middle = matrix(y_max * runif(ncell*2, 0.1, 0.9), ncell, 2)
  
  # x_half = x_max * 0.5
  # y_half = y_max * 0.5
  
  # center_point <- which.min(sapply(1:npoints,function(p){(x[p]-x_half)^2+(y[p]-y_half)^2}))
  # 
  # all_dist <- sapply(1:npoints,function(p){(x[p]-x[center_point])^2+(y[p]-y[center_point])^2})
  # all_dist_order <- all_dist[order(all_dist)]
  # 
  
  marx = as.data.frame(matrix(NA, nrow = npoints, ncol = nmarks))
  for (j_gene in 1:ncell) {
    marx[, j_gene] <- low_marks
    if (j_gene <= nMarkers){
      # for (i in 1:50){
      set.seed(j_gene)
      #   message(runif(1, cell_proportion-0.15, cell_proportion+0.15))
      # }
      center_point1 <- which.min(sapply(1:npoints,function(p){(x[p]- x_middle[j_gene, 1])^2+(y[p]-y_middle[j_gene, 1])^2}))
      center_point2 <- which.min(sapply(1:npoints,function(p){(x[p]- x_middle[j_gene, 2])^2+(y[p]-y_middle[j_gene, 2])^2}))
      
      cell_proportion1  <- runif(1, cell_proportion/2, cell_proportion*3/2)
      
      all_dist1 <- sapply(1:npoints,function(p){(x[p]-x[center_point1])^2+(y[p]-y[center_point1])^2})
      all_dist_order1 <- all_dist1[order(all_dist1)]
      
      high_ind1 <- which(all_dist1 %in% all_dist_order1[1:round(npoints*cell_proportion1)])
      # low_ind1 = setdiff(1:npoints, high_ind1)
      
      
      all_dist2 <- sapply(1:npoints,function(p){(x[p]-x[center_point2])^2+(y[p]-y[center_point2])^2})
      all_dist_order2 <- all_dist2[order(all_dist2)]
      
      high_ind2 <- which(all_dist2 %in% all_dist_order2[1:round(npoints*cell_proportion1)])
      # low_ind2 = setdiff(1:npoints, high_ind2)
      
      # marx[low_ind1, j_gene]   = low_marks
      marx[high_ind1, j_gene]  = high_marks
      
      # marx[low_ind2, j_gene]   = low_marks
      marx[high_ind2, j_gene]  = high_marks
    }  
  }
  pp[["markformat"]] = marx
  
  return(pp)
}







## Simulate data with two hotspots
##============================
add_markdist_tri_hotspot <- function(pp,low_marks,high_marks,nMarkers,ncell,cell_proportion = 0.2){
  x       <- pp[["x"]]
  y       <- pp[["y"]]
  x_max   <- pp[["window"]][["xrange"]][2]
  y_max   <- pp[["window"]][["yrange"]][2]
  
  nmarks  <- length(low_marks)
  npoints <- pp[['n']]
  
  set.seed(2021)
  x_middle = matrix(rep(rep(x_max * c(0.25, 0.75), 2), each=ncell), ncell, 4)
  set.seed(2022)
  y_middle = matrix(rep(rep(x_max * c(0.25, 0.75), each=2), each=ncell), ncell, 4)
  
  marx = as.data.frame(matrix(NA, nrow = npoints, ncol = nmarks))
  for (j_gene in 1:ncell) {
    marx[, j_gene] <- low_marks
    if (j_gene <= nMarkers){
      # for (i in 1:50){
      set.seed(j_gene)
      #   message(runif(1, cell_proportion-0.15, cell_proportion+0.15))
      # }
      cell_proportion1  <- runif(1, cell_proportion/2, cell_proportion*3/2)
      
      for (r in c(1:ncol(x_middle))){
        center_point1 <- which.min(sapply(1:npoints,function(p){(x[p]- x_middle[j_gene, r])^2+(y[p]-y_middle[j_gene, r])^2}))
        all_dist1 <- sapply(1:npoints,function(p){(x[p]-x[center_point1])^2+(y[p]-y[center_point1])^2})
        all_dist_order1 <- all_dist1[order(all_dist1)]
        
        high_ind1 <- which(all_dist1 %in% all_dist_order1[1:round(npoints*cell_proportion1)])
        marx[high_ind1, j_gene]  = high_marks
      }
    }  
  }
  pp[["markformat"]] = marx
  
  return(pp)
}






## Simulate data with streak
##============================
add_markdist_streak <- function(pp,low_marks,high_marks,nMarkers,ncell,cell_proportion = 0.2){
    x       <- pp[["x"]]
    y       <- pp[["y"]]
    x_max   <- pp[["window"]][["xrange"]][2]
    y_max   <- pp[["window"]][["yrange"]][2]

    nmarks  <- length(low_marks)
    npoints <- pp[['n']]

    set.seed(2021)
    x_middle = matrix(x_max * runif(ncell*1, 0.15, 0.85), ncell, 1)
    set.seed(2022)
    y_middle = matrix(y_max * runif(ncell*1, 0.15, 0.85), ncell, 1)
    
    # border_point  <- which.max(x)
    # all_dist <- sapply(1:npoints,function(p){(x[p]-x[center_point])^2+(y[p]-y[center_point])^2})
    # all_dist_order <- all_dist[order(all_dist)]

    marx = as.data.frame(matrix(NA, nrow = npoints, ncol = nmarks))
    for (j_gene in 1:ncell) {
      marx[, j_gene] <- low_marks
      if (j_gene <= nMarkers){
        # for (i in 1:50){
        set.seed(j_gene)
        #   message(runif(1, cell_proportion-0.15, cell_proportion+0.15))
        #  }
        center_point <- which.min(sapply(1:npoints,function(p){(x[p]- x_middle[j_gene, 1])^2+(y[p]-y_middle[j_gene, 1])^2}))
        cell_proportion1  <- runif(1, cell_proportion/2, cell_proportion*3/2)
        
        x_dist          <- abs(x[center_point] - x)
        x_dist_order    <- x_dist[order(x_dist)]
        
        high_ind <- which(x_dist %in% x_dist_order[1:round(npoints*cell_proportion1)])
        low_ind = setdiff(1:npoints, high_ind)
      
        marx[low_ind, j_gene]   = low_marks
        marx[high_ind, j_gene]  = high_marks
      }
    }
    pp[["markformat"]] = marx
    
    return(pp)
}




add_markdist_streak2 <- function(pp,low_marks,high_marks,nMarkers,ncell,cell_proportion = 0.2, seed){
  x       <- pp[["x"]]
  y       <- pp[["y"]]
  x_max   <- pp[["window"]][["xrange"]][2]
  y_max   <- pp[["window"]][["yrange"]][2]
  
  nmarks  <- length(low_marks)
  npoints <- pp[['n']]
  
  set.seed(seed)
  x_middle = matrix(x_max * runif(ncell*1, 0.15, 0.85), ncell, 1)
  set.seed(seed+1)
  y_middle = matrix(y_max * runif(ncell*1, 0.15, 0.85), ncell, 1)
  
  # border_point  <- which.max(x)
  # all_dist <- sapply(1:npoints,function(p){(x[p]-x[center_point])^2+(y[p]-y[center_point])^2})
  # all_dist_order <- all_dist[order(all_dist)]
  
  marx = as.data.frame(matrix(NA, nrow = npoints, ncol = nmarks))
  for (j_gene in 1:ncell) {
    marx[, j_gene] <- low_marks
    if (j_gene <= nMarkers){
      # for (i in 1:50){
      set.seed(j_gene)
      #   message(runif(1, cell_proportion-0.15, cell_proportion+0.15))
      #  }
      center_point <- which.min(sapply(1:npoints,function(p){(x[p]- x_middle[j_gene, 1])^2+(y[p]-y_middle[j_gene, 1])^2}))
      cell_proportion1  <- runif(1, cell_proportion/2, cell_proportion*3/2)
      
      x_dist          <- abs(x[center_point] - x)
      x_dist_order    <- x_dist[order(x_dist)]
      
      high_ind <- which(x_dist %in% x_dist_order[1:round(npoints*cell_proportion1)])
      low_ind = setdiff(1:npoints, high_ind)
      
      marx[low_ind, j_gene]   = low_marks
      marx[high_ind, j_gene]  = high_marks
    }
  }
  pp[["markformat"]] = marx
  
  return(pp)
}



add_markdist_streak_h <- function(pp,low_marks,high_marks,nMarkers,ncell,cell_proportion = 0.2){
  x       <- pp[["x"]]
  y       <- pp[["y"]]
  x_max   <- pp[["window"]][["xrange"]][2]
  y_max   <- pp[["window"]][["yrange"]][2]
  
  nmarks  <- length(low_marks)
  npoints <- pp[['n']]
  
  set.seed(2021)
  x_middle = matrix(x_max * runif(ncell*1, 0.15, 0.85), ncell, 1)
  set.seed(2022)
  y_middle = matrix(y_max * runif(ncell*1, 0.15, 0.85), ncell, 1)
  
  # border_point  <- which.max(x)
  # all_dist <- sapply(1:npoints,function(p){(x[p]-x[center_point])^2+(y[p]-y[center_point])^2})
  # all_dist_order <- all_dist[order(all_dist)]
  
  marx = as.data.frame(matrix(NA, nrow = npoints, ncol = nmarks))
  for (j_gene in 1:ncell) {
    marx[, j_gene] <- low_marks
    if (j_gene <= nMarkers){
      # for (i in 1:50){
      set.seed(j_gene)
      #   message(runif(1, cell_proportion-0.15, cell_proportion+0.15))
      #  }
      center_point <- which.min(sapply(1:npoints,function(p){(x[p]- x_middle[j_gene, 1])^2+(y[p]-y_middle[j_gene, 1])^2}))
      cell_proportion1  <- runif(1, cell_proportion/2, cell_proportion*3/2)
      
      y_dist          <- abs(y[center_point] - y)
      y_dist_order    <- y_dist[order(y_dist)]
      
      high_ind <- which(y_dist %in% y_dist_order[1:round(npoints*cell_proportion1)])
      low_ind = setdiff(1:npoints, high_ind)
      
      marx[low_ind, j_gene]   = low_marks
      marx[high_ind, j_gene]  = high_marks
    }
  }
  pp[["markformat"]] = marx
  
  return(pp)
}




## Simulate data with streak
##============================
add_markdist_bi_streak <- function(pp,low_marks,high_marks,nMarkers,ncell,cell_proportion = 0.2){
  x       <- pp[["x"]]
  y       <- pp[["y"]]
  x_max   <- pp[["window"]][["xrange"]][2]
  y_max   <- pp[["window"]][["yrange"]][2]
  
  nmarks  <- length(low_marks)
  npoints <- pp[['n']]
  
  set.seed(2021)
  x_middle = matrix(x_max * runif(ncell*2, 0.1, 0.9), ncell, 2)
  set.seed(2022)
  y_middle = matrix(y_max * runif(ncell*2, 0.1, 0.9), ncell, 2)
  
  marx = as.data.frame(matrix(NA, nrow = npoints, ncol = nmarks))
  for (j_gene in 1:ncell) {
    marx[, j_gene] <- low_marks
    if (j_gene <= nMarkers){
      
      set.seed(j_gene)
      center_point1 <- which.min(sapply(1:npoints,function(p){(x[p]- x_middle[j_gene, 1])^2+(y[p]-y_middle[j_gene, 1])^2}))
      center_point2 <- which.min(sapply(1:npoints,function(p){(x[p]- x_middle[j_gene, 2])^2+(y[p]-y_middle[j_gene, 2])^2}))
      
      cell_proportion1  <- runif(1, cell_proportion/2, cell_proportion*3/2)
      
      x_dist1 <- abs(x[center_point1] - x)
      x_dist_order1 <- x_dist1[order(x_dist1)]
      
      high_ind1 <- which(x_dist1 %in% x_dist_order1[1:round(npoints*cell_proportion1)])
      
      x_dist2 <- abs(x[center_point2] - x)
      x_dist_order2 <- x_dist2[order(x_dist2)]
      
      high_ind2 <- which(x_dist2 %in% x_dist_order2[1:round(npoints*cell_proportion1)])
      
      marx[high_ind1, j_gene]  = high_marks
      marx[high_ind2, j_gene]  = high_marks
    }
  }
  pp[["markformat"]] = marx
  
  return(pp)
}






## Simulate data with streak
##============================
add_markdist_tri_streak <- function(pp,low_marks,high_marks,nMarkers,ncell,cell_proportion = 0.2){
  x       <- pp[["x"]]
  y       <- pp[["y"]]
  x_max   <- pp[["window"]][["xrange"]][2]
  y_max   <- pp[["window"]][["yrange"]][2]
  
  nmarks  <- length(low_marks)
  npoints <- pp[['n']]
  
  set.seed(2021)
  x_middle = matrix(rep(x_max * c(1:4) * 0.2, each=ncell), ncell, 4)
  set.seed(2022)
  y_middle = matrix(rep(x_max * c(1:4) * 0.2, each=ncell), ncell, 4)
  
  marx = as.data.frame(matrix(NA, nrow = npoints, ncol = nmarks))
  for (j_gene in 1:ncell) {
    marx[, j_gene] <- low_marks
    if (j_gene <= nMarkers){
      
      set.seed(j_gene)
      cell_proportion1  <- runif(1, cell_proportion/2, cell_proportion*3/2)

      for (r in 1:ncol(x_middle)){
        center_point1 <- which.min(sapply(1:npoints,function(p){(x[p]- x_middle[j_gene, r])^2+(y[p]-y_middle[j_gene, r])^2}))
        x_dist1 <- abs(x[center_point1] - x)
        x_dist_order1 <- x_dist1[order(x_dist1)]
        
        high_ind1 <- which(x_dist1 %in% x_dist_order1[1:round(npoints*cell_proportion1)])
        marx[high_ind1, j_gene]  = high_marks
      }
    }
  }
  pp[["markformat"]] = marx
  
  return(pp)
}






## Simulate data with gradient
##============================
add_markdist_lin <- function(expr,ranidx,decrease=F){
    out             <- expr
    sub_expr        <- expr[ranidx]
    order_sub_expr  <- sub_expr[order(sub_expr,decreasing=decrease)] 
    out[ranidx[order(ranidx)]] <- order_sub_expr
    return(out)
}



## Simulate data with gradient
##============================
pattern_plot_sparkx <- function(pltdat,igene,xy=T,main=F,titlesize=2,pointsize=3,min.pand=0.99,max.pand=1.01,title=NULL,pal=NULL,expand_par=0.05,ncolors=5){
    if(!xy){
        xy              <- matrix(as.numeric(do.call(rbind,strsplit(as.character(pltdat[,1]),split="x"))),ncol=2)
        rownames(xy)    <- as.character(pltdat[,1])
        colnames(xy)    <- c("x","y")
        pd              <- cbind.data.frame(xy,pltdat[,2:ncol(pltdat)])
    }else{
        pd              <- pltdat
    }
    # pal             <- colorRampPalette(c("mediumseagreen","lightyellow2","deeppink"))
    library(viridis)
    if(is.null(pal)){
       # pal             <- colorRampPalette(c("antiquewhite",viridis_pal(alpha=0.8)(10)[c(10,6,2)]))

      pal             <- colorRampPalette(c("antiquewhite",viridis_pal()(10)[c(6,5,4,2)]))

      # pal             <- colorRampPalette(c("antiquewhite",viridis_pal()(10)[c(6,2)]))
    }
    
    pd1 <- pd[,c(1,2,igene+2)]
    colnames(pd1)[3] <- "RE"
    gpt             <- ggplot(pd1,aes(x=x, y=y,color=RE)) + 
                          geom_point(size=pointsize) +
                        scale_color_gradientn(colours=pal(ncolors))+ 
                        # viridis::scale_color_viridis(alpha=0.8,option = opt, direction = direct)+
                        scale_x_discrete(expand = c(0, expand_par))+ 
                        scale_y_discrete(expand = c(0, expand_par))+
                        # guides(fill = guide_legend(title="Relative Expression"))  +        
                        # scale_fill_manual("Relative Expression")+
                        expand_limits(x=c(min(pd$x)*min.pand,max(pd$x)*max.pand),y=c(min(pd$y)*min.pand,max(pd$y)*max.pand))+
                        theme_bw()
    if(main){
        if(is.null(title)){
          title=colnames(pd)[igene+2]
        }
        out = gpt + labs(title = title, x = NULL, y = NULL)+
              theme(plot.title = element_text(hjust = 0.5,size=rel(titlesize),face="italic"))
      }else{
        out = gpt + labs(title = NULL, x = NULL, y = NULL)
    }
    return(out)
}




relative_func <- function(expres){
    maxd    = max(expres)-min(expres)
    rexpr   = (expres-min(expres))/maxd
    return(rexpr)
}

