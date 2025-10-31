#' @name plot.batss
#' @title Plot function for 'BATSS' outputs
#' @description Plot for objects of class 'batss' 
#' @param x An object of class 'batss' (i.e., output of the function [batss.glm]).
#' @param which An integer vector indicating the list of desired plots. If a subset of the plots is required, specify a subset of the numbers \code{1:4}. By default, all plots are provided, i.e., `which=1:4`. Plot `1` displays the boxplot of the total sample size as well as the boxplot of the sample sizes per group, Plot `2` displays a barplot of the probability of stopping at each look, plot `3` displays the violin plot of the sample size per group, and plot `4` displays the Monte Carlo trial target estimates as a function of the sample size.
#' @param ask A logical. If `TRUE`, the user is prompted to hit the \kbd{Enter} key before each plot. See \code{\link{par}(ask=.)}. The default is \code{ask=TRUE}.
#' @param hypothesis A character string indicating which alternative hypothesis to use for analyses considering both "H0" and "H1", with options "H1" (default) and "H0".
#' @param title Either a \link[base]{logical} indicating if a title should be added or a string (of class \link[base]{character}) indicating the title to be added. If `title` equals `TRUE` (default), the title 'Under 'H1' or 'Under 'H0' (depending on the argument `hypothesis`) is added to the outer margin of the plot. No outer margin space is added if `title = FALSE`.
#' @param legend A \link[base]{logical} (with default set to `TRUE`) if a legend should be added at the bottom of plots 3 and 4, or a list with names `height`, `cex` and `pt` respectively indicating i/ the fraction of the plot to be used for the legend as a numeric (default is `.15`), ii/ the character expansion factor relative to current `par("cex")` as a numeric (default to `1.25`), and iii/ the expansion factor(s) for the points as a numeric (default to `2`). The input `legend = TRUE` is equivalent to `legend = c(height=.15, cex=1, pt=2)`.
#' @param ess A \link[base]{logical} (with default set to `TRUE`) indicating if the expected sample size should be displayed in plots 2, 3 and 4, or a list with names `col`, `cex` and `bg` respectively indicating i/ the colour of the label as a character (plots 2, 3, and 4), ii/ the text expension level as a numerical value (plots 2, 3 and 4) and iii/ the background colour as a character (plot 3). The input `ess = TRUE` is equivalent to `ess = list(col="blue", cex=1, bg="#FFD70070")`.
#' @param percentage A \link[base]{logical} (with default set to `TRUE`) indicating if the probability of stopping at each look should be displayed in plots 2 (as a percentage), or a list with names `col`, and `cex` indicating i/ the colour of the label as a character, ii/ the text expension level as a numerical value. The input `percentage = TRUE` is equivalent to `percentage = list(col="violet", cex=1)`.
#' @param beta A \link[base]{logical} (with default set to `TRUE`) indicating if the assumed target parameter value(s) should be displayed in plot 4, or a list with names `col`, `lwd`, and `lty` respectively indicating i/ the colour of the target parameter value horizontal line(s) as a character, ii/ the thickness of the line(s) as a numerical value  and iii/ the type of line as an integer (see \link[graphics]{par} for details). The input `beta = TRUE` is equivalent to `beta = list(col="blue", lwd=1.5, lty=1)`.
#' @param col A vector of length 4 specifying the colours to be used. Default to `c("#8B897040","#008B0040","#8B3A3A40","#FF990075")` where the last two digits of the hexadecimal strings specify the level of transparency. Refer to the Section 'colour specification' in \link[graphics]{par} for details. If the length of `col` equals 1, the same colour is used for all cases. For plots `1` and `2`, the 3rd colour of the vector `col` is used to display the barplot and boxplots. 
#' @param smooth A numerical (>0) indictating the level of smoothing of the violin plots (for plot `3`). Default to `1`. When `smooth=NULL`, the smoothing value is optimised in the \link[sm]{sm.density} function.
#' @param ... Additional arguments affecting the plot produced, like ylim and ylab.
#' @returns Generates graphical displays of results for objects of class 'batss'.
#' @seealso [batss.glm()], [batss.surv()], the function generating S3 objects of class 'batss'. 
#' @export

plot.batss = function(x, which=if (x$type!="surv") 1:4 else 1:5, ask=TRUE, hypothesis=ifelse(is.null(x$H1),"H0","H1"), 
                      title=TRUE, legend=TRUE, ess=TRUE, percentage=TRUE, beta=TRUE, 
                      col = c("#008B0040","#8B3A3A40","#8B897040","#FF990075"),
                      smooth = 1, ...){
  which.plot=rep(TRUE,5)
  if(!is.null(which)){which.plot[-which]=FALSE}    
  if(sum(which.plot)==1){ask=FALSE}    
  if(ask){oask <- devAskNewPage(TRUE)
  on.exit(devAskNewPage(oask))
  }
  oldpar <- par(no.readonly = TRUE) 
  on.exit(par(oldpar)) 
  mc = match.call()    
  # hypothesis 
  if(is.na(match(hypothesis,names(x)))){
    stop("'",hypothesis,"' not available in 'batss' object")  
  }
  # title
  if(is.character(title)){
    titlew = title
    title  = TRUE
  }else{
    if(is.logical(title)){
      if(title){titlew = paste0("Under '",hypothesis,"'")}
    }else{
      warning("title should be a logical or a character string")
      title = FALSE
    }
  }     
  # col
  if(length(col)==1){col = rep(col,4)
  }else{
    if(length(col)<4){
      col[(length(col)+1):4] = c("#008B0040","#8B3A3A40","#8B897040","#FF990075")[(length(col)+1):4]
    }
  }
  # legend
  if(is.list(legend)){# legend = list(col="green", cex=1.25)
    legend.height = ifelse(!is.null(legend[["height"]]),legend[["height"]],.15)        
    legend.cex    = ifelse(!is.null(legend[["cex"]]),legend[["cex"]],1)    
    legend.pt     = ifelse(!is.null(legend[["pt"]]),legend[["pt"]],2)    
    legend        = TRUE
  }else{
    if(!is.logical(legend)){
      stop("legend should be a list or a logical")
    }else{
      if(legend){
        legend.height = .15      
        legend.cex    = 1
        legend.pt     = 2
      }
    }
  }
  # ess
  if(is.list(ess)){# ess = list(col="green", cex=1.25)
    ess.col = ifelse(!is.null(ess[["col"]]),ess[["col"]],"blue")        
    ess.cex = ifelse(!is.null(ess[["cex"]]),ess[["cex"]],1)    
    ess.bg  = ifelse(!is.null(ess[["bg"]]),ess[["bg"]],"#FFD70070")    
    ess     = TRUE
  }else{
    if(!is.logical(ess)){
      stop("ess should be a list or a logical")
    }else{
      if(ess){
        ess.col = "blue"       
        ess.cex = 1
        ess.bg  = "#FFD70070"
      }
    }
  }
  # percentage
  if(is.list(percentage)){# perc = list(col="green", cex=1.25)
    percentage.col = ifelse(!is.null(percentage[["col"]]),percentage[["col"]],"violet")        
    percentage.cex = ifelse(!is.null(percentage[["cex"]]),percentage[["cex"]],1)    
    percentage     = TRUE
  }else{
    if(!is.logical(percentage)){
      stop("perc should be a list or a logical")
    }else{
      if(percentage){
        percentage.col = "violet"       
        percentage.cex = 1
      }
    }
  }    
  # beta
  if(is.list(beta)){# 
    beta.col = ifelse(!is.null(beta[["col"]]),beta[["col"]],"blue")
    beta.lwd = ifelse(!is.null(beta[["lwd"]]),beta[["lwd"]],1.5)    
    beta.lty = ifelse(!is.null(beta[["lty"]]),beta[["lty"]],1)    
    percentage     = TRUE
  }else{
    if(!is.logical(beta)){
      stop("beta should be a list or a logical")
    }else{
      if(beta){
        beta.col = "blue"       
        beta.lwd = 1.5
        beta.lty = 1
      }
    }
  }        
  # prep
  id.targetw = x[[hypothesis]]$target$par
  n.targetw  = nrow(id.targetw)
  if (x$type=="surv") {
    id.targetw$beta = x$hr[id.targetw$id,grepl(hypothesis,colnames(x$hr))]
  } else {
    id.targetw$beta = x$beta[id.targetw$id,grepl(hypothesis,colnames(x$beta))]
  }
  ar.inf.rt4 = array(dim=c(dim(x[[hypothesis]]$sample)[1],n.targetw,4),
                     dimnames=list(1:dim(x[[hypothesis]]$sample)[1],
                                   x[[hypothesis]]$target$id,
                                   c("n","est","type","last")))
  ar.inf.rt4[,,"n"]    = if (x$type=="surv") as.matrix(x[[hypothesis]]$sample[,paste0("n(",x[[hypothesis]]$target$par$group,")")]) else as.matrix(x[[hypothesis]]$sample[,x[[hypothesis]]$target$par$group])
  ar.inf.rt4[,,"est"]  = t(x[[hypothesis]]$estimate[,"mid",])
  ar.inf.rt4[,,"type"] = t(x[[hypothesis]]$estimate[,"type",])
  ar.inf.rt4[,,"last"] = t(x[[hypothesis]]$estimate[,"look",]==nrow( x$look))
  #ar.inf.rt4 = ar.inf.rt4[order(x[[3]]$seed),,]
  n.group = nrow(x$par$group)
  mx.n.rg = x[[hypothesis]]$sample[,1:n.group] 
  if (x$type == "surv") times <- c(0,rowMeans(sapply(x[[hypothesis]]$trial,function(x) as.matrix(x$look[,c("t")]),simplify="array"),dims=2,na.rm=TRUE))
  
  ##
  ## boxplots of all sample sizes
  ##
  if(which.plot[1]){
    cex  = if(!any(names(mc)=="cex")){1}else{eval(mc$cex)}        
    ylim = if(!any(names(mc)=="ylim")){c(0,max(x$look$n))}else{eval(mc$ylim)}
    ylab = if(!any(names(mc)=="ylab")){"Sample size"}else{mc$ylab}
    layout(matrix(1:2,ncol=2),widths=c(2,n.group),heights=1)
    par(mar=c(2,4,3,0.25),omi=c(0,0,ifelse(title,0.25,0),.05))        
    #
    boxplot(apply(mx.n.rg,1,sum),ylab=ylab,xlab="",main="Total",
            ylim=ylim,axes=FALSE,col=col[3],cex=cex)
    box()       
    axis(2,las=2,cex.axis=.75)   
    axis(2,at=max(x$look$n),las=2)       
    #
    boxplot(mx.n.rg,ylab=ylab,xlab="",main="Per group",
            axes=FALSE,col=col[3])
    box()       
    axis(2,las=2)       
    axis(1,at=1:n.group,x$par$group$id)       
    #
    if(title){mtext(titlew,3,outer=TRUE,font=2,cex=1.5)}             
  }
  ##
  ## barplot of total sample size
  ##
  if(which.plot[2]){
    # 
    if (x$type=="glm") {
      total.r = apply(mx.n.rg,1,sum)        
      x$look$pi.hat = 0
      x$look$pi.hat[match(names(table(total.r)),x$look$n)] = table(total.r)/nrow(mx.n.rg)*100
      n.look = nrow(x$look)
      shift  = mean(x$look$n[-1]-x$look$n[-n.look])*0.2
      # options
      cex  = if(!any(names(mc)=="cex")){1}else{eval(mc$cex)}        
      ylim = if(!any(names(mc)=="ylim")){max(x$look$pi.hat)*c(0,1.05)}else{eval(mc$ylim)}
      xlim = if(!any(names(mc)=="xlim")){range(x$look$n)+c(-shift,shift)}else{eval(mc$xlim)}
      ylab = if(!any(names(mc)=="ylab")){"Frequency (%)"}else{mc$ylab}
      xlab = if(!any(names(mc)=="xlab")){"Cumulative total sample size per look (n)"}else{mc$xlab}        
      font = if(!any(names(mc)=="font")){2}else{eval(mc$font)}        
      lty  = if(!any(names(mc)=="lty")){2}else{eval(mc$lty)}                
      #
      layout(matrix(1:1,ncol=1))
      par(mar=c(4,4,1,0.25),omi=c(0,0,ifelse(title,0.25,0),.05))     
      plot(1, 1, pch = "", axes = FALSE, xlab = "", ylab = "",
           main = "", ylim = ylim, xlim = xlim)
      abline(h=0)
      #
      pos_ticks <- pretty(x$look$pi.hat, n=5)
      pos_ticks <- pos_ticks[pos_ticks >= 0 & pos_ticks< ylim[2]]
      axis(2,at=pos_ticks, las=2)        
      axis(2,at=ylim[2]/2,ylab,tick=FALSE,font=font,padj=-3)
      #
      axis(1,at=x$look$n, paste0(x$look$n,"\n(#",x$look$pos,")"), pos=0,padj=0.5)        
      axis(1,at=mean(xlim),xlab,tick=FALSE,font=font,padj=3.5,pos=0)        
      #
      for(lw in 1:n.look){
        rect(x$look$n[lw]-shift,0,
             x$look$n[lw]+shift,x$look$pi.hat[lw],col=col[3])
      }
      if(percentage){
        for(lw in 1:n.look){
          text(x$look$n[lw],x$look$pi.hat[lw],
               paste0(format(round(x$look$pi.hat[lw],2)),"%"),
               col=percentage.col,pos=3,cex=percentage.cex)
        }
      }
      #
      if(ess){
        segments(mean(total.r),0,
                 mean(total.r),ylim[2],col=ess.col,lty=lty)   
        axis(3,at=mean(total.r),paste0("ESS = ",round(mean(total.r),1)),
             col.axis=ess.col, tick=FALSE, pos=ylim[2], cex=ess.cex) 
      }
      # 
      if(title){mtext(titlew,3,outer=TRUE,font=2,cex=1.5)}    
    }
    if (x$type=="surv" && (is.null(x$par$interim$event.type) || x$par$interim$event.type!="cplusone")) {
      total.r = apply(mx.n.rg,1,sum)
      n.look = nrow(x$look)  
      looks = apply(x[[hypothesis]]$estimate[,"look",],2,max)
      x$look$pi.hat = 0
      x$look$pi.hat[match(names(table(looks)),as.character(1:n.look))] = table(looks)/nrow(x[[hypothesis]]$sample)*100
      
      scaling_factor <- 5
      
      # options
      cex  = if(!any(names(mc)=="cex")){1}else{eval(mc$cex)}        
      ylim = if(!any(names(mc)=="ylim")){c(0,100)*1.05}else{eval(mc$ylim)}
      xlim = if(!any(names(mc)=="xlim")){range(times)*1.1}else{eval(mc$xlim)}
      ylab = if(!any(names(mc)=="ylab")){"Stopping frequency (%)"}else{mc$ylab}
      xlab = if(!any(names(mc)=="xlab")){"Time"}else{mc$xlab}        
      font = if(!any(names(mc)=="font")){2}else{eval(mc$font)}        
      lty  = if(!any(names(mc)=="lty")){2}else{eval(mc$lty)}                
      #
      
      layout(matrix(1:1,ncol=1))
      par(mar=c(4,5,5.5,1),omi=c(0,0,ifelse(title,0.25,0),.05))  
      plot(1, 1, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, xlab = "", ylab = "", main = "")
      
      abline(v=round(rowMeans(sapply(x[[hypothesis]]$trial,function(x) as.matrix(x$look[,c("t")]),simplify="array"),dims=2,na.rm=TRUE),2),col="lightgrey",lty=2)
      
      width <- 1/(n.look+1)
      
      rectbase <- rep(0,n.look)
      wl <- times[-1]-width
      wr <- times[-1]+width
      rect(wl,rectbase,wr,x$look$pi.hat,col=col[3])
      
      # barplot(x$look$pi.hat, axes = FALSE, xlab = xlab, ylab = "",
      #         main = "", space=diff(times)/sum(diff(times))*scaling_factor, font.lab = 2,ylim=ylim)
      abline(h=0)
      #
      pos_ticks <- pretty(ylim, n=6)
      pos_ticks <- pos_ticks[pos_ticks >= 0 & pos_ticks<= ylim[2]]
      axis(2,at=c(0,20,40,60,80,100), las=2,pos=0)        
      axis(2,at=ylim[2]/2,ylab,tick=FALSE,font=font,padj=-3, pos=0)
      #
      x_ticks <- pretty(xlim)
      axis(1,at=x_ticks, pos=0,padj=0.5)        
      axis(1,at=mean(xlim),xlab,tick=FALSE,font=font,padj=3.5,pos=0)

      axis(3,at=times[-1],paste0("(#",x$look$pos,")\n\n\n"),pos=100,lwd=0)
      axis(3,at=times[-1],paste0("\n",
                                 round(rowMeans(sapply(x[[hypothesis]]$trial,function(x) as.matrix(x$look[,c("t")]),simplify="array"),dims=2,na.rm=TRUE),2),
                                 "\n",
                                 round(rowMeans(sapply(x[[hypothesis]]$trial,function(x) as.matrix(x$look[,c("n")]),simplify="array"),dims=2,na.rm=TRUE),0)),
           pos=100,line=0,cex.axis=.9,lwd=0)
      axis(3,at=0,paste0("\n\n\nt = \nn = "),lwd=0,line=0,pos=100,cex.axis=.9)
      mtext("Average timing of looks and sample sizes",side=3,font=2,line=4)        

      
      if(percentage){
        for(lw in 1:n.look){
          text(times[lw+1],x$look$pi.hat[lw],
               paste0(format(round(x$look$pi.hat[lw],2)),"%"),
               col=percentage.col,pos=3,cex=percentage.cex)
        }
      }
      #
      # if(ess){
      #   segments(mean(total.r),0,
      #            mean(total.r),ylim[2],col=ess.col,lty=lty)   
      #   axis(3,at=mean(total.r),paste0("ESS = ",round(mean(total.r),1)),
      #        col.axis=ess.col, tick=FALSE, pos=ylim[2], cex=ess.cex) 
      # }
      # # 
      if(title){mtext(titlew,3,outer=TRUE,font=2,cex=1.5)}    
    }
  }
  ##
  ## violin plots
  ##
  if(which.plot[3]){
    # options
    cex  = if(!any(names(mc)=="cex")){1}else{eval(mc$cex)}        
    ylim = if(!any(names(mc)=="ylim")){range(mx.n.rg)*c(0,1)}else{eval(mc$ylim)}
    xlim = if(!any(names(mc)=="xlim")){c(1,n.group)+c(-.5,.5)}else{eval(mc$xlim)}
    ylab = if(!any(names(mc)=="ylab")){"Sample size (n)"}else{mc$ylab}
    #xlab = if(!any(names(mc)=="xlab")){""}else{mc$xlab}        
    font = if(!any(names(mc)=="font")){2}else{eval(mc$font)}        
    lty  = if(!any(names(mc)=="lty")){2}else{eval(mc$lty)}                
    pch  = if(!any(names(mc)=="pch")){c(19,4)}else{eval(mc$pch)}        
    h    = if(!any(names(mc)=="h")){c(1)}else{eval(mc$h)}        
    pos_ticks <- pretty(c(ylim[1],as.matrix(mx.n.rg)), n=8)
    ylim = if(!any(names(mc)=="ylim")){range(pos_ticks)}
    lwd  = if(!any(names(mc)=="lwd")){1.25}
    
    if(legend){
      layout(matrix(1:2,ncol=1),heights=c(1,legend.height))
    }else{
      layout(matrix(1:1,ncol=1))
    }
    par(mar=c(1.5,4,0.25,0.25),omi=c(0,0,ifelse(title,0.25,0),.05))  
    plot(1, 1, pch = "", axes = FALSE, xlab = "", ylab = "",
         main = "", ylim = ylim, xlim = xlim)
    abline(h=pos_ticks,col=gray(.85),lwd=.85,lty=2)
    axis(2,las=2,at=pos_ticks)
    axis(2,at=mean(ylim),ylab,tick=FALSE,font=font,padj=-3)        
    axis(1,at=1:n.group,colnames(mx.n.rg),pos=0)
    
    for(gw in 1:n.group){   # gw=1        
      if(any(id.targetw$group==colnames(mx.n.rg)[gw])){
        colw  = col[c(3,1,2,4)][ar.inf.rt4[,which(id.targetw$group==colnames(mx.n.rg)[gw]),"type"]+1]
      }else{
        colw  = col[3]
      }
      violin    = violin.density(mx.n.rg[,gw], h=smooth)
      # points
      epsilon.x = rnorm(nrow(mx.n.rg),rep(0,nrow(mx.n.rg)),
                        sqrt(violin$pw.freq/max(violin$pw.freq))*.15)
      epsilon.y = runif(nrow(mx.n.rg),-.25,.25)
      points(rep(gw,nrow(mx.n.rg))+epsilon.x,mx.n.rg[,gw]+epsilon.y,
             col=colw,pch=pch)
      # violin
      xx = c(-violin$sm.density, violin$sm.density[length(violin$sm.density):1])+gw
      yy = c(violin$sm.points, violin$sm.points[length(violin$sm.points):1])
      polygon(xx, yy, col = NA, border = 1, lwd=1.75)
      # boxplot
      segments(gw,violin$quantile[1],gw,violin$quantile[2],lwd=1)
      segments(gw,violin$quantile[4],gw,violin$quantile[5],lwd=1)
      rect(gw-.025, violin$quantile[2], gw+.025, violin$quantile[4], 
           col=paste0(gray(.99),99), border=1)
      segments(gw-.025,violin$quantile[3],
               gw+.025,violin$quantile[3],col=1, lwd=2)
      # ess
      if(ess){
        points(gw,violin$mean,pch=15,col=ess.col)
        #segments(gw-.025,violin$mean,
        #         gw+.025,violin$mean,col=col[5], lwd=5)
        label        = paste0("ESS=",round(violin$mean,1))
        label_width  = strwidth(label,cex=ess.cex)
        label_height = strheight(label,cex=ess.cex)
        rect(gw+.05, violin$mean - label_height*3/4, 
             gw+.05+label_width*11/10, violin$mean + label_height*3/4, 
             col = ess.bg, border = NA)
        text(gw,violin$mean,pos=4,label,
             col=ess.col, cex=ess.cex)
      }
    }
    #
    if(title){mtext(titlew,3,outer=TRUE,font=2,cex=1.5)}     
    if(legend){            
      par(mar=c(0,4,0,0))
      plot(1, 1, pch = "", axes = FALSE, xlab = "", ylab = "",
           main = "", ylim = c(0,1), xlim = c(0,1))
      if(!any(ar.inf.rt4[,,"type"]==3)){
        if(all(nchar(col)[1:3]==9)){colw = substr(col[1:3],1,7)}else{colw=col}
        legend("top",legend=c("Efficacy","Futility","Neither"),
               col=colw,pch=15,box.lwd=NA,ncol=3,pt.cex=legend.pt,cex=legend.cex)
      }else{
        if(all(nchar(col)[1:4]==9)){colw = substr(col[1:4],1,7)}else{colw=col}
        legend("top",legend=c("Efficacy","Futility","Neither","Both"),
               col=colw,pch=15,box.lwd=NA,ncol=4,pt.cex=legend.pt,cex=legend.cex)                
      }
      legend("bottom",legend=c("Early stopping","Stopping at last look"),
             col=1,pch=pch,box.lwd=NA,ncol=2,pt.cex=legend.pt,cex=legend.cex)
    }
  }
  # estimates
  if(which.plot[4]){
    cex   = if(!any(names(mc)=="cex")){1}else{eval(mc$cex)}
    pch   = if(!any(names(mc)=="pch")){c(19,4)}else{eval(mc$pch)}        
    ylimw = if(!any(names(mc)=="ylim")){range(ar.inf.rt4[,,"est"])}else{eval(mc$ylim)}
    xlimw = range(ar.inf.rt4[,,"n"])
    if(legend){
      layout(rbind(matrix(1:n.targetw,ncol=n.targetw),n.targetw+1),
             heights=c(1,legend.height))
    }else{
      layout(matrix(1:n.targetw,ncol=n.targetw))
    }
    par(mar=c(4,4,3,1),omi=c(0,0,ifelse(title,0.25,0),.05))
    for(tw in 1:n.targetw){
      dataw = data.frame(n=ar.inf.rt4[,tw,"n"],
                         est=ar.inf.rt4[,tw,"est"],
                         type=ar.inf.rt4[,tw,"type"],
                         last=ar.inf.rt4[,tw,"last"],
                         col=col[c(3,1,2,4)][ar.inf.rt4[,tw,"type"]+1],
                         pch=pch[ar.inf.rt4[,tw,"last"]+1])
      plot(dataw$n,dataw$est,
           col=dataw$col,ylim=ylimw,xlim=xlimw,
           main=id.targetw$id[tw],xlab="Sample size",
           ylab = "Estimates",pch=dataw$pch,cex=cex)
      abline(h= if (x$type=="glm") id.targetw$beta[tw] else log(id.targetw$beta[tw]),col=beta.col,lty=beta.lty, lwd=beta.lwd)  
      #if(!is.null(beta.label){
      #    text(xlimw[1]+(xlimw[2]-xlimw[1]*.05),id.targetw$beta[tw],
      #         beta.label, col=beta.col, cex=beta.cex, pos=3)
      #}
      if(ess){   
        abline(v= mean(ar.inf.rt4[,tw,"n"]),col=ess.col,lty=3)
        axis(3,mean(ar.inf.rt4[,tw,"n"]),tick=FALSE,
             paste0("ESS = ",round(mean(ar.inf.rt4[,tw,"n"]),1)),
             col.axis=ess.col,padj=2,cex.axis=ess.cex)
      }
    }
    if(title){mtext(titlew,3,outer=TRUE,font=2,cex=1.5)}     
    if(legend){            
      par(mar=c(0,4,0,0))
      plot(1, 1, pch = "", axes = FALSE, xlab = "", ylab = "",
           main = "", ylim = c(0,1), xlim = c(0,1))
      if(!any(dataw$type==3)){
        if(all(nchar(col)[1:3]==9)){colw = substr(col[1:3],1,7)}else{colw=col}
        legend("top",legend=c("Efficacy","Futility","Neither"),
               col=colw,pch=15,box.lwd=NA,ncol=3,pt.cex=legend.pt,cex=legend.cex)
      }else{
        if(all(nchar(col)[1:4]==9)){colw = substr(col[1:4],1,7)}else{colw=col}
        legend("top",legend=c("Efficacy","Futility","Neither","Both"),
               col=colw,pch=15,box.lwd=NA,ncol=4,pt.cex=legend.pt,cex=legend.cex)                
      }
      legend("bottom",legend=c("Early stopping","Stopping at last look"),
             col=1,pch=pch,box.lwd=NA,ncol=2,pt.cex=legend.pt,cex=legend.cex)
    }
  }
  if(which.plot[5]){
    if (x$type!="surv") {
      stop("which = 5 currently only available for 'batss' object of type 'surv'")
    } else {
      sum.x <- summary(x)
      cum.eff <- sum.x[[hypothesis]]$cum.eff
      cum.fut <- sum.x[[hypothesis]]$cum.fut
      
      # options
      cex  = if(!any(names(mc)=="cex")){1}else{eval(mc$cex)}        
      ylim = if(!any(names(mc)=="ylim")){c(0,1)}else{eval(mc$ylim)}
      xlim = if(!any(names(mc)=="xlim")){c(0,max(times))*1.1}else{eval(mc$xlim)}
      ylab = if(!any(names(mc)=="ylab")){"Frequency (%)"}else{mc$ylab}
      xlab = if(!any(names(mc)=="xlab")){"Cumulative total sample size per look (n)"}else{mc$xlab}        
      font = if(!any(names(mc)=="font")){2}else{eval(mc$font)}        
      lty  = if(!any(names(mc)=="lty")){2}else{eval(mc$lty)}                
      #
      cols.eff <- rev(RColorBrewer::brewer.pal(n.group+3,"Greens")[-c(1:3,n.group+2)])
      cols.fut <- rev(RColorBrewer::brewer.pal(n.group+3,"OrRd")[-c(1:3,n.group+2)])
      
      if (legend) layout(matrix(c(1,2,3,3),2,2,byrow = TRUE), heights=c(1,legend.height)) else layout(matrix(c(1,2),ncol=2))
      par(mar=c(4,4,1,0.25),omi=c(0,0,ifelse(title,0.25,0),.05)) 
      # cumulative efficacy
      plot(times[-1],cum.eff[1,],axes=FALSE,type="n",xlab="Time",ylab="",xlim=xlim,ylim=ylim,main="Cumulative efficacy")
      axis(1,pos=0)
      axis(2,las=2,pos=0)
      abline(h=c(0.05,0.8),col="lightgrey",lty=2)
      for (i in 1:(n.group-1)) {
        lines(times[-1],cum.eff[i,],type="b",col=cols.eff[i],lwd=2,lty=i)
      }
      # cumulative futility
      plot(times[-1],cum.fut[1,],axes=FALSE,type="n",xlab="Time",xlim=xlim,ylab="",ylim=ylim,main="Cumulative futility")
      axis(1,pos=0)
      axis(2,las=2,pos=0)
      abline(h=0.5,col="lightgrey",lty=2)
      for (i in 1:(n.group-1)) {
        lines(times[-1],cum.fut[i,],type="b",col=cols.fut[i],lwd=2,lty=i)
      }
      if(title){mtext(titlew,3,outer=TRUE,font=2,cex=1.5)}  
      if (legend) {
        par(mar=c(0,4,0,0))
        plot(1, 1, pch = "", axes = FALSE, xlab = "", ylab = "",
             main = "", ylim = c(0,1), xlim = c(0,1))
        cols.gray <- gray.colors(n.group-1,end=0.8)
        legend("bottom",legend=x$par$group$id[-1],
               col=cols.gray,box.lwd=NA,ncol=n.group-1,pt.cex=legend.pt,cex=legend.cex,lwd=3,lty=1:(n.group-1))
      }
    }
  }
}


#' @name violin.density
#' @title Violin plot density estimates
#' @description Extracts density estimates for a violin plot. 
#' @param data data vector.
#' @param h a scalar corresponding to the smoothing parameter.
#' @returns a list of coordinates and corresponding density values.
#' @seealso \link[sm]{sm.density} 
violin.density = function(data, h=NULL){
  data      = data[!is.na(data)]
  quantilew = quantile(data, prob=c(0,.25,.5,.75,1))
  iqrw      = quantilew[4] - quantilew[2]
  upper     = min(quantilew[4] + 1.5 * iqrw, quantilew[5])
  lower     = max(quantilew[2] - 1.5 * iqrw, quantilew[1])
  rangew    = c(min(lower, quantilew[1]), max(upper, quantilew[5]))
  argw      = list(display="none",h=h)
  # sm
  if(is.null(h)){
    smoothw   = do.call(sm::sm.density, c(list(data, xlim = rangew, display="none")))
  }else{
    smoothw   = do.call(sm::sm.density, c(list(data, xlim = rangew, display="none", h=h)))
  }
  # block
  n.lim   = min(100,length(unique(data)))
  lim     = seq(quantilew[1],quantilew[5],length=n.lim)
  pw.data = rep(n.lim, length(data))
  for(lw in n.lim:1){
    pw.data[data<=lim[lw]] = lw
  }
  pw.freq = tabulate(pw.data)/length(data)
  # 
  list(sm.points  = smoothw$eval.points,
       sm.density = smoothw$estimate * 0.4/max(smoothw$estimate),
       quantile   = c(lower,quantilew[2:4],upper), range=rangew, mean=mean(data),
       data       = data, pw.freq = pw.freq[pw.data])
}
