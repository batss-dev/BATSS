#' @name plot.batss
#' @title Plot function for BATSS outputs
#' @description Plot for objects of class 'batss' 
#' @param x An object of class 'batss' (i.e., output of the function [batss.glm]).
#' @param focus A character string indicating the type of plot with options 'size' (default) to display the total and per group sample size observed in the Monte Carlo trials, and 'estimates' to display the Monte Carlo trial target estimates as a function of the sample size.   
#' @param hypothesis A character string indicating which alternative hypothesis to use for analyses considering both "H0" and "H1", with options "H1" (default) and "H0".
#' @param ... Additional arguments affecting the plot produced.
#' @seealso [batss.glm()], the function generating S3 objects of class 'batss'. 
#' @export
plot.batss = function(x,focus="size",hypothesis="H1", ...){
    # hypothesis 
    pos = match(names(x),c("H0","H1"))
    if(sum(!is.na(pos))!=2){
        hypothesis = names(x)[!is.na(pos)]
        }
    # col
    colw = c("#EC008C","#00B6ED","#2E008B")
    # sample size
    if(focus=="size"){
        mx.n.rg = x[[hypothesis]]$sample[,-(ncol(x[[hypothesis]]$sample)+(-1:0))]
        n.g  = ncol(mx.n.rg)
        ylim = range(x$look$n)
        layout(matrix(1:2,ncol=2),widths=c(2,n.g),heights=1)
        par(mar=c(2,4,3,0),omi=c(0,0,0.25,.05))        
        #
        boxplot(apply(mx.n.rg,1,sum),ylab="Sample size",xlab="",main="Total",
               ylim=ylim,axes=FALSE,col=paste0(colw[3],50))
        box()       
        axis(2,las=2,cex.axis=.75)   
        axis(2,at=max(x$look$n),las=2)       
        abline(h=ylim[2],lty=2,col=colw[1])
        #
        boxplot(mx.n.rg,ylab="Sample size",xlab="Groups",main="Per group",
                axes=FALSE,col=paste0(colw[3],50))
        box()       
        axis(2,las=2)       
        axis(1,at=1:n.g,colnames(mx.n.rg))       
        abline(h=ylim[2],lty=2,col=colw[1])
        #
        mtext(paste0("Under '",hypothesis,"'"),3,outer=TRUE)        
        }
    # estimates
    if(focus=="estimates"){
        id.targetw = x[[hypothesis]]$target$par
        n.targetw  = nrow(id.targetw)
        id.targetw$beta = x$beta[id.targetw$id,grepl(hypothesis,colnames(x$beta))]
        ar.inf.rt3 = array(dim=c(dim(x[[hypothesis]]$sample)[1],n.targetw,3),
                           dimnames=list(1:dim(x[[hypothesis]]$sample)[1],
                                         x[[hypothesis]]$target$id,
                                         c("n","est","type")))
        ar.inf.rt3[,,"n"]    = as.matrix(x[[hypothesis]]$sample[,x[[hypothesis]]$target$par$group]) 
        ar.inf.rt3[,,"est"]  = t(x[[hypothesis]]$estimate[,"mid",])
        ar.inf.rt3[,,"type"] = t(x[[hypothesis]]$estimate[,"type",])
        ylimw = range(ar.inf.rt3[,,"est"])
        xlimw = range(ar.inf.rt3[,,"n"])
        par(mfrow=c(1,n.targetw),mar=c(4,4,3,1),omi=c(0,0,0.25,.05))
        for(tw in 1:n.targetw){
            dataw = data.frame(n=ar.inf.rt3[,tw,"n"],
                               est=ar.inf.rt3[,tw,"est"],
                               type=ar.inf.rt3[,tw,"type"],
                               col=c("#8B8970","#008B00","#8B3A3A")[ar.inf.rt3[,tw,"type"]+1])
            plot(dataw$n,dataw$est,
                 col=paste0(dataw$col,50),pch=19,ylim=ylimw,xlim=xlimw,
                 main=id.targetw$id[tw],xlab="Sample size",
                 ylab = "Estimates")
            abline(h= id.targetw$beta[tw],col="blue",lty=1)     
            abline(v= mean(ar.inf.rt3[,tw,"n"]),col="blue",lty=3)
            axis(3,mean(ar.inf.rt3[,tw,"n"]),tick=FALSE,
                 paste0("ESS = ",round(mean(ar.inf.rt3[,tw,"n"]),1)),
                 col.axis="blue",padj=2,cex.axis=.75)
            legend("topright",title="Stop:",legend=c("Efficacy","None","Futility"),
                   col=c("#008B00","#8B8970","#8B3A3A"),pch=16,box.lwd=NA)
            #spl <- with(dataw[], smooth.spline(n, est,penalty=2))
            #lines(spl,col=colw[2],lwd=1.5)
        }
        mtext(paste0("Under '",hypothesis,"'"),3,outer=TRUE)
    }
}
