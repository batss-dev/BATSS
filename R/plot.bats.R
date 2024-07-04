#' @name plot.batss
#' @title Plot function for BATSS outputs
#' @description Plot for objects of class 'batss' 
#' @param x An object of class 'batss' (i.e., output of the function [batss.glm]).
#' @param type A character string indicating the type of plot with options 'size' (default) to display the total and per group sample size observed in the Monte Carlo trials, and 'estimates' to display the Monte Carlo trial target estimates as a function of the sample size.   
#' @param hypothesis A character string indicating which alternative hypothesis to use for analyses considering both "H0" and "H1", with options "H1" (default) and "H0".
#' @param title Either a \link[base]{logical} indicating if a title should be added or a string (of class \link[base]{character}) indicating the title to be added. If `title` equals `TRUE` (default), the title 'Under 'H1' or 'Under 'H0' (depending on the argument `hypothesis`) is added to the outer margin of the plot. No outer margin space is added if `title = FALSE`.
#' @param legend a \link[base]{logical} (with default set to `TRUE`) indicating if, when '`type = estimates`', a legend should be added at the bottom of the plot. 
#' @param ... Additional arguments affecting the plot produced, like ylim and ylab.
#' @seealso [batss.glm()], the function generating S3 objects of class 'batss'. 
#' @export
plot.batss = function(x,type="size",hypothesis="H1", title=TRUE, legend=TRUE, ...){
    # hypothesis 
    pos = match(names(x),c("H0","H1"))
    if(sum(!is.na(pos))!=2){
        hypothesis = names(x)[!is.na(pos)]
        }
    # col
    colw  = c("#EC008C","#00B6ED","#2E008B")
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
    mc = match.call() 
    # sample size
    if(type=="size"){
        mx.n.rg = x[[hypothesis]]$sample[,-(ncol(x[[hypothesis]]$sample)+(-1:0))]
        n.g  = ncol(mx.n.rg)
        ylim = if(!any(names(mc)=="ylim")){range(x$look$n)}else{eval(mc$ylim)}
        ylab = if(!any(names(mc)=="ylab")){"Sample size"}else{mc$ylab}
        layout(matrix(1:2,ncol=2),widths=c(2,n.g),heights=1)
        par(mar=c(2,4,3,0.25),omi=c(0,0,ifelse(title,0.25,0),.05))        
        #
        boxplot(apply(mx.n.rg,1,sum),ylab=ylab,xlab="",main="Total",
                ylim=ylim,axes=FALSE,col=paste0(colw[3],50))
        box()       
        axis(2,las=2,cex.axis=.75)   
        axis(2,at=max(x$look$n),las=2)       
        abline(h=ylim[2],lty=2,col=colw[1])
        #
        boxplot(mx.n.rg,ylab=ylab,xlab="",main="Per group",
                axes=FALSE,col=paste0(colw[3],50))
        box()       
        axis(2,las=2)       
        axis(1,at=1:n.g,colnames(mx.n.rg))       
        abline(h=ylim[2],lty=2,col=colw[1])
        #
        if(title){mtext(titlew,3,outer=TRUE,font=2,cex=1.5)}     
        }
    # estimates
    if(type=="estimates"){
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
        ylimw = if(!any(names(mc)=="ylim")){range(ar.inf.rt3[,,"est"])}else{eval(mc$ylim)}        
        xlimw = range(ar.inf.rt3[,,"n"])
        if(legend){
            layout(rbind(matrix(1:n.targetw,ncol=n.targetw),n.targetw+1),
                   heights=c(1,.1))
        }else{
            layout(matrix(1:n.targetw,ncol=n.targetw))
        }
        par(mar=c(4,4,3,1),omi=c(0,0,ifelse(title,0.25,0),.05))
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
        }
        if(title){mtext(titlew,3,outer=TRUE,font=2,cex=1.5)}     
        if(legend){
            par(mar=c(0,4,0,0))
            plot(1, 1, pch = "", axes = FALSE, xlab = "", ylab = "",
                 main = "", ylim = c(0,1), xlim = c(0,1))
            legend("bottom",title="Reason for early stopping:",legend=c("Efficacy","None","Futility"),
                   col=c("#008B00","#8B8970","#8B3A3A"),pch=15,box.lwd=NA,ncol=3,pt.cex=2.5,cex=1.25)
        }
    }
}
