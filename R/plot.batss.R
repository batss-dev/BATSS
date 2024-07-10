#' @name plot.batss
#' @title Plot function for 'BATSS' outputs
#' @description Plot for objects of class 'batss' 
#' @param x An object of class 'batss' (i.e., output of the function [batss.glm]).
#' @param type A character string indicating the type of plot with options 'size' (default) to display the total and per group sample size observed in the Monte Carlo trials, and 'estimates' to display the Monte Carlo trial target estimates as a function of the sample size.   
#' @param hypothesis A character string indicating which alternative hypothesis to use for analyses considering both "H0" and "H1", with options "H1" (default) and "H0".
#' @param title Either a \link[base]{logical} indicating if a title should be added or a string (of class \link[base]{character}) indicating the title to be added. If `title` equals `TRUE` (default), the title 'Under 'H1' or 'Under 'H0' (depending on the argument `hypothesis`) is added to the outer margin of the plot. No outer margin space is added if `title = FALSE`.
#' @param legend a \link[base]{logical} (with default set to `TRUE`) indicating if, when '`type = estimates`', a legend should be added at the bottom of the plot. 
#' @param col a vector of length 5 specifiying the colour respectively assigned to i/ efficacy, ii/ futility, iii/ neither or iv/ both when color-coding trials when `type = "estimates"`. The 5th colour is used for lines. Default to c("#8B897040","#008B0040","#8B3A3A40","#FF990075","blue") where the 2 last digits of the long hexadecimal strings of colours 1 to 4 specify the level of transluency. Refer to the Section 'colour specification' in \link[graphics]{par} for details. If the length of `col` equals 1, the same colour is used for all cases. When `type = "size"`, the 3rd and 5th colours of the vector `col` are used to display boxplots and lines.
#' @param ... Additional arguments affecting the plot produced, like ylim and ylab.
#' @returns Generates graphical displays of results for objects of class 'batss'.
#' @seealso [batss.glm()], the function generating S3 objects of class 'batss'. 
#' @export
plot.batss = function(x, type="size", hypothesis="H1", title=TRUE, legend=TRUE, 
                      col = c("#008B0040","#8B3A3A40","#8B897040","#FF990075","blue"),
                      ...){
    oldpar <- par(no.readonly = TRUE) 
    on.exit(par(oldpar)) 
    mc = match.call()    
    # hypothesis 
    pos = match(names(x),c("H0","H1"))
    if(sum(!is.na(pos))!=2){
        hypothesis = names(x)[!is.na(pos)]
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
    if(length(col)==1){col = rep(col,5)
    }else{
        if(length(col)<5){
            col[(length(col)+1):5] = c("#008B0040","#8B3A3A40","#8B897040","#FF990075","blue")[(length(col)+1):5]
        }
    }
    # sample size
    if(type=="size"){
        mx.n.rg = x[[hypothesis]]$sample[,-(ncol(x[[hypothesis]]$sample)+(-1:0))]
        n.g  = ncol(mx.n.rg)
        cex  = if(!any(names(mc)=="cex")){1}else{eval(mc$cex)}        
        ylim = if(!any(names(mc)=="ylim")){range(x$look$n)}else{eval(mc$ylim)}
        ylab = if(!any(names(mc)=="ylab")){"Sample size"}else{mc$ylab}
        layout(matrix(1:2,ncol=2),widths=c(2,n.g),heights=1)
        par(mar=c(2,4,3,0.25),omi=c(0,0,ifelse(title,0.25,0),.05))        
        #
        boxplot(apply(mx.n.rg,1,sum),ylab=ylab,xlab="",main="Total",
                ylim=ylim,axes=FALSE,col=col[3],cex=cex)
        box()       
        axis(2,las=2,cex.axis=.75)   
        axis(2,at=max(x$look$n),las=2)       
        abline(h=ylim[2],lty=2,col=col[5])
        #
        boxplot(mx.n.rg,ylab=ylab,xlab="",main="Per group",
                axes=FALSE,col=col[3])
        box()       
        axis(2,las=2)       
        axis(1,at=1:n.g,colnames(mx.n.rg))       
        abline(h=ylim[2],lty=2,col=col[5])
        #
        if(title){mtext(titlew,3,outer=TRUE,font=2,cex=1.5)}     
        }
    # estimates
    if(type=="estimates"){
        id.targetw = x[[hypothesis]]$target$par
        n.targetw  = nrow(id.targetw)
        id.targetw$beta = x$beta[id.targetw$id,grepl(hypothesis,colnames(x$beta))]
        ar.inf.rt4 = array(dim=c(dim(x[[hypothesis]]$sample)[1],n.targetw,4),
                           dimnames=list(1:dim(x[[hypothesis]]$sample)[1],
                                         x[[hypothesis]]$target$id,
                                         c("n","est","type","last")))
        ar.inf.rt4[,,"n"]    = as.matrix(x[[hypothesis]]$sample[,x[[hypothesis]]$target$par$group]) 
        ar.inf.rt4[,,"est"]  = t(x[[hypothesis]]$estimate[,"mid",])
        ar.inf.rt4[,,"type"] = t(x[[hypothesis]]$estimate[,"type",])
        ar.inf.rt4[,,"last"] = t(x[[hypothesis]]$estimate[,"look",]==nrow( x$look))
        cex   = if(!any(names(mc)=="cex")){1}else{eval(mc$cex)}
        pch   = if(!any(names(mc)=="pch")){c(19,4)}else{eval(mc$pch)}        
        ylimw = if(!any(names(mc)=="ylim")){range(ar.inf.rt4[,,"est"])}else{eval(mc$ylim)}
        xlimw = range(ar.inf.rt4[,,"n"])
        if(legend){
            layout(rbind(matrix(1:n.targetw,ncol=n.targetw),n.targetw+1),
                   heights=c(1,.1))
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
            abline(h= id.targetw$beta[tw],col=col[5],lty=1)     
            abline(v= mean(ar.inf.rt4[,tw,"n"]),col=col[5],lty=3)
            axis(3,mean(ar.inf.rt4[,tw,"n"]),tick=FALSE,
                 paste0("ESS = ",round(mean(ar.inf.rt4[,tw,"n"]),1)),
                 col.axis=col[5],padj=2,cex.axis=.75)
        }
        if(title){mtext(titlew,3,outer=TRUE,font=2,cex=1.5)}     
        if(legend){            
            par(mar=c(0,4,0,0))
            plot(1, 1, pch = "", axes = FALSE, xlab = "", ylab = "",
                 main = "", ylim = c(0,1), xlim = c(0,1))
            if(!any(dataw$type==3)){
                if(all(nchar(col)[1:3]==9)){colw = substr(col[1:3],1,7)}else{colw=col}
                legend("top",legend=c("Efficacy","Futility","Neither"),
                       col=colw,pch=15,box.lwd=NA,ncol=3,pt.cex=2.5,cex=1.25)
                }else{
                if(all(nchar(col)[1:4]==9)){colw = substr(col[1:4],1,7)}else{colw=col}
                legend("top",legend=c("Efficacy","Futility","Neither","Both"),
                       col=colw,pch=15,box.lwd=NA,ncol=4,pt.cex=2.5,cex=1.25)                
            }
            legend("bottom",legend=c("Early stopping","Stopping at last look"),
                   col=1,pch=pch,box.lwd=NA,ncol=2,pt.cex=2,cex=1.25)
        }
    }
}
