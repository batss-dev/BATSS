#' @name summary.batss
#' @title Summary function for 'BATSS' outputs
#' @description Summary method function for objects of class 'batss'.
#' @param object An object of class 'batss' (i.e., output of the function [batss.glm]).
#' @param ... For future use
#' @returns Object of class 'summary.batss'.
#' @seealso [batss.glm()], the function generating S3 objects of class 'batss'. 
#' @export
summary.batss = function(object, ...){
  res <- list()
  res$call <- object$call
  res$par <- object$par
  res$type <- object$type
  if (object$type=="surv") {
    res$hr <- object$hr
  } else {
    res$beta <- object$beta
    res$sample <- object$look
    colnames(res$sample)[1] = ""
  }
  if (object$par$H0) {
    
    res$H0$sample.sizes <- object$H0$sample[,1:ifelse(object$type=="surv",(sum(object$hr$target)+1),sum(object$beta$target))]
    
    if (object$type=="surv") {
      res$H0$events <- object$H0$sample[,(2*(sum(object$hr$target)+1)+1):(3*(sum(object$hr$target)+1))]
      res$H0$t.obs <- object$H0$sample[,(sum(object$hr$target)+2):(2*(sum(object$hr$target)+1))]
      res$H0$t.trial <- object$H0$sample[,"t"]
    }
    
    temp = rbind(object$H0$efficacy$par,object$H0$efficacy$global)
    temp[,1][-(1:nrow(object$H0$efficacy$par))] = ""
    
    res$H0$efficacy <- temp
    
    temp = rbind(object$H0$futility$par,object$H0$futility$global)
    temp[,1][-(1:nrow(object$H0$futility$par))] = ""
    
    res$H0$futility <- temp
    
    temp_col = (ifelse(object$type=="surv",nrow(object$hr)+1,nrow(object$beta)))
    temp = object$H0$sample[,1:temp_col]
    temp = cbind(temp,Total=apply(temp,1,sum))
    temp = data.frame(pos=1:ncol(temp),
                      id=colnames(temp),
                      'ESS'=round(apply(temp,2,mean),2),
                      'St.Dev'=round(sqrt(apply(temp,2,var)),2),
                      'q10'=round(apply(temp,2,quantile,probs=0.1),2),
                      'q50'=round(apply(temp,2,quantile,probs=0.5),2),
                      'q90'=round(apply(temp,2,quantile,probs=0.9),2))
    colnames(temp)[c(1,5:7)] = c("",paste0("q(",c(0.1,0.5,0.9),")"))
    temp[,1][nrow(temp)] = ""
    
    res$H0$summary.sample.sizes <- temp
    
    if (object$type=="surv") {
      temp = object$H0$sample[,(2*temp_col+1):(3*temp_col)]
      temp = cbind(temp,Total=apply(temp,1,sum))
      temp = data.frame(pos=1:ncol(temp),
                        id=colnames(temp),
                        'Mean'=round(apply(temp,2,mean),2),
                        'St.Dev'=round(sqrt(apply(temp,2,var)),2),
                        'q10'=round(apply(temp,2,quantile,probs=0.1),2),
                        'q50'=round(apply(temp,2,quantile,probs=0.5),2),
                        'q90'=round(apply(temp,2,quantile,probs=0.9),2))
      colnames(temp)[c(1,5:7)] = c("",paste0("q(",c(0.1,0.5,0.9),")"))
      temp[,1][nrow(temp)] = ""
      
      res$H0$summary.events <- temp
      
      temp = object$H0$sample[,(temp_col+1):(2*temp_col)]
      temp = cbind(temp,Total=apply(temp,1,sum))
      temp = data.frame(pos=1:ncol(temp),
                        id=colnames(temp),
                        'Mean'=round(apply(temp,2,mean),2),
                        'St.Dev'=round(sqrt(apply(temp,2,var)),2),
                        'q10'=round(apply(temp,2,quantile,probs=0.1),2),
                        'q50'=round(apply(temp,2,quantile,probs=0.5),2),
                        'q90'=round(apply(temp,2,quantile,probs=0.9),2))
      colnames(temp)[c(1,5:7)] = c("",paste0("q(",c(0.1,0.5,0.9),")"))
      temp[,1][nrow(temp)] = ""
      
      res$H0$summary.t.obs <- temp
    }
    
    res$H0$scenario <- object$H0$scenario
    
    if (object$type=="surv" && !is.null(object$H0$trial)) {
      if (is.null(object$par$interim$event) || (!is.null(object$par$interim$event.type) && object$par$interim$event.type!="cplusone")) {
        objectw = rowMeans(sapply(object$H0$trial,function(x) as.matrix(cbind(x$look[,c("n","m","t","t(n)","ev(n)")],!is.na(x$look[,c("t")]))),
                                        simplify="array"),dims=2,na.rm=TRUE)
      
      objectw = cbind(objectw[,-6],abs(diff(c(objectw[,6]*length(object$par$seed),0))),objectw[,6]*length(object$par$seed),objectw[,6])
      colnames(objectw)[6:8] <- c("stops","n.trial","prcnt.trial")
      } else {
        objectw <- NULL
      }
      
      res$H0$avg.sample <- objectw
    }
    
    if (object$type=="surv") {
      res$H0$est.error <- t(sweep(exp(object$H0$estimate[,3,]),1,object$hr[,4]))
    } else {
      res$H0$est.error <- t(sweep(object$H0$estimate[,3,],1,object$beta[,4]))
    }
    
  } else {
    res$H0 <- NULL 
  } 
  
  if(object$par$H1){
    res$H1$sample.sizes <- object$H1$sample[,1:ifelse(object$type=="surv",(sum(object$hr$target)+1),sum(object$beta$target))]
    
    if (object$type=="surv") {
      res$H1$events <- object$H1$sample[,(2*(sum(object$hr$target)+1)+1):(3*(sum(object$hr$target)+1))]
      res$H1$t.obs <- object$H1$sample[,(sum(object$hr$target)+2):(2*(sum(object$hr$target)+1))]
      res$H1$t.trial <- object$H1$sample[,"t"]
    }
    
    temp = rbind(object$H1$efficacy$par,object$H1$efficacy$global)
    temp[,1][-(1:nrow(object$H1$efficacy$par))] = ""
    
    res$H1$efficacy <- temp
    
    temp = rbind(object$H1$futility$par,object$H1$futility$global)
    temp[,1][-(1:nrow(object$H1$futility$par))] = ""
    
    res$H1$futility <- temp
    
    temp_col = (ifelse(object$type=="surv",nrow(object$hr)+1,nrow(object$beta)))
    temp = object$H1$sample[,1:temp_col]
    temp = cbind(temp,Total=apply(temp,1,sum))
    temp = data.frame(pos=1:ncol(temp),
                      id=colnames(temp),
                      'ESS'=round(apply(temp,2,mean),2),
                      'St.Dev'=round(sqrt(apply(temp,2,var)),2),
                      'q10'=round(apply(temp,2,quantile,probs=0.1),2),
                      'q50'=round(apply(temp,2,quantile,probs=0.5),2),
                      'q90'=round(apply(temp,2,quantile,probs=0.9),2))
    colnames(temp)[c(1,5:7)] = c("",paste0("q(",c(0.1,0.5,0.9),")"))
    temp[,1][nrow(temp)] = ""
    
    res$H1$summary.sample.sizes <- temp
    
    if (object$type=="surv") {
      temp = object$H1$sample[,(2*temp_col+1):(3*temp_col)]
      temp = cbind(temp,Total=apply(temp,1,sum))
      temp = data.frame(pos=1:ncol(temp),
                        id=colnames(temp),
                        'Mean'=round(apply(temp,2,mean),2),
                        'St.Dev'=round(sqrt(apply(temp,2,var)),2),
                        'q10'=round(apply(temp,2,quantile,probs=0.1),2),
                        'q50'=round(apply(temp,2,quantile,probs=0.5),2),
                        'q90'=round(apply(temp,2,quantile,probs=0.9),2))
      colnames(temp)[c(1,5:7)] = c("",paste0("q(",c(0.1,0.5,0.9),")"))
      temp[,1][nrow(temp)] = ""
      
      res$H1$summary.events <- temp
      
      temp = object$H1$sample[,(temp_col+1):(2*temp_col)]
      temp = cbind(temp,Total=apply(temp,1,sum))
      temp = data.frame(pos=1:ncol(temp),
                        id=colnames(temp),
                        'Mean'=round(apply(temp,2,mean),2),
                        'St.Dev'=round(sqrt(apply(temp,2,var)),2),
                        'q10'=round(apply(temp,2,quantile,probs=0.1),2),
                        'q50'=round(apply(temp,2,quantile,probs=0.5),2),
                        'q90'=round(apply(temp,2,quantile,probs=0.9),2))
      colnames(temp)[c(1,5:7)] = c("",paste0("q(",c(0.1,0.5,0.9),")"))
      temp[,1][nrow(temp)] = ""
      
      res$H1$summary.t.obs <- temp
    }
    
    res$H1$scenario <- object$H1$scenario
    
    if (object$type=="surv" && !is.null(object$H1$trial)) {
      if (is.null(object$par$interim$event) || (!is.null(object$par$interim$event.type) && object$par$interim$event.type!="cplusone"))  {
        objectw = rowMeans(sapply(object$H1$trial,function(x) as.matrix(cbind(x$look[,c("n","m","t","t(n)","ev(n)")],!is.na(x$look[,c("t")]))),simplify="array"),dims=2,na.rm=TRUE)
        objectw = cbind(objectw[,-6],abs(diff(c(objectw[,6]*length(object$par$seed),0))),objectw[,6]*length(object$par$seed),objectw[,6])
        colnames(objectw)[6:8] <- c("stops","n.trial","prcnt.trial")
      } else {
        objectw <- NULL
      }
      res$H1$avg.sample <- objectw
    }
    
    if (object$type=="surv") {
      res$H1$est.error <- if (object$par$H0 && object$par$H1) t(sweep(exp(object$H1$estimate[,3,]),1,object$hr[,5])) else t(sweep(exp(object$H1$estimate[,3,]),1,object$hr[,4]))
    } else {
      res$H1$est.error <- if (object$par$H0 && object$par$H1) t(sweep(object$H1$estimate[,3,],1,object$beta[,5])) else t(sweep(object$H1$estimate[,3,],1,object$beta[,4]))
    }
    
  } else {
    res$H1 <- NULL 
  } 
  class(res) = "summary.batss"
  return(res)
}
  
#' @name print.summary.batss
#' @title Print the summary function for 'BATSS' outputs
#' @description Print method function for objects of class 'summary.batss'.
#' @param x An object of class 'summary.batss' (i.e., output of the function [batss.glm]).
#' @param full set to 'TRUE' the function prints a long version of the summary
#' @param ... For future use
#' @seealso [batss.glm()], the function generating S3 objects of class 'batss'. 
#' @export 
print.summary.batss = function(x, full=FALSE, ...){
  # common part
  cat("\n")
  if(!is.null(x$par$RAR)){
    cli_h1("Bayesian Adaptive Design with Laplace Approx.")
  }else{
    cli_h1("MAMS with Laplace Approx.")
  }
  cat("  (",length(x$par$seed)," Monte Carlo samples)\n",sep="")
  cat("\n")
  cli_h3("Variables:")
  for(i in 2:length(x$call$var)){
    cat("  *",names(x$call$var)[i],":",as.character(x$call$var)[i],"\n")
  }
  if(!is.null(x$par$RAR)){
    cat("\n")
    cli_h3("Group randomisation:")
    cat("  *",x$call$RAR,"\n")
  }
  if(!(is.null(x$par$eff.arm) && is.null(x$par$fut.arm))){
    cat("\n")
    cli_h3("Decision rules:")
    if(!is.null(x$par$eff.arm)) cat("  * Efficacy: ",format(x$call$eff.arm),"\n")
    if(!is.null(x$par$fut.arm)) cat("  * Futility: ",format(x$call$fut.arm),"\n")
  }
  cat("\n")
  cli_h3("Model: ")
  cat("  *",format(x$call$model),if (x$type=="surv") "\n" else paste("(with",ifelse(is.null(x$call$link),"identity",x$call$link), "link)\n"))
  cat("\n")
  cli_h3("Fixed effect parameters:\n")
  if (x$type=="surv") xw = x$hr else xw = x$beta
  colnames(xw)[1] = ""
  print(xw,row.names=FALSE)
  if (x$type!="surv"){
    cat("\n")
    cli_h3("Sample size per interim analyis:\n")
    print(x$sample,row.names=FALSE)
  }
  # H0
  if(x$par$H0){
    cat("\n\n")
    #cat(paste0(rep("-",floor(options()$width/2)),collapse=""))
    cli_h2("\n H0: Under the null hypothesis\n")
    #cat(paste0(rep("-",floor(options()$width/2)),collapse=""))
    #
    if (x$type=="surv") {
      cat("\n")
      cli_h3("Average sample size per interim analyis:\n")
      print(x$H0$avg.sample,row.names=FALSE,digits=4)
    }
    cat("\n")
    cli_h3("Target parameters:\n")
    temp = cbind(x$H0$efficacy[,-(5:6)],x$H0$futility[,7])
    colnames(temp)[c(1,5,6)] = c("","efficacy","futility")
    print(temp,row.names=FALSE)
    if (x$type=="surv") {
      cat("\n")
      cli_h3(paste0( "Trial duration:\n"))
      print(summary(x$H0$t.trial))
    }
    #
    if(full){
      cat("\n")
      cli_h3("Efficacy:\n")
      print(x$H0$efficacy,row.names=FALSE)
      #
      cat("\n")
      cli_h3("Futility:\n")
      print(x$H0$futility,row.names=FALSE)
      #
      cat("\n")
      cli_h3("Sample size per group:\n")
      print(x$H0$summary.sample.sizes,row.names=FALSE)
      #cat("\n")
      if (x$type=="surv") {
        cli_h3("Events per group:\n")
        print(x$H0$summary.events,row.names=FALSE)
        #cat("\n")
        cli_h3("Time observed per group:\n")
        print(x$H0$summary.t.obs,row.names=FALSE)
      }
      #
      cat("\n")
      cli_h3("Scenarios:\n")
      print(x$H0$scenario,row.names=FALSE)
      cat(" where 0 = no stop, 1 = efficacy stop, 2 = futility stop\n")
      if(any(x$H0$scenario[,x$H0$target$par$id]==3)){
        cat(",\n       3 = simultaneous efficacy and futility stops")
      }else{cat("\n")}
    }
  }
  # H1
  if(x$par$H1){
    cat("\n\n")
    #cat(paste0(rep("-",floor(options()$width/2)),collapse=""))
    cli_h2("\n H1: Under the alternative hypothesis\n")#cat("\n H1: Under the alternative hypothesis\n")
    #cat(paste0(rep("-",floor(options()$width/2)),collapse=""))
    #
    if (x$type=="surv") {
      cat("\n")
      cli_h3("Average sample size per interim analyis:\n")
      print(x$H1$avg.sample,row.names=FALSE,digits=4)
    }
    cat("\n")
    cli_h3("Target parameters:\n")
    temp = cbind(x$H1$efficacy[,-(5:6)],x$H1$futility[,7])
    colnames(temp)[c(1,5,6)] = c("","efficacy","futility")
    print(temp,row.names=FALSE)
    if (x$type=="surv") {
      cat("\n")
      cli_h3(paste0( "Trial duration:\n"))
      print(summary(x$H1$t.trial))
    }
    #
    if(full){
      cat("\n")
      cli_h3("Efficacy:\n")
      print(x$H1$efficacy,row.names=FALSE)
      #
      cat("\n")
      cli_h3("Futility:\n")
      print(x$H1$futility,row.names=FALSE)
      #
      cat("\n")
      cli_h3("Sample size per group:\n")
      print(x$H1$summary.sample.sizes,row.names=FALSE)
      #cat("\n")
      if (x$type=="surv") {
        cli_h3("Events per group:\n")
        print(x$H1$summary.events,row.names=FALSE)
        #cat("\n")
        cli_h3("Time observed per group:\n")
        print(x$H1$summary.t.obs,row.names=FALSE)
      }
      #
      cat("\n")
      cli_h3("Scenarios:\n")
      print(x$H1$scenario,row.names=FALSE)
      cat(" where 0 = no stop, 1 = efficacy stop, 2 = futility stop")
      if(any(x$H1$scenario[,x$H1$target$par$id]==3)){
        cat(",\n       3 = simultaneous efficacy and futility stops")
      }else{cat("\n")}
    }
  }
  cli_h1("")#cat(paste0(rep("-",options()$width),collapse=""))
  cat("\n")
}

