#' @name summary.batss
#' @title Summary function for 'BATSS' outputs
#' @description Summary method function for objects of class 'batss'.
#' @param object An object of class 'batss' (i.e., output of the function [batss.glm]).
#' @param extended A logical indicating if a standard (extended = FALSE, default) or extended output (extended = TRUE) should be returned. Default to `NULL` in which case the input of the argument `extended` chosen when generating `object` with [batss.glm()] or [batss.surv()] is used.
#' @param ... For future use
#' @returns Object of class 'summary.batss'.
#' @returns The function [summary.batss] returns an S3 list of class 'summary.batss' with available print functions. The list elements are
#' \itemize{
#'   \item beta - A data frame providing information related to the beta parameter vector, such as parameter names and values, for example.
#'   \item look - A data frame providing information related to looks, like sample size of a given interim (m) and cumulative sample size at a given interim (n), for example.
#'   \item par - A list providing different information, like the used seeds (seed) and the groups (group), for example.
#'   \item H1 - A list providing trial aggregated results under the alternative, like the probability of efficacy, futility, or both, per arm or globally (`object$H1$target`), the probability of stopping early for efficacy (`object$H1$efficacy`) and futility (`object$H1$futility`), the sample size expectation, standard deviation, and quantiles 0.1, 0.5 and 0.9, per group and overall (`object$H1$summary.sample.sizes`), the probabilities associated to each combination of efficacy and futility per group (scenario).
#'   \item H0 - A list providing trial aggregated results under the global null hypothesis (same structure as H1).
#'   \item call - The matched call.
#'   \item type - The type of 'BATSS' analysis (currently either 'glm' or 'surv').
#' }
#' @seealso [batss.glm()], [batss.surv()], the functions generating S3 objects of class 'batss'. 
#' @export
summary.batss = function(object, extended=NULL, ...){
  res          <- list()
  res$call     <- object$call
  res$par      <- object$par
  res$type     <- object$type
  if (object$type=="surv") {
    res$hr     <- object$hr
  } else {
    res$beta   <- object$beta
    res$sample <- object$look
    colnames(res$sample)[1] = ""
  }
  if(is.null(extended)){
    res$extended <- ifelse(is.null(object$call$extended),FALSE,object$call$extended>0)
  }else{
    res$extended <- as.numeric(extended) 
  }
  if (object$par$H0) {
    
    # sample size
    res$H0$sample.sizes <- object$H0$sample[,1:ifelse(object$type=="surv",(sum(object$hr$target)+1),sum(object$beta$target))]
    
    if (object$type=="surv") {
      # events
      res$H0$events <- object$H0$sample[,(2*(sum(object$hr$target)+1)+1):(3*(sum(object$hr$target)+1))]
      # time observed
      res$H0$t.obs <- object$H0$sample[,(sum(object$hr$target)+2):(2*(sum(object$hr$target)+1))]
      # trial duration 
      res$H0$t.trial <- object$H0$sample[,"t"]
    }
    
    # target parameters
    res$H0$target <- object$H0$target
    temp = rbind(object$H0$target$par,object$H0$target$global)
    if(all(temp$both==0)){temp = temp[,colnames(temp)!="both"]}       
    res$H0$target <- temp
    
    # efficacy
    temp = rbind(object$H0$efficacy$par,object$H0$efficacy$global)
    temp[,1][-(1:nrow(object$H0$efficacy$par))] = ""
    res$H0$efficacy <- temp
    
    # futility
    temp = rbind(object$H0$futility$par,object$H0$futility$global)
    temp[,1][-(1:nrow(object$H0$futility$par))] = ""
    res$H0$futility <- temp
    
    #cumulative
    tmp.length <- ifelse(is.null(dim(object$H0$estimate[,,1])),1,dim(object$H0$estimate[,,1])[1])
    cum.eff <- cum.fut <- matrix(NA,nrow=tmp.length+2,ncol=dim(object$look)[1])
    for (i in 1:dim(object$look)[1]){
      #arms
      cum.eff[1:tmp.length,i] <- apply((object$H0$estimate[,"type",,drop=FALSE]==1)*(object$H0$estimate[,"look",,drop=FALSE]<=i),1,mean)
      cum.fut[1:tmp.length,i] <- apply((object$H0$estimate[,"type",,drop=FALSE]==2)*(object$H0$estimate[,"look",,drop=FALSE]<=i),1,mean)
      #any
      cum.eff[tmp.length+1,i] <- mean(apply(object$H0$estimate[,"type",,drop=FALSE]==1&(object$H0$estimate[,"look",,drop=FALSE]<=i),3,sum)>0)
      cum.fut[tmp.length+1,i] <- mean(apply(object$H0$estimate[,"type",,drop=FALSE]==2&(object$H0$estimate[,"look",,drop=FALSE]<=i),3,sum)>0)
      #all
      cum.eff[tmp.length+2,i] <- mean(apply(object$H0$estimate[,"type",,drop=FALSE]==1&(object$H0$estimate[,"look",,drop=FALSE]<=i),3,all)>0)
      cum.fut[tmp.length+2,i] <- mean(apply(object$H0$estimate[,"type",,drop=FALSE]==2&(object$H0$estimate[,"look",,drop=FALSE]<=i),3,all)>0)
    }
    rownames(cum.eff) <- rownames(cum.fut) <- c(object$H0$efficacy$par$group,"At least one","All")
    colnames(cum.eff) <- colnames(cum.fut) <- c(paste("Look", 1:(dim(object$look)[1]-1)), "   Final")
    res$H0$cum.eff <- cum.eff
    res$H0$cum.fut <- cum.fut

    if(res$extended>0){
      # ss
      temp_col = nrow(object$par$group)
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
        # summary events
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
        
        # summary time observed
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
      
      # scenarios
      res$H0$scenario <- object$H0$scenario
      
      if (object$type=="surv") {
        # average sample size at interims
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
    }
    
    # estimation error
    if (object$type=="surv") {
      if (!is.null(dim(object$H0$estimate[,3,]))) {
        res$H0$est.error <- t(sweep(exp(object$H0$estimate[,3,]),1,object$hr[object$hr$target,4]))
      } else {
        res$H0$est.error <- exp(object$H0$estimate[,3,]) - object$hr[object$hr$target,4]
      }
    } else {
      if (!is.null(dim(object$H0$estimate[,3,]))) {
        res$H0$est.error <- t(sweep(object$H0$estimate[,3,],1,object$beta[object$beta$target,4]))
      } else {
        res$H0$est.error <- object$H0$estimate[,3,] - object$beta[object$beta$target,4]
      }
    }
    
  } else {
    res$H0 <- NULL 
  } 
  
  if(object$par$H1){
    #sample size
    res$H1$sample.sizes <- object$H1$sample[,1:ifelse(object$type=="surv",(sum(object$hr$target)+1),sum(object$beta$target))]
    
    if (object$type=="surv") {
      # events
      res$H1$events <- object$H1$sample[,(2*(sum(object$hr$target)+1)+1):(3*(sum(object$hr$target)+1))]
      # time observed
      res$H1$t.obs <- object$H1$sample[,(sum(object$hr$target)+2):(2*(sum(object$hr$target)+1))]
      # trial duration
      res$H1$t.trial <- object$H1$sample[,"t"]
    }
    
    # target parameters
    res$H1$target <- object$H1$target
    temp = rbind(object$H1$target$par,object$H1$target$global)
    if(all(temp$both==0)){temp = temp[,colnames(temp)!="both"]}       
    res$H1$target <- temp
    
    # efficacy
    temp = rbind(object$H1$efficacy$par,object$H1$efficacy$global)
    temp[,1][-(1:nrow(object$H1$efficacy$par))] = ""
    res$H1$efficacy <- temp
    
    # futility
    temp = rbind(object$H1$futility$par,object$H1$futility$global)
    temp[,1][-(1:nrow(object$H1$futility$par))] = ""
    res$H1$futility <- temp
    
    #cumulative
    tmp.length <- ifelse(is.null(dim(object$H1$estimate[,,1])),1,dim(object$H1$estimate[,,1])[1])
    cum.eff <- cum.fut <- matrix(NA,nrow=tmp.length+2,ncol=dim(object$look)[1])
    for (i in 1:dim(object$look)[1]){
      #arms
      cum.eff[1:tmp.length,i] <- apply((object$H1$estimate[,"type",,drop=FALSE]==1)*(object$H1$estimate[,"look",,drop=FALSE]<=i),1,mean)
      cum.fut[1:tmp.length,i] <- apply((object$H1$estimate[,"type",,drop=FALSE]==2)*(object$H1$estimate[,"look",,drop=FALSE]<=i),1,mean)
      #any
      cum.eff[tmp.length+1,i] <- mean(apply(object$H1$estimate[,"type",,drop=FALSE]==1&(object$H1$estimate[,"look",,drop=FALSE]<=i),3,sum)>0)
      cum.fut[tmp.length+1,i] <- mean(apply(object$H1$estimate[,"type",,drop=FALSE]==2&(object$H1$estimate[,"look",,drop=FALSE]<=i),3,sum)>0)
      #all
      cum.eff[tmp.length+2,i] <- mean(apply(object$H1$estimate[,"type",,drop=FALSE]==1&(object$H1$estimate[,"look",,drop=FALSE]<=i),3,all)>0)
      cum.fut[tmp.length+2,i] <- mean(apply(object$H1$estimate[,"type",,drop=FALSE]==2&(object$H1$estimate[,"look",,drop=FALSE]<=i),3,all)>0)
    }
    rownames(cum.eff) <- rownames(cum.fut) <- c(object$H1$efficacy$par$group,"At least one","All")
    colnames(cum.eff) <- colnames(cum.fut) <- c(paste("look", 1:(dim(object$look)[1]-1)), "   Final")
    res$H1$cum.eff <- cum.eff
    res$H1$cum.fut <- cum.fut
    
    if (res$extended>0){
      # ss
      temp_col = nrow(object$par$group)
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
        # summary events
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
        
        #summary time observed
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
      
      # scenarios
      res$H1$scenario <- object$H1$scenario
      
      if (object$type=="surv") {
        # average sample size at interims
        if (is.null(object$par$interim$event) || (!is.null(object$par$interim$event.type) && object$par$interim$event.type!="cplusone"))  {
          objectw = rowMeans(sapply(object$H1$trial,function(x) as.matrix(cbind(x$look[,c("n","m","t","t(n)","ev(n)")],1-(is.na(x$look[,c("t")])&!is.na(x$look[,c("pos")])),is.na(x$look[,c("pos")]))),
                                    simplify="array"),dims=2,na.rm=TRUE)
          objectw = cbind(objectw[,-(6:7)],abs(diff(c(objectw[,6]*length(object$par$seed),0))),objectw[,6]*length(object$par$seed),objectw[,6],objectw[,7]*length(object$par$seed),objectw[,7])
          colnames(objectw)[6:10] <- c("stops","n.trial","prcnt.trial","n.skipped","prcnt.skipped")
        } else {
          objectw <- NULL
        }
        res$H1$avg.sample <- objectw
      }
    }
    
    # estimation error
    if (object$type=="surv") {
      if (!is.null(dim(object$H1$estimate[,3,]))) {
        res$H1$est.error <- if (object$par$H0 && object$par$H1) t(sweep(exp(object$H1$estimate[,3,]),1,object$hr[object$hr$target,5])) else t(sweep(exp(object$H1$estimate[,3,]),1,object$hr[object$hr$target,4]))
      } else {
        res$H1$est.error <- if (object$par$H0 && object$par$H1) exp(object$H1$estimate[,3,]) - object$hr[object$hr$target,5] else object$H1$estimate[,3,] - object$hr[object$hr$target,4]
      }
    } else {
      if (!is.null(dim(object$H1$estimate[,3,]))) {
        res$H1$est.error <- if (object$par$H0 && object$par$H1) t(sweep(object$H1$estimate[,3,],1,object$beta[object$beta$target,5])) else t(sweep(object$H1$estimate[,3,],1,object$beta[object$beta$target,4]))
      } else {
        res$H1$est.error <- if (object$par$H0 && object$par$H1) object$H1$estimate[,3,] - object$beta[object$beta$target,5] else object$H1$estimate[,3,] - object$beta[object$beta$target,4]
      }
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
#' @param ... For future use
#' @seealso [batss.glm()], [batss.surv()], the functions generating S3 objects of class 'batss'. 
#' @export 
print.summary.batss = function(x, ...){
  # common part
  if(!is.null(x$par$RAR)){
    cli_h1("Bayesian Adaptive Design with Laplace Approx.")
  }else{
    cli_h1("MAMS with Laplace Approx.")
  }
  cat("  (",length(x$par$seed)," Monte Carlo samples)\n",sep="")
  cli_h3("Variables:")
  for(i in 2:length(x$call$var)){
    cat("  *",names(x$call$var)[i],":",as.character(x$call$var)[i],"\n")
  }
  if(!is.null(x$par$RAR)){
    cli_h3("Group randomisation:")
    cat("  *",x$call$RAR,"\n")
  }
  if(!(is.null(x$par$eff.arm) && is.null(x$par$fut.arm))){
    cli_h3("Decision rules:")
    if(!is.null(x$par$eff.arm)) cat("  * Efficacy: ",format(x$call$eff.arm),"\n")
    if(!is.null(x$par$fut.arm)) cat("  * Futility: ",format(x$call$fut.arm),"\n")
  }
  cli_h3("Model: ")
  cat("  *",format(x$call$model),if (x$type=="surv") "\n" else paste("(with",ifelse(is.null(x$call$link),"identity",x$call$link), "link)\n"))
  cli_h3("Fixed effect parameters:\n")
  if (x$type=="surv") xw = x$hr else xw = x$beta
  colnames(xw)[1] = ""
  print(xw,row.names=FALSE)
  if (x$type!="surv"){
    cli_h3("Sample size per interim analyis:\n")
    print(x$sample,row.names=FALSE)
  }
  # H0
  if(x$par$H0){
    cat("\n")
    #cat(paste0(rep("-",floor(options()$width/2)),collapse=""))
    cli_h2("\n H0: Under the null hypothesis\n")
    #cat(paste0(rep("-",floor(options()$width/2)),collapse=""))
    #
    if (x$type=="surv" && x$extended>0) {
      cli_h3("Average sample size per interim analyis:\n")
      print(x$H0$avg.sample,row.names=FALSE,digits=4)
    }
    cli_h3("Target parameters:\n")
    print(x$H0$target,row.names=FALSE)
    if (x$type=="surv") {
      cli_h3(paste0( "Trial duration:\n"))
      print(summary(x$H0$t.trial))
    }
    #
    if(x$extended>0){
      cli_h3("Efficacy:\n")
      print(x$H0$efficacy,row.names=FALSE,digits=4)
      #
      cli_h3("Cumulative efficacy:\n")
      print(x$H0$cum.eff,row.names=FALSE,digits=4)
      #
      if (!is.null(x$par$fut.arm)){
        cli_h3("Futility:\n")
        print(x$H0$futility,row.names=FALSE,digits=4)
        #
        cli_h3("Cumulative futility:\n")
        print(x$H0$cum.fut,row.names=FALSE,digits=4)
      }
      #
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
    #cat(paste0(rep("-",floor(options()$width/2)),collapse=""))
    cli_h2("\n H1: Under the alternative hypothesis\n")#cat("\n H1: Under the alternative hypothesis\n")
    #cat(paste0(rep("-",floor(options()$width/2)),collapse=""))
    #
    if (x$type=="surv" && x$extended>0) {
      cli_h3("Average sample size per interim analyis:\n")
      print(x$H1$avg.sample,row.names=FALSE,digits=4)
    }
    cli_h3("Target parameters:\n")
    print(x$H1$target,row.names=FALSE)
    if (x$type=="surv") {
      cli_h3(paste0( "Trial duration:\n"))
      print(summary(x$H1$t.trial))
    }
    #
    if(x$extended>0){
      cli_h3("Efficacy:\n")
      print(x$H1$efficacy,row.names=FALSE,digits=4)
      #
      cli_h3("Cumulative efficacy:\n")
      print(x$H1$cum.eff,row.names=FALSE,digits=4)
      #
      if (!is.null(x$par$fut.arm)){
        cli_h3("Futility:\n")
        print(x$H1$futility,row.names=FALSE,digits=4)
        #
        cli_h3("Cumulative futility:\n")
        print(x$H1$cum.fut,row.names=FALSE,digits=4)
      }
      #
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

