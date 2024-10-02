#' @name summary.batss
#' @title Summary function for 'BATSS' outputs
#' @description Summary method function for objects of class 'batss'.
#' @param object An object of class 'batss' (i.e., output of the function [batss.glm]).
#' @param full A logical indicating if a standard (full = FALSE, default) or extended output (full = TRUE) should be returned.
#' @param ... Additional arguments affecting the summary produced.
#' @returns Prints a summary for objects of class 'batss'.
#' @seealso [batss.glm()], the function generating S3 objects of class 'batss'. 
#' @export
summary.batss = function(object, full=FALSE, ...){
    # common part
    cat("\n")
    if(!is.null(object$par$RAR)){
      cli_h1("Bayesian Adaptive Design with Laplace Approx.")
    }else{
      cli_h1("MAMS with Laplace Approx.")
    }
    cat("  (",length(object$par$seed)," Monte Carlo samples)\n",sep="")
    cat("\n")
    cli_h3("Variables:")
    for(i in 2:length(object$call$var)){
      cat("  *",names(object$call$var)[i],":",as.character(object$call$var)[i],"\n")
    }
    if(!is.null(object$par$RAR)){
        cat("\n")
        cli_h3("Group randomisation:")
        cat("  *",object$call$RAR,"\n")
    }
    cat("\n")
    cli_h3("Model: ")
    cat("  *",format(object$call$model),if (object$type=="surv") "\n" else paste("(with",ifelse(is.null(object$call$link),"identity",object$call$link), "link)\n"))
    cat("\n")
    cli_h3("Fixed effect parameters:\n")
    if (object$type=="surv") objectw = object$hr else objectw = object$beta
    colnames(objectw)[1] = ""
    print(objectw,row.names=FALSE)
    if (object$type!="surv"){
      cat("\n")
      cli_h3("Sample size per interim analyis:\n")
      objectw = object$look
      colnames(objectw)[1] = ""
      print(objectw,row.names=FALSE)
    }
    # H0
    if(object$par$H0){
      cat("\n\n")
      #cat(paste0(rep("-",floor(options()$width/2)),collapse=""))
      cli_h2("\n H0: Under the null hypothesis\n")
      #cat(paste0(rep("-",floor(options()$width/2)),collapse=""))
      #
      if (object$type=="surv") {
        cat("\n")
        cli_h3("Average sample size per interim analyis:\n")
        objectw = round(rowMeans(sapply(object$H0$trial,function(x) as.matrix(x$look[,c("n","m","t","t(n)","ev(n)")]),simplify="array"),dims=2,na.rm=TRUE),2)
        print(objectw,row.names=FALSE)
      }
      cat("\n")
      cli_h3("Target parameters:\n")
      temp = rbind(object$H0$target$par,object$H0$target$global)
      if(all(temp$both==0)){temp = temp[,colnames(temp)!="both"]}       
      temp[,1][-(1:nrow(object$H0$target$par))] = ""
      print(temp,row.names=FALSE)
      if (object$type=="surv") {
        cat("\n")
        cli_h3(paste0( "Trial duration:\n"))
        print(summary(object$H0$sample[,"t"]))
      }
      #
      if(full){
          cat("\n")
          cli_h3("Efficacy:\n")
          temp = rbind(object$H0$efficacy$par,object$H0$efficacy$global)
          temp[,1][-(1:nrow(object$H0$efficacy$par))] = ""
          print(temp,row.names=FALSE)
          #
          cat("\n")
          cli_h3("Futility:\n")
          temp = rbind(object$H0$futility$par,object$H0$futility$global)
          temp[,1][-(1:nrow(object$H0$futility$par))] = ""
          print(temp,row.names=FALSE)
          #
          cat("\n")
          cli_h3("Sample size per group:\n")
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
          print(temp,row.names=FALSE)
          #cat("\n")
          if (object$type=="surv") {
            cli_h3("Events per group:\n")
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
            print(temp,row.names=FALSE)
            #cat("\n")
            cli_h3("Time observed per group:\n")
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
            print(temp,row.names=FALSE)
          }
          #
          cat("\n")
          cli_h3("Scenarios:\n")
          print(object$H0$scenario,row.names=FALSE)
          cat(" where 0 = no stop, 1 = efficacy stop, 2 = futility stop\n")
          if(any(object$H0$scenario[,object$H0$target$par$id]==3)){
            cat(",\n       3 = simultaneous efficacy and futility stops")
          }else{cat("\n")}
          }
    }
    # H1
    if(object$par$H1){
      cat("\n\n")
      #cat(paste0(rep("-",floor(options()$width/2)),collapse=""))
      cli_h2("\n H1: Under the alternative hypothesis\n")#cat("\n H1: Under the alternative hypothesis\n")
      #cat(paste0(rep("-",floor(options()$width/2)),collapse=""))
      #
      if (object$type=="surv") {
        cat("\n")
        cli_h3("Average sample size per interim analyis:\n")
        objectw = round(rowMeans(sapply(object$H1$trial,function(x) as.matrix(x$look[,c("n","m","t","t(n)","ev(n)")]),simplify="array"),dims=2,na.rm=TRUE),2)
        print(objectw,row.names=FALSE)
      }
      cat("\n")
      cli_h3("Target parameters:\n")
      temp = rbind(object$H1$target$par,object$H1$target$global)
      if(all(temp$both==0)){temp = temp[,colnames(temp)!="both"]} 
      temp[,1][-(1:nrow(object$H1$target$par))] = ""
      print(temp,row.names=FALSE)
      if (object$type=="surv") {
        cat("\n")
        cli_h3(paste0( "Trial duration:\n"))
        print(summary(object$H1$sample[,"t"]))
      }
      #
      if(full){
          cat("\n")
          cli_h3("Efficacy:\n")
          temp = rbind(object$H1$efficacy$par,object$H1$efficacy$global)
          temp[,1][-(1:nrow(object$H1$efficacy$par))] = ""
          print(temp,row.names=FALSE)
          #
          cat("\n")
          cli_h3("Futility:\n")
          temp = rbind(object$H1$futility$par,object$H1$futility$global)
          temp[,1][-(1:nrow(object$H1$futility$par))] = ""
          print(temp,row.names=FALSE)
          #
          cat("\n")
          cli_h3("Sample size per group:\n")
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
          print(temp,row.names=FALSE)
          #cat("\n")
          if (object$type=="surv") {
            cli_h3("Events per group:\n")
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
            print(temp,row.names=FALSE)
            #cat("\n")
            cli_h3("Time observed per group:\n")
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
            print(temp,row.names=FALSE)
          }
          #
          cat("\n")
          cli_h3("Scenarios:\n")
          print(object$H1$scenario,row.names=FALSE)
          cat(" where 0 = no stop, 1 = efficacy stop, 2 = futility stop")
          if(any(object$H1$scenario[,object$H1$target$par$id]==3)){
            cat(",\n       3 = simultaneous efficacy and futility stops")
          }else{cat("\n")}
          }
    }
    cli_h1("")#cat(paste0(rep("-",options()$width),collapse=""))
    cat("\n")
}
