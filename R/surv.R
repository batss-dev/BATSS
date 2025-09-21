#' @name batss.surv
#' @title Bayesian adaptive trial simulations for time-to-event data
#' @description Simulation of Bayesian adaptive trials with time-to-event endpoint using Integrated Nested Laplace Approximation (INLA).  
#' @param model an object of class '\link[stats]{formula}' indicating a symbolic description of the model to be fitted
#' @param surv Time-to-event data simulator function. Currently, '\link[simsurv]{simsurv}' is the only option and used by default.
#' @param surv.control A list of parameters for the call to '\link[simsurv]{simsurv}'.
#' @param var A list. Each entry corresponds to a variable described under '`model`' and indicates the name of a function allowing to generate variates (like \link[stats]{rnorm} and \link[stats]{rexp}, for example). The list names have to match the variable names used in '`model`'. The grouping variable corresponding to the target parameters has to be of class '\link[base]{factor}' with levels corresponding to the names indicated in argument `prob0` (see below).
#' @param var.control An optional list of control parameters for the functions indicated in '`var`'. The names of the list items need to correspond to the names used in '`var`'. Each element is another list with names of the elements corresponding to the parameter names of the functions specified in '`var`'. 
#' @param family A character string indicating the log-survival likelihood as described in the package INLA. Default set to '`weibullsurv`', other options include '`exponential.surv`', '`loglogistic.surv`','`lognormal.surv`', '`coxph`'. See \link[INLA]{inla.list.models}.
#' @param hr A numerical vector of parameter values for the hazard ratios. Its length has to match the number of treatments (excluding control).
#' @param which A numerical vector indicating the position of the target parameters.
#' @param cens Random number generating function for censoring, e.g. '\link[stats]{rexp}'
#' @param cens.control Additional parameters for call to censoring function.
#' @param accr Random number generating function for entry times, e.g. '\link[stats]{rexp}'
#' @param accr.control Additional parameters for call to entry times generating function.
#' @param accr.type One of 'random' or 'fixed', the former will draw entry times from the distribution defined in '`accr`', adding the resulting times results in a random study duration. Setting this parameter to 'fixed' will allow setting a fixed study length and '`accr`' will be used to describe the distribution within the interval from '`0`' to study duration.
#' @param R a vector of natural numbers to be used as seeds (check \link[base]{set.seed}) for the different Monte Carlo trials (the vector length will thus correspond to the number of Monte Carlo trials). When `R` is a scalar, seeds `1` to `R` are used, where `R` corresponds to the number of Monte Carlo trials. 
#' @param alternative A vector of strings providing the one-sided direction of the alternative hypothesis corresponding to each target parameter indicated under '`which`' (in the same order). Possibilities are 'greater' or 'less'. If the vector is of length 1, the same direction will be used for all target parameter tests.
#' @param RAR A character string indicating the how the response-adaptive randomised probabilities corresponding to each group, reference group included, are defined. If RAR=NULL, equal probabilities of being attributed to each (active) group are used. Check ?BATS::prob.fun for examples.
#' @param RAR.control An optional list of control parameters the function in 'RAR'.
#' @param N A scalar indicating the total sample size.
#' @param n.max A vector indicating the maximum sample size per arm, if a scalar is given the same size will be used for all arms.
#' @param interim A list of parameters related to interim analyses. Possible list items include, '`recruited`' a vector of integers indicating the number of recruited participants at each look, last excluded, in increasing order, '`time`' a vector of integers indicating the time points of interim analyses, and '`event`' a vector of integers indicating the number of events. '`event.type`' is a character string describing the groups that are counted towards the numbers in '`event`' with options options '`all`', counting all events,'`control`', counting events in the control group only, '`cplusone`', counting the control plus treatments individually (this will result in different time points for the interims of the treatment groups), and '`min`', counting the last group to reach the threshold. See details for explanations and combinations.
#' @param fup A scalar specifying the follow-up time of the last patient recruited and hence the final analysis, this only applies if 'maxt=NULL' in 'surv.control'.
#' @param prob0 A named vector with initial allocation probabilities. Names need to correspond to the levels of the grouping variable. If `RAR = NULL`, these probabilities/ratios will be used throughout (fixed allocation probabilities).
#' @param delta.eff A vector (of length equal to the number of looks (i.e., number of interims + 1)) of clinically meaningful treatment effect values (on the linear predictor scale) to be used to define the efficacy-related posterior probabilities for each target parameter at each look. If a scalar is provided, the same value is used at each look. The default is `delta.eff = 0`. 
#' @param delta.fut A vector (of length equal to the number of looks (i.e., number of interims + 1)) of clinically meaningful treatment effect values (on the linear predictor scale) to be used to define the futility-related posterior probabilities for each target parameter at each look. If a scalar is provided, the same value is used at each look. The default is `delta.fut = delta.eff`. 
#' @param delta.RAR A vector (of length equal to the number of looks (i.e., number of interims + 1)) of clinically meaningful treatment effect values (on the linear predictor scale) to be used to define the RAR-related posterior probabilities for each target parameter at each look. If a scalar is provided, the same value is used at each interim analysis. The default is `delta.RAR = 0`. Note that, when a vector is provided, its last value is ignored as no randomisation is made at the last look.
#' @param eff.arm A function defining if efficacy has been achieved at a given look given the information available at that stage a given target parameter. The output of this function must be a \link[base]{logical} (of length 1). Arguments of this function will typically consider 'BATSS' ingredients. Check [eff.arm.simple] and [eff.arm.infofract] for examples. 
#' @param eff.arm.control An optional list of parameters for the function indicated in '`eff.arm`'.
#' @param eff.trial A function defining if the trial can be stopped for efficacy given the output of the function indicated in '`eff.arm`'. The output of this function must be a \link[base]{logical} of length one. Arguments of this function will typically only consider the 'BATSS' ingredient `eff.target`. Check [eff.trial.all] and [eff.trial.any] for examples. When `eff.trial = NULL` (default), the trial stops for efficacy when *all* target parameters are found to be effective (like in [eff.trial.all]). 
#' @param eff.trial.control An optional list of parameters for the function indicated in '`eff.trial`'.
#' @param fut.arm A function defining if futility has been achieved at a given look given the information available at that stage for each target parameter. The output of this function must be a \link[base]{logical} (of length 1). Arguments of this function will typically consider 'BATSS' ingredients. Check [fut.arm.simple] to see an example of such a function. 
#' @param fut.arm.control An optional list of parameters for the function indicated in '`fut.arm`'.
#' @param fut.trial A function defining if the trial can be stopped for futility given the output of the function indicated in '`fut.arm`'. The output of this function must be a \link[base]{logical} of length one. Arguments of this function will typically only consider the 'BATSS' ingredient `fut.target`. Check [fut.trial.all] for an example of such a function. When `fut.trial = NULL` (default), the trial stops for futility when *all* target parameters are found to be futile (like in [fut.trial.all]).
#' @param fut.trial.control An optional list of parameters for the function indicated in '`fut.trial`'.
#' @param H0 A logical indicating whether the simulation should also consider the case with all target parameters set to 0 to check the probability of rejecting the hypothesis that the target parameter value is equal to 0 individually (pairwise type I error) or globally (family-wise error rate). Default set to `H0=TRUE`.
#' @param computation A character string indicating how the computation should be performed. Possibilities are '`parallel`' or '`sequential`' with default `computation="parallel"` meaning that the computation is split between `mc.cores`. 
#' @param mc.cores An integer indicating the number of CPUs to be used when `computation="parallel"` (Default to 3 if no global '`mc.cores`' global option is available via \link[base]{getOption}).
#' @param extended an integer indicating the type of results to be returned. 0 (default) provides summary statistics, 1 adds the results of each Monte Carlo trial and 2 additionally returns each Monte Carlo dataset. [batss.combine] requires extended > 0 as the function needs to merge results of different sets of seeds.
#' @param ... Additional arguments to control fitting in \link[INLA]{inla}.
#' @returns The function [batss.surv] returns an S3 object of class 'batss' with available print/summary/plot functions
#' \itemize{
#'   \item beta - A data frame providing information related to the beta parameter vector, like parameter names and values, for example.
#'   \item look - A data frame providing information related to looks, like sample size of a given interim (m) and cumulative sample size at a given interim (n), for example.
#'   \item par - A list providing different information, like the used seeds (seed) and the groups (group), for example.
#'   \item H1 - A list providing trial results under the alternative, like the estimates per target parameter when the corresponding arm was stopped (estimate), the efficacy and futility probabilites per target parameter and overall (target, efficacy and futility), the sample size per group and trial (sample), the probabilities associated to each combination of efficacy and futility per group (scenario), the detailed results per trial (trial), for example.
#'   \item H0 - A list providing trial results under the global null hypothesis (same structure as H1).
#'   \item call - The matched call.
#'   \item type - The type of 'BATSS' analysis, i.e. '`surv`'
#' }
#' @export
#' @seealso [summary.batss] and [plot.batss] for detailed summaries and plots, and [batss.combine] to combine different evaluations of [batss.surv] considering the same trial design but different sets of seeds (useful for cluster computation). 
#' @examples
#' # Example 1: TBC
#' @export
batss.surv = function(
    model,family="weibullsurv",surv=simsurv::simsurv,surv.control,var,var.control=NULL, fup=NULL,
    cens=NULL,cens.control=NULL,accr=runif,accr.control=list(min = 0, max = 3), accr.type="fixed",
    hr,which,R=1e+4,N,n.max=NULL,
    alternative = "less",RAR=NULL,RAR.control=NULL,
    interim,
    prob0,delta.eff=0,delta.fut=delta.eff,delta.RAR=0,
    eff.arm,eff.trial=NULL,
    eff.arm.control=NULL,eff.trial.control=NULL,
    fut.arm,fut.trial=NULL,
    fut.arm.control=NULL,fut.trial.control=NULL,
    H0=TRUE,computation="parallel",
    mc.cores=getOption("mc.cores", 3L),
    extended = 0,...){
  
  #---
#  require(plyr); require(rlang); require(R.utils); require(simsurv); require(foreach); require(INLA)
  
  call <- match.call()                       #save call
  model <- as.formula(model)                 #allow for string and formula input
  
  #---

  ##
  ## dataset structure and useful definitions
  ##
 
  #some checks
  n = m = prob = NULL
  #error messages
  if (!is.character(model) && !is.formula(model))
    stop("invalid 'model' argument")
  if (!is.list(var) || !all(sapply(var,is.function)))
    stop("'var' must be a list of functions")
  if (is.null(intersect(names(var),intersect(setdiff(unlist(strsplit(attr(terms(model),"term.labels"),"1")),attr(terms(model),"term.labels")),
                                             setdiff(unlist(strsplit(attr(terms(model),"term.labels"),"2")),attr(terms(model),"term.labels"))))) &&
      !setequal(attr(terms(model),"term.labels"),names(var)))
    stop("all variables on the right side of the model formula must have a generating function in the 'var' list")
  if (!(family %in% names(INLA::inla.models()$likelihood)))
    stop("invalid 'family' argument, see help files and inla documentation for available families")
  if (!is.null(interim$recruited) && !is.numeric(unlist(interim$recruited)))
    stop("'interim$recruited' must be a (list of) numeric vector(s)")
  if (!is.null(interim$time) && any(na.omit(interim$time) < 0))
    stop("negative times for interims not allowed")
  if (!is.null(interim$event) && any(interim$event < 0))
    stop("negative interim event numbers not allowed")
  if (!is.null(interim$recruited) && any(interim$recruited < 0))
    stop("negative interim recruitment numbers not allowed")
  if ((N <= 0) || length(N)>1)
    stop("total sample size 'N' must be a positive scalar")
  if (!is.null(n.max) && (any(n.max <= 0) || length(n.max)!=length(prob0)))
    stop("vector of maximum sample sizes must be positive and the same length as 'prob0'")
  if (length(which)>length(hr))
    stop("number of targets greater than number of parameters")
  if ((!is.null(RAR) && !(is.function(RAR))) ||
      (!is.null(eff.arm) && !(is.function(eff.arm))) || (!is.null(eff.trial) && !(is.function(eff.trial))) ||
      (!is.null(fut.arm) && !(is.function(fut.arm))) || (!is.null(fut.trial) && !(is.function(fut.trial))))
    stop("'RAR', 'eff.arm', 'eff.trial', 'fut.arm', 'fut.trial' must be functions or NULL")
  if (!is.null(interim$event)) {
      if (is.null(interim$event.type)){
        stop("'event.type' must be specified if interim based on events")
      } else if (interim$event.type=="cplusone" && !is.null(RAR)) {
        stop("RAR and event-based interims based on control plus one intervention group are not compatible")
      }
  }
  
  #warnings
  if (is.null(interim$recruited) && is.null(interim$time) && is.null(interim$event))
    warning("no interim analyses specified")
  if (all(prob0<0) && sum(prob0)!=1)
    warning("sum of 'prob0' not equal to 1")
  if (!is.null(interim$event) && (interim$event.type %in% c("min","control")) && !is.null(n.max))
    warning("low intervention maximum sample sizes may lead to error as event numbers may not be reached in all (active) groups")
  
  # size per look
  if (!is.null(interim$recruited)){
    if(sum(!is.na(match(c("m0","m"),names(interim$recruited))))==2){
      size_look = seq(interim$recruited$m0,N,interim$recruited$m)
      if(any(size_look > N)) {
        size_look <- c(size_look[size_look<N],N)
        warning("some interim analyses are outside the maximum sample size and will be ignored")
        }
    }else{
      if (any(interim$recruited > N)) warning("some interim analyses are outside the maximum sample size and will be ignored")
      size_look = interim$recruited[interim$recruited<=N]
    }
  }
  
  message("    Initialisation") 
  
  #number of looks
  n.look  = ifelse(is.null(interim$time),
                   ifelse(is.null(interim$event),length(size_look)+1,
                          ifelse(interim$event.type=="cplusone",length(interim$event)*(length(prob0)-1)+1,length(interim$event)+1)),
                   length(interim$time)+1)
  
  id.look = data.frame(pos = 1:n.look,
                       id  = if (!is.null(interim$time) && !is.na(interim$time[1])) {c(paste0("t=",interim$time),"final")}
                       else {
                         if (!is.null(interim$event)) {
                           if (!is.null(interim$time)) {
                             c(paste0("n(e)=",interim$event[1]),paste0("t=t(e=",interim$event[1],")+",interim$time[-1]),"final")
                           } else {
                             if (interim$event.type!="cplusone") {
                               c(paste0("n(e)=",interim$event),"final")
                             } else {
                               c(paste0("n(e)=",rep(interim$event,length(prob0)-1)),"final")
                             }
                           }
                         } else {c(paste0("n=",size_look),paste0("n=",N," (final)"))}},
                       n   = if (!is.null(interim$time) | !is.null(interim$event)) rep(N,n.look) else c(size_look,N),
                       m   = if (!is.null(interim$time) | !is.null(interim$event)) c(N,rep(0,n.look-1)) else c(size_look[1],size_look[-1]-size_look[-(n.look-1)],N-size_look[n.look-1]))

  # generate predictors
  env0 = new.env()
  n.var  = length(labels(terms(model)))
  id.var = names(var)
  assign("m", ifelse(!is.null(interim$time),N,id.look$m[1]), envir = env0)
  assign("prob",prob0,envir = env0)
  assign("var.control", var.control, envir = env0)
  
  covar <- vector("list",length(var))
  for (ii in 1:length(var)) {
    tmp_nam <- names(var)[ii]
    args_ <- plyr::.(n=n,m=m,prob=prob)
    if (tmp_nam %in% names(var.control)) args_ <- c(args_, var.control[[tmp_nam]])
    covar[[ii]] <- R.utils::doCall(var[[ii]], envir = env0, args = args_)   #call functions directly
    if (is.matrix(covar[[ii]])) colnames(covar[[ii]]) <- paste0(tmp_nam,1:dim(covar[[ii]])[2])
  }
  
  if (any(sapply(covar,is.matrix))) {
    where.mat <- which(sapply(covar,is.matrix))
    tmp_var <- covar[where.mat]
    id.var <- names(var)[1]
    for (ii in 1:length(covar)) {
      if (ii %in% (where.mat)) {
        id.var <- c(id.var,colnames(covar[[ii]]))
      } else {
        id.var <- c(id.var,names(var[ii]))
      }
    }
  }
  
  m0     = length(covar[[1]])
  if(length(m0)>1|m0[1]!=id.look$m[1]){stop("different predictor length")}
  
  data   = as.data.frame(matrix(NA,m0,n.var+2,
                                dimnames=list(paste0("1-",1:m0),c("time","status",id.var))))
  pos.col <- 3
  for (var.count in 1:length(var)){
    if (!is.matrix(covar[[var.count]])) {
      data[,pos.col] = covar[[var.count]]
      pos.col <- pos.col+1
    } else {
      for (jj in 1:dim(covar[[var.count]])[2]) {
        data[,pos.col] = covar[[var.count]][,jj]
        pos.col <- pos.col+1
      }
    }
  }
  
  # group
  groupvar <- names(var)[1]
  n.group  = nlevels(data[,groupvar])
  id.group = data.frame(pos = 1:n.group,
                        id = levels(data[,groupvar]),
                        reference = levels(data[,groupvar])==levels(data[,groupvar])[1],
                        active = TRUE,
                        row.names = levels(data[,groupvar]))
  
  # define covariate name corresponding to group
  tmp = labels(terms(model))
  whichw = rep(FALSE,length(tmp))
  for(pw in 1:length(tmp)){
    if(is.factor(data[,tmp[pw]])){
      whichw[pw] = all(!is.na(match(id.group$id,levels(data[,tmp[pw]]))))&
        all(!is.na(match(levels(data[,tmp[pw]]),id.group$id)))
    }
  }
  if(sum(whichw)==1){
    groupvar = tmp[whichw]
  }else{
    if(sum(whichw)==0){
      stop("the variable corresponding to the treatment isn't the expected factor")
    }else{
      stop("2 factors share the same levels")
    }
  }
  # look
  id.look = cbind(id.look,
                  matrix(NA,n.look,n.group,dimnames=list(id.look$id,id.group$id)),
                  matrix(NA,n.look,1,dimnames=list(id.look$id,"t")),
                  matrix(NA,n.look,n.group+1,dimnames=list(id.look$id,c("t(n)",paste0("t(",id.group$id,")")))),     #space for observed times
                  matrix(NA,n.look,n.group+1,dimnames=list(id.look$id,c("ev(n)",paste0("ev(",id.group$id,")")))))     #space for observed events
  
  # generate X matrix for names
  X_tmp <- model.matrix(model[-2], data = data)
  X <- as.matrix(X_tmp[,-1])
  colnames(X) <- colnames(X_tmp)[-1]
  
  if(ncol(X)!=length(hr)){stop("length of 'hr' not compatible with X matrix")}
  names(hr) = colnames(X)
  # targets (only defined once)
  n.target  = length(which)
  if(length(alternative)==1){alternative=rep(alternative,n.target)
  }else{if(length(alternative)!=n.target){stop("length(alternative)!=n.target")}}
  id.target = data.frame(pos = NA,
                         id  = colnames(X)[which],
                         alternative = alternative,
                         group = NA, active = TRUE,
                         look = NA, efficacy = NA, futility = NA,
                         low = NA, mid = NA, high = NA,
                         row.names =  colnames(X)[which])
  id.target$group = sapply(id.target$id,function(x){
    levels(data[,groupvar])[which(sapply(split(X[,x]!=0,data[,groupvar]),any))]
    
  })
  id.target = id.target[order(id.target$group),]
  id.target$pos = 1:n.target
  
  # # delta vector(s)
  # mw = c(match("delta",names(eff.arm.control)),
  #        match("delta",names(fut.arm.control)))
  # none
  if(identical(delta.fut,delta.eff)){
    twodelta = FALSE
    if (length(delta.eff)==1) delta.eff = delta.fut = rep(delta.eff,n.look)
    if (length(delta.eff)!=n.look) stop("length of delta not equal to number of looks")
  }else{
    # both
    if(!(is.null(eff.arm) || is.null(fut.arm))){
      twodelta = TRUE
      # efficacy
      if(length(delta.eff)==1){
        delta.eff = rep(delta.eff,n.look)
      }else{if(length(delta.eff)!=n.look){
        stop("length of delta not equal to number of looks")
      }}
      # futility
      if(length(delta.fut)==1){
        delta.fut = rep(delta.fut,n.look)
      }else{if(length(delta.fut)!=n.look){
        stop("length of delta not equal to number of looks")
      }}
      # one
    }else{
      twodelta = FALSE
      # unique
      #tmp = par[[mw[!is.na(mw)]]]
      tmp <- if (is.null(fut.arm)) delta.eff else delta.fut
      if(length(tmp)==1){
        tmp = rep(tmp,n.look)
      }else{if(length(tmp)!=n.look){
        stop("length of delta not equal to number of looks")
      }}
      # assign
      delta.eff = delta.fut = tmp
    }
  }
  
  if (length(delta.RAR==1)) delta.RAR = rep(delta.RAR,n.look)
  if (length(delta.RAR)!=n.look) stop("length of delta.RAR not equal to number of looks")
  
  # trial stopping rules
  if(is.null(eff.trial) && !is.null(eff.arm)){
    eff.trial = function(eff.target){all(eff.target)}             #use function directly
    #---
  }
  if(is.null(fut.trial) && !is.null(fut.arm)){
    fut.trial = function(fut.target){all(fut.target)}             #use function directly
    #---
  }
  
  # seeds
  if(length(R)==1){id.seed=1:R}else{id.seed=R}
  
  #############################
  # H1
  #############################
  
  if(!all(hr[which]==1)){
    
    message("    Evaluation of H1")   
    H1 = TRUE
    # sequential
    if(computation!="parallel"){
      trial_r = lapply(id.seed,batss.surv.trial,
                       data=data,model=model,family=family,hr=hr,
                       RAR=RAR,RAR.control=RAR.control,twodelta=twodelta,delta.eff=delta.eff,delta.fut=delta.fut,delta.RAR=delta.RAR,
                       eff.arm=eff.arm,eff.trial=eff.trial,
                       eff.arm.control=eff.arm.control,eff.trial.control=eff.trial.control,
                       fut.arm=fut.arm,fut.trial=fut.trial,
                       fut.arm.control=fut.arm.control,fut.trial.control=fut.trial.control,
                       id.target=id.target,n.target=n.target,
                       id.look=id.look,n.look=n.look,prob0=prob0,n.max=n.max,
                       id.group=id.group,n.group=n.group,groupvar=groupvar,
                       surv=surv,surv.control=surv.control,
                       cens=cens,cens.control=cens.control,
                       accr=accr,accr.control=accr.control,accr.type=accr.type,
                       fup=fup,interim=interim,
                       var=var,var.control=var.control,id.var=id.var,n.var=n.var,
                       #linux.os=linux.os,
                       extended=extended,...)
    # parallel
    }else{if(computation=="parallel"){
      # unix via forking
      if(Sys.info()[[1]]!="Windows"){
        trial_r = parallel::mclapply(id.seed,batss.surv.trial,
                                     data=data,model=model,family=family,hr=hr,
                                     RAR=RAR,RAR.control=RAR.control,twodelta=twodelta,delta.eff=delta.eff,delta.fut=delta.fut,delta.RAR=delta.RAR,
                                     eff.arm=eff.arm,eff.trial=eff.trial,
                                     eff.arm.control=eff.arm.control,eff.trial.control=eff.trial.control,
                                     fut.arm=fut.arm,fut.trial=fut.trial,
                                     fut.arm.control=fut.arm.control,fut.trial.control=fut.trial.control,
                                     id.target=id.target,n.target=n.target,
                                     id.look=id.look,n.look=n.look,prob0=prob0,n.max=n.max,
                                     id.group=id.group,n.group=n.group,groupvar=groupvar,
                                     surv=surv,surv.control=surv.control,
                                     cens=cens,cens.control=cens.control,
                                     accr=accr,accr.control=accr.control,accr.type=accr.type,
                                     fup=fup,interim=interim,
                                     var=var,var.control=var.control,id.var=id.var,n.var=n.var,
                                     #linux.os=linux.os,
                                     extended=extended,
                                     mc.cores=mc.cores,mc.set.seed = FALSE,...)
        # windows without forking
      }else{
        cl = parallel::makeCluster(getOption("cl.cores", mc.cores))
        parallel::clusterEvalQ(cl, c(library(INLA),library(foreach)))
        #parallel::clusterExport(cl, transfer, envir = .GlobalEnv)
        #parallel::clusterExport(cl, varlist=c(".expit"), envir = environment())
        trial_r = parallel::parLapply(cl=cl,id.seed,batss.surv.trial,
                                      data=data,model=model,family=family,hr=hr,
                                      RAR=RAR,RAR.control=RAR.control,twodelta=twodelta,delta.eff=delta.eff,delta.fut=delta.fut,delta.RAR=delta.RAR,
                                      eff.arm=eff.arm,eff.trial=eff.trial,
                                      eff.arm.control=eff.arm.control,eff.trial.control=eff.trial.control,
                                      fut.arm=fut.arm,fut.trial=fut.trial,
                                      fut.arm.control=fut.arm.control,fut.trial.control=fut.trial.control,
                                      id.target=id.target,n.target=n.target,
                                      id.look=id.look,n.look=n.look,prob0=prob0,n.max=n.max,
                                      id.group=id.group,n.group=n.group,groupvar=groupvar,
                                      surv=surv,surv.control=surv.control,
                                      cens=cens,cens.control=cens.control,
                                      accr=accr,accr.control=accr.control,accr.type=accr.type,
                                      fup=fup,interim=interim,
                                      var=var,var.control=var.control,id.var=id.var,n.var=n.var,
                                      #linux.os=linux.os,
                                      extended=extended,...)
        stopCluster(cl)
      }
    }}
    
    
    ##
    ## results
    ##

    estimate = batss.res.e(trial_r,id.target)
    tar.p    = batss.res.tp(estimate,id.target)
    tar.g    = batss.res.tg(estimate,id.target)
    eff.p    = batss.res.ep(estimate,id.target,n.look)
    eff.g    = batss.res.eg(estimate,id.target,n.look)
    fut.p    = batss.res.fp(estimate,id.target,n.look)
    fut.g    = batss.res.fg(estimate,id.target,n.look)
    sample   = batss.surv.res.s1(trial_r,group=id.group$id,
                           type=c(apply(estimate[,"type",,drop=FALSE],2:3,paste0,collapse="")),
                           early=sapply(trial_r,function(x){all(x$target$look<dim(x$look)[1])}))
    scenario = batss.res.s2(sample,target=id.target$id)
    res_H1   = list(estimate = estimate,
                    target   = list(par=tar.p,global=tar.g),
                    efficacy = list(par=eff.p,global=eff.g),
                    futility = list(par=fut.p,global=fut.g),
                    sample   = sample,
                    scenario = scenario
                    )
    trial_H1 = trial_r
  }else{
    H1 = FALSE
  }
  
  #############################
  # H0
  #############################
  
  if(H0==TRUE | all(hr[which]==1)){
    
    message("    Evaluation of H0")
    H0 = TRUE
    hr0 = hr
    hr0[which] = 1
    # sequential
    if(computation!="parallel"){
      trial_r = lapply(id.seed,batss.surv.trial,
                       data=data,model=model,family=family,hr=hr0,
                       RAR=RAR,RAR.control=RAR.control,twodelta=twodelta,delta.eff=delta.eff,delta.fut=delta.fut,delta.RAR=delta.RAR,
                       eff.arm=eff.arm,eff.trial=eff.trial,
                       eff.arm.control=eff.arm.control,eff.trial.control=eff.trial.control,
                       fut.arm=fut.arm,fut.trial=fut.trial,
                       fut.arm.control=fut.arm.control,fut.trial.control=fut.trial.control,
                       id.target=id.target,n.target=n.target,
                       id.look=id.look,n.look=n.look,prob0=prob0,n.max=n.max,
                       id.group=id.group,n.group=n.group,groupvar=groupvar,
                       surv=surv,surv.control=surv.control,
                       cens=cens,cens.control=cens.control,
                       accr=accr,accr.control=accr.control,accr.type=accr.type,
                       fup=fup,interim=interim,
                       var=var,var.control=var.control,id.var=id.var,n.var=n.var,
                       #linux.os=linux.os,
                       extended=extended,...)
    # parallel 
    }else{if(computation=="parallel"){
      # unix via forking
      if(Sys.info()[[1]]!="Windows"){
        trial_r = parallel::mclapply(id.seed,batss.surv.trial,
                                     data=data,model=model,family=family,hr=hr0,
                                     RAR=RAR,RAR.control=RAR.control,twodelta=twodelta,delta.eff=delta.eff,delta.fut=delta.fut,delta.RAR=delta.RAR,
                                     eff.arm=eff.arm,eff.trial=eff.trial,
                                     eff.arm.control=eff.arm.control,eff.trial.control=eff.trial.control,
                                     fut.arm=fut.arm,fut.trial=fut.trial,
                                     fut.arm.control=fut.arm.control,fut.trial.control=fut.trial.control,
                                     id.target=id.target,n.target=n.target,
                                     id.look=id.look,n.look=n.look,prob0=prob0,n.max=n.max,
                                     id.group=id.group,n.group=n.group,groupvar=groupvar,
                                     surv=surv,surv.control=surv.control,
                                     cens=cens,cens.control=cens.control,
                                     accr=accr,accr.control=accr.control,accr.type=accr.type,
                                     fup=fup,interim=interim,
                                     var=var,var.control=var.control,id.var=id.var,n.var=n.var,
                                     #linux.os=linux.os,
                                     extended=extended,
                                     mc.cores=mc.cores,mc.set.seed = FALSE,...)
        # windows without forking
      }else{
        cl = parallel::makeCluster(getOption("cl.cores", mc.cores))
        parallel::clusterEvalQ(cl, c(library(INLA),library(foreach)))
        #parallel::clusterExport(cl, transfer, envir = .GlobalEnv)
        #parallel::clusterExport(cl, c(".expit"), envir = environment())
        trial_r = parallel::parLapply(cl=cl,id.seed,batss.surv.trial,
                                      data=data,model=model,family=family,hr=hr0,
                                      RAR=RAR,RAR.control=RAR.control,twodelta=twodelta,delta.eff=delta.eff,delta.fut=delta.fut,delta.RAR=delta.RAR,
                                      eff.arm=eff.arm,eff.trial=eff.trial,
                                      eff.arm.control=eff.arm.control,eff.trial.control=eff.trial.control,
                                      fut.arm=fut.arm,fut.trial=fut.trial,
                                      fut.arm.control=fut.arm.control,fut.trial.control=fut.trial.control,
                                      id.target=id.target,n.target=n.target,
                                      id.look=id.look,n.look=n.look,prob0=prob0,n.max=n.max,
                                      id.group=id.group,n.group=n.group,groupvar=groupvar,
                                      surv=surv,surv.control=surv.control,
                                      cens=cens,cens.control=cens.control,
                                      accr=accr,accr.control=accr.control,accr.type=accr.type,
                                      fup=fup,interim=interim,
                                      var=var,var.control=var.control,id.var=id.var,n.var=n.var,
                                      #linux.os=linux.os,
                                      extended=extended,...)
        stopCluster(cl)
      }

    }}
    
    ##
    ## results
    ##
    estimate = batss.res.e(trial_r,id.target)
    tar.p    = batss.res.tp(estimate,id.target)
    tar.g    = batss.res.tg(estimate,id.target)
    eff.p    = batss.res.ep(estimate,id.target,n.look)
    eff.g    = batss.res.eg(estimate,id.target,n.look)
    fut.p    = batss.res.fp(estimate,id.target,n.look)
    fut.g    = batss.res.fg(estimate,id.target,n.look)
    sample   = batss.surv.res.s1(trial_r,group=id.group$id,
                                 type=c(apply(estimate[,"type",,drop=FALSE],2:3,paste0,collapse="")),
                                 early=sapply(trial_r,function(x){all(x$target$look<dim(x$look)[1])}))
    scenario = batss.res.s2(sample,target=id.target$id)
    res_H0   = list(estimate = estimate,
                    target   = list(par=tar.p,global=tar.g),
                    efficacy = list(par=eff.p,global=eff.g),
                    futility = list(par=fut.p,global=fut.g),
                    sample   = sample,
                    scenario = scenario)
    trial_H0 = trial_r
    
  }else{
    H0 = FALSE
  }
  
  ##
  ## output
  ##
  
  message("    Results") 
  look        = id.look[,c("pos","id","n","m")]
  FE  = data.frame(pos=1:ncol(X),id=colnames(X),target=FALSE,
                   row.names = colnames(X))
  FE[id.target$id,"target"] = TRUE
  if(H0){FE[,'HR (H0)'] = hr0}
  if(H1){FE[,'HR (H1)'] = hr}
  #
  
  par = list(model=model,family=family,var=var,
             eff.arm=eff.arm,eff.trial=eff.trial,
             eff.arm.control=eff.arm.control,eff.trial.control=eff.trial.control,
             fut.arm=fut.arm,fut.trial=fut.trial,
             fut.arm.control=fut.arm.control,fut.trial.control=fut.trial.control,
             RAR=RAR, RAR.control=RAR.control,
             surv=surv,surv.control=surv.control,
             cens=cens,cens.control=cens.control,
             accr=accr,accr.control=accr.control,accr.type=accr.type,
             fup=fup,interim=interim,
             seed=id.seed,H0=H0,H1=H1,prob0=prob0, 
             group=id.group[,c("pos","id","reference")], version=utils::packageVersion("BATSS"))
  out = list(hr = FE, look = look, par = par)
 
  if(H0){
    out$H0 = res_H0
    if(extended>0){out$H0$trial = trial_H0}
  }
  if(H1){
    out$H1 = res_H1
    if(extended>0){out$H1$trial = trial_H1}
  }
 
  #---
  out$call <- call
  out$type <- "surv"
  #---
  class(out) = "batss"
  out
}






