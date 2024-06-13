#' @name bats.glm
#' @title BATS for generalised linear models
#' @description Simulation of Bayesian adaptive trials with GLM endpoint using INtegrated Laplace Approximation (INLA).  
#' @param model A character string indicating a symbolic description of the model to be fitted (as in the lm() and glm() functions).
#' @param var A list. Each entry corresponds to a variable described under 'model' and indicates, by means of a character string, how each variable is generated. The list names have to match the variable names as indicated under 'model' and the first list element corresponds to the outcome. The variable corresponding to the 'treatment' has to be a factor with levels corresponding to the names indicated for prob0 (see below).
#' @param var.control An optional list of control parameters for functions in 'var'. The names of the list items need to correspond with the names of 'var'. Each element is another list with names of the elements corresponding to the parameter names of the functions specified in 'var'. 
#' @param family A character string indicating the name of the conditional distribution as described the package INLA. Default set to 'gaussian'.
#' @param link A character string describing the link function to be used in the model to relate the outcome to the set of predictors: "identity", "logit" and "log" are the currently available options. Default set to 'identity'.
#' @param beta A numerical vector of parameter values for the linear predictor. Its length has to match the number of column of the X matrix induced by the formula indicated under 'model' (check ?model.matrix).
#' @param which A numerical vector indicating the position of the target parameters.
#' @param R a scalar corresponding to the number of Monte Carlo trials or a vector of natural numbers to be used as seeds (?set.seed) for the different Monte Carlo trial (the vector length will thus correspond to the number of Monte Carlo trials).
#' @param alternative The direction of the alternative hypothesis corresponding to each target parameter indicated under 'which'. Possibilities are 'greater' (default) or 'less'. If the vector is of length 1, the same direction will be used for all target parameter tests.
#' @param RAR A character string indicating the how the response-adaptive randomised probabilities corresponding to each group, reference group included, are defined. If RAR=NULL, equal probabilities of being attributed to each (active) group are used. Check ?BATS::prob.fun for examples.
#' @param RAR.control An optional list of control parameters the function in 'RAR'. 
#' @param N A scalar indicating the total sample size.
#' @param interim A list of parameters related to interim analyses. Currently, only the parameter 'recruited' is available.  It consists in a vector of integers indicating the number of completed observations at each look, last excluded, in increasing order.
#' @param prob0 A named vector with initial allocation probabilities. Names need to correspond with the levels of the grouping variable.
#' @param delta.eff A scalar or vector defining the boundaries for efficacy.
#' @param delta.fut A scalar or vector defining the boundaries for futility.
#' @param delta.RAR A scalar or vector defining the boundaries for RAR.
#' @param eff.arm A character string showing an expression indicating how to evaluate if efficacy has been achieved at a given look given the information (ingredients) available at that stage for a given target parameter. The output of this expression must be logical for each target parameter. Check ?BATS::eff.arm_example for examples. 
#' @param eff.arm.control A list of parameters for functions in 'eff.arm'
#' @param eff.trial A character string showing an expression (a function of ingredient eff.target) indicating if the trial can be stopped for efficacy given the output of the function specified in eff.arm on all parameters. The output of this expression must be a logical of length one. When eff.trial=NULL (default), the trial stops for efficacy when all target parameters are found to be effective. Check ?BATS::eff.trial_example for examples. 
#' @param eff.trial.control A list of parameters for functions in 'eff.trial'
#' @param fut.arm A character string showing an expression indicating how to evaluate if futility has been achieved at a given look given the information (ingredients) available at that stage for a given target parameter. The output of this expression must be a logical for each target parameter. Check ?BATS::fut.arm_example for examples. 
#' @param fut.arm.control A list of parameters for functions in 'fut.arm'
#' @param fut.trial A character string showing an expression (a function of ingredient fut.target) indicating if the trial can be stopped for futility given the output of the function specified in fut.arm on all parameters. The output of this expression must be logical of length one. When fut.trial=NULL (default), the trial stops when all target parameters are found to be futile. Check ?BATS::fut.trial_example for examples. 
#' @param fut.trial.control A list of parameters for functions in 'fut.trial'
#' @param H0 A logical indicating whether the simulation should also consider the case with all target parameters set to 0 (or to the first value of par$delta) to check the probability of rejecting the hypothesis that the target parameter value is equal to par$delta individually (type I error when par$delta=0) or globally (FWER when par$delta=0).
#' @param computation A character string indicating how to parallelise the task, with possibilities 'parallel' when running this task on a single computer (the option 'cluster' - using the Slurm workload manager - is temporarily removed due to issues with cluster variants of INLA shared objects not being available).
#' @param mc.cores An integer indicating how many CPUs to use when computation='parallel' (Default to 3 if no global option was set).
#' @param extended an integer indicating the type of results to be returned. 0 (default) provides summary statistics, 1 adds the results of each Monte Carlo trial and 2 additionally returns the each Monte Carlo datasets. BATS::bats.combine requires extended > 0 as the function needs to merge results of different sets of seeds.
#' @param ... Additional arguments to control fitting in INLA
#' @returns The function [bats.glm] returns an S3 object of class 'bats.glm' with available print/summary/plot functions.
#' @export
#' @examples
#'\dontrun{
#' # Example: 
#' # fixed efficacy (0.975) and futility (0.05) stops, no RAR, 
#' # 3 groups with group means C = 1 (ref), T1 = 2, T2 = 3, 
#' # Gaussian conditional distribution with sigma = 2, 6 looks.
#' 
#' sim = bats.glm(  
#'                model            = y ~ group,   
#'                var              = list(y     = rnorm,
#'                                        group = alloc.balanced),
#'                var.control      = list(y = list(sd = 5)),
#'                beta             = c(1, 1, 2),
#'                which            = c(2:3),
#'                R                = 25,
#'                alternative      = "greater",
#'                N                = 200,
#'                interim          = list(recruited = seq(100, 180, 20)),
#'                prob0            = c(C = 1/3, T1 = 1/3, T2 = 1/3),
#'                eff.arm          = eff.arm.simple,
#'                eff.arm.control  = list(b = 0.975),
#'                fut.arm          = fut.arm.simple,
#'                fut.arm.control  = list(b = 0.05),
#'                computation      = "parallel",
#'                H0               = TRUE,
#'                mc.cores         = parallel::detectCores()-1)
#' }
#' @export
bats.glm = function(
model,var,var.control=NULL,family="gaussian",link="identity",
beta,which,alternative = "greater",R=1e+4,N,interim,prob0,
delta.eff=0,delta.fut=delta.eff, delta.RAR=0,
eff.arm,eff.arm.control=NULL,
eff.trial=NULL,eff.trial.control=NULL,
fut.arm,fut.arm.control=NULL,
fut.trial=NULL,fut.trial.control=NULL,
RAR=NULL,RAR.control=NULL,
H0=TRUE,computation="parallel",
mc.cores=getOption("mc.cores", 3L),
#linux.os = NA,
extended=0, ...){

#---    
call <- match.call()                       # save call
model <- as.formula(model)                 # allow for string and formula input
#--- 
                
##
## dataset structure and useful definitions 
##

cat("\n\tInitialisation\n")    

#some checks
n = m = prob = NULL
#error messages
if (!is.character(model) && !is.formula(model)) 
  stop("invalid 'model' argument")
if (!is.list(var) || !all(sapply(var,is.function)))
  stop("'var' must be a list of functions")
if (is.null(intersect(names(var),intersect(setdiff(unlist(strsplit(all.vars(model),"1")),all.vars(model)),
                                           setdiff(unlist(strsplit(all.vars(model),"2")),all.vars(model))))) && 
    !setequal(all.vars(model),names(var)))
  stop("all variables in the model formula must have a generating function in the 'var' list")
if (!(family %in% names(inla.models()$likelihood))){
  stop("invalid 'family' argument, see help files and inla documentation for available families")
} else {
  if (!(family %in% c("gaussian","binomial","nbinomial","poisson")))
    warning("functionality only tested for gaussian, binomial, negative binomial and poisson distributions")
}
if (!(link %in% inla.models()$likelihood[[family]]$link) || !(link %in% c("identity","log","logit","probit","robit","cauchit","loglog","cloglog")))
  stop("'link' not supported, see help files and inla documentation for available link functions")
if (!is.null(interim)){
  if(!inherits(interim,"list")){stop("'interim' should be a list")}
  interim.recruited = interim$recruited
}else{
    stop("'interim' should be provided")
  } 
if (!is.null(interim.recruited) && !is.numeric(unlist(interim.recruited))) 
  stop("'interim.recruited' must be a (list of) numeric vector(s)")
if (!is.null(interim.recruited) && any(interim.recruited < 0)) 
  stop("negative interim recruitment numbers not allowed")
if ((N < 0) || length(N)>1) {
  stop("total sample size 'N' must be a positive scalar")
  N <-  floor(N)
}
if (length(which)>length(beta))
  stop("number of targets greater than number of parameters")
if ((!is.null(RAR) && !(is.function(RAR))) || 
    (!is.null(eff.arm) && !(is.function(eff.arm))) || (!is.null(eff.trial) && !(is.function(eff.trial))) ||
    (!is.null(fut.arm) && !(is.function(fut.arm))) || (!is.null(fut.trial) && !(is.function(fut.trial))))
  stop("'RAR', 'eff.arm', 'eff.trial', 'fut.arm', 'fut.trial' must be functions or NULL")

#warnings
if (!is.null(interim.recruited) && any(interim.recruited > N)) {
  warning("some interim analyses are outside the maximum sample size and will be ignored")
  interim.recruited <- interim.recruited[interim.recruited<N]
}
if (all(prob0<0) && sum(prob0)!=1)
  warning("sum of 'prob0' not equal to 1")



# size per look
if(sum(!is.na(match(c("m0","m"),names(interim.recruited))))==2){
   size_look = seq(interim.recruited$m0,N,interim.recruited$m)  
   size_look[length(size_look)] = N
}else{
  interim.recruited <- sort(interim.recruited)        #sort interim recruited to make sure the values are ordered
   size_look = c(interim.recruited,N)  
   }
n.look  = length(size_look)   
id.look = data.frame(pos = 1:n.look,
                     id  = paste0("n=",size_look),
                     n   = size_look,
                     m   = c(size_look[1],size_look[-1]-size_look[-n.look]))
# generate predictors
env0 = new.env()
assign("m", id.look$m[1], envir = env0)
assign("n", id.look$n[1], envir = env0)
assign("prob",prob0,envir = env0)
    
assign("var.control", var.control, envir = env0)  
covar <- vector("list",length(var[-1]))
for (ii in 1:length(var[-1])) {
  tmp_nam <- names(var)[ii+1]
  args_ <- plyr::.(n=n,m=m,prob=prob)
  if (tmp_nam %in% names(var.control)) args_ <- c(args_, var.control[[tmp_nam]])
  covar[[ii]] <- R.utils::doCall(var[[ii+1]], envir = env0, args = args_)   #call functions directly
  if (is.matrix(covar[[ii]])) colnames(covar[[ii]]) <- paste0(tmp_nam,1:dim(covar[[ii]])[2])
}

n.var <- length(all.vars(model))
id.var <- names(var)
if (any(sapply(covar,is.matrix))) {
  where.mat <- which(sapply(covar,is.matrix))
  tmp_var <- covar[where.mat]
  id.var <- names(var)[1]
  for (ii in 1:length(covar)) {
    if (ii %in% (where.mat)) {
      id.var <- c(id.var,colnames(covar[[ii]]))
    } else {
      id.var <- c(id.var,names(var[ii+1]))
    }
  }
}

m0     = length(covar[[1]])
if(length(m0)>1|m0[1]!=id.look$m[1]){stop("different predictor length")}

data <- as.data.frame(matrix(NA,m0,n.var,
                     dimnames=list(paste0("1-",1:m0),id.var)))
pos.col <- 2
for (var.count in 1:(length(var)-1)){
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
groupvar <- names(var)[2]
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
                matrix(NA,n.look,n.group,
                dimnames=list(id.look$id,id.group$id)))
# generate X matrix for names
X = model.matrix(model[-2], data = data)                               #create model matrix straight from formula

if(ncol(X)!=length(beta)){stop("length of 'beta' not compatible with X matrix")}
names(beta) = colnames(X)
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


# delta vector(s)
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
}}  

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

if(!all(beta[which]==0)){

    cat("\tEvaluation of H1\n")           
    H1 = TRUE
    # lapply
    if(computation!="parallel"&computation!="cluster"){
        trial_r = lapply(id.seed,bats.trial,
               data=data,model=model,link=link,family=family,beta=beta,
               RAR=RAR,RAR.control=RAR.control,twodelta=twodelta,delta.eff=delta.eff,delta.fut=delta.fut,
               eff.arm=eff.arm,eff.trial=eff.trial,delta.RAR=delta.RAR,
               eff.arm.control=eff.arm.control,eff.trial.control=eff.trial.control,
               fut.arm=fut.arm,fut.trial=fut.trial,
               fut.arm.control=fut.arm.control,fut.trial.control=fut.trial.control,
               id.target=id.target,n.target=n.target,
               id.look=id.look,n.look=n.look,prob0=prob0,
               id.group=id.group,n.group=n.group,groupvar=groupvar,
               var=var,var.control=var.control,id.var=id.var,n.var=n.var,
               #linux.os=linux.os,
               extended = extended, ...)
    }else{if(computation=="parallel"){           
        # unix via forking
        if(Sys.info()[[1]]!="Windows"){
        trial_r = parallel::mclapply(id.seed,bats.trial,
                                      data=data,model=model,link=link,family=family,beta=beta,
                                      RAR=RAR,RAR.control=RAR.control,twodelta=twodelta,delta.eff=delta.eff,delta.fut=delta.fut,
                                      eff.arm=eff.arm,eff.trial=eff.trial,delta.RAR=delta.RAR,
                                      eff.arm.control=eff.arm.control,eff.trial.control=eff.trial.control,
                                      fut.arm=fut.arm,fut.trial=fut.trial,
                                      fut.arm.control=fut.arm.control,fut.trial.control=fut.trial.control,
                                      id.target=id.target,n.target=n.target,
                                      id.look=id.look,n.look=n.look,prob0=prob0,
                                      id.group=id.group,n.group=n.group,groupvar=groupvar,
                                      var=var,var.control=var.control,id.var=id.var,n.var=n.var,
                                      # linux.os=linux.os,
                                      extended = extended, mc.cores=mc.cores,mc.set.seed = FALSE,...)
        # windows without forking
        }else{
        cl = parallel::makeCluster(getOption("cl.cores", mc.cores))
        parallel::clusterEvalQ(cl, c(library(INLA)))
        #parallel::clusterExport(cl, transfer, envir = .GlobalEnv)           
#        parallel::clusterExport(cl, c(".expit"), envir = environment())           
        trial_r = parallel::parLapply(cl=cl,id.seed,bats.trial,
                                      data=data,model=model,link=link,family=family,beta=beta,
                                      RAR=RAR,RAR.control=RAR.control,twodelta=twodelta,delta.eff=delta.eff,delta.fut=delta.fut,
                                      eff.arm=eff.arm,eff.trial=eff.trial,delta.RAR=delta.RAR,
                                      eff.arm.control=eff.arm.control,eff.trial.control=eff.trial.control,
                                      fut.arm=fut.arm,fut.trial=fut.trial,
                                      fut.arm.control=fut.arm.control,fut.trial.control=fut.trial.control,
                                      id.target=id.target,n.target=n.target,
                                      id.look=id.look,n.look=n.look,prob0=prob0,
                                      id.group=id.group,n.group=n.group,groupvar=groupvar,
                                      var=var,var.control=var.control,id.var=id.var,n.var=n.var,
                                      #linux.os=linux.os,
                                      extended = extended,...)
        stopCluster(cl)    
        }
        # cluster via rslurm (>=0.6)
        }else{
        stop("'cluster' option temporarily not available")            
    }}     
    
    
    ##
    ## results
    ##
    estimate = bats.res.e(trial_r,id.target)
    tar.p    = bats.res.tp(estimate,id.target)
    tar.g    = bats.res.tg(estimate,id.target)
    eff.p    = bats.res.ep(estimate,id.target,n.look)
    eff.g    = bats.res.eg(estimate,id.target,n.look)
    fut.p    = bats.res.fp(estimate,id.target,n.look)
    fut.g    = bats.res.fg(estimate,id.target,n.look)
    sample   = bats.res.s1(trial_r,group=id.group$id,
                           type=c(apply(estimate[,"type",,drop=FALSE],2:3,paste0,collapse="")),
                           early=apply(estimate[,"look",,drop=TRUE]<n.look,2,all))
    scenario = bats.res.s2(sample,target=id.target$id)
    res_H1   = list(estimate = estimate,
                    target   = list(par=tar.p,global=tar.g),
                    efficacy = list(par=eff.p,global=eff.g),
                    futility = list(par=fut.p,global=fut.g),
                    sample=sample,scenario=scenario)
    trial_H1    = trial_r
}else{
    H1 = FALSE
}

#############################
# H0
#############################

if(H0==TRUE | all(beta[which]==0)){

    cat("\tEvaluation of H0\n")
    H0 = TRUE
    beta0 = beta
    beta0[which] = 0
    # lapply
    if(computation!="parallel"&computation!="cluster"){
        trial_r = lapply(id.seed,bats.trial,
                         data=data,model=model,link=link,family=family,beta=beta0,
                         RAR=RAR,RAR.control=RAR.control,twodelta=twodelta,delta.eff=delta.eff,delta.fut=delta.fut,
                         eff.arm=eff.arm,eff.trial=eff.trial,delta.RAR=delta.RAR,
                         eff.arm.control=eff.arm.control,eff.trial.control=eff.trial.control,
                         fut.arm=fut.arm,fut.trial=fut.trial,
                         fut.arm.control=fut.arm.control,fut.trial.control=fut.trial.control,
                         id.target=id.target,n.target=n.target,
                         id.look=id.look,n.look=n.look,prob0=prob0,
                         id.group=id.group,n.group=n.group,groupvar=groupvar,
                         var=var,var.control=var.control,id.var=id.var,n.var=n.var,
                         #linux.os=linux.os,
                         extended = extended, ...)
    }else{if(computation=="parallel"){           
        # unix via forking
        if(Sys.info()[[1]]!="Windows"){
        trial_r = parallel::mclapply(id.seed,bats.trial,
                                     data=data,model=model,link=link,family=family,beta=beta0,
                                     RAR=RAR,RAR.control=RAR.control,twodelta=twodelta,delta.eff=delta.eff,delta.fut=delta.fut,
                                     eff.arm=eff.arm,eff.trial=eff.trial,delta.RAR=delta.RAR,
                                     eff.arm.control=eff.arm.control,eff.trial.control=eff.trial.control,
                                     fut.arm=fut.arm,fut.trial=fut.trial,
                                     fut.arm.control=fut.arm.control,fut.trial.control=fut.trial.control,
                                     id.target=id.target,n.target=n.target,
                                     id.look=id.look,n.look=n.look,prob0=prob0,
                                     id.group=id.group,n.group=n.group,groupvar=groupvar,
                                     var=var,var.control=var.control,id.var=id.var,n.var=n.var,
                                     #linux.os=linux.os,
                                     extended = extended, 
                                     mc.cores=mc.cores,mc.set.seed = FALSE,...)
        # windows without forking
        }else{
        cl = parallel::makeCluster(getOption("cl.cores", mc.cores))
        parallel::clusterEvalQ(cl, c(library(INLA)))
        #parallel::clusterExport(cl, transfer, envir = .GlobalEnv)           
#        parallel::clusterExport(cl, c(".expit"), envir = environment())           
        trial_r = parallel::parLapply(cl=cl,id.seed,bats.trial,
                                      data=data,model=model,link=link,family=family,beta=beta0,
                                      RAR=RAR,RAR.control=RAR.control,twodelta=twodelta,delta.eff=delta.eff,delta.fut=delta.fut,
                                      eff.arm=eff.arm,eff.trial=eff.trial,delta.RAR=delta.RAR,
                                      eff.arm.control=eff.arm.control,eff.trial.control=eff.trial.control,
                                      fut.arm=fut.arm,fut.trial=fut.trial,
                                      fut.arm.control=fut.arm.control,fut.trial.control=fut.trial.control,
                                      id.target=id.target,n.target=n.target,
                                      id.look=id.look,n.look=n.look,prob0=prob0,
                                      id.group=id.group,n.group=n.group,groupvar=groupvar,
                                      var=var,var.control=var.control,id.var=id.var,n.var=n.var,
                                      #linux.os=linux.os,
                                      extended=extended, ...)
        stopCluster(cl)    
        }
        # cluster via rslurm (>=0.6)
        }else{
        stop("'cluster' option temporarily not available")            
    }}     
    
    ##
    ## results
    ##
    #browser()
    estimate = bats.res.e(trial_r,id.target)
    tar.p    = bats.res.tp(estimate,id.target)
    tar.g    = bats.res.tg(estimate,id.target)
    eff.p    = bats.res.ep(estimate,id.target,n.look)
    eff.g    = bats.res.eg(estimate,id.target,n.look)
    fut.p    = bats.res.fp(estimate,id.target,n.look)
    fut.g    = bats.res.fg(estimate,id.target,n.look)
    sample   = bats.res.s1(trial_r,group=id.group$id,
                           type=c(apply(estimate[,"type",,drop=FALSE],2:3,paste0,collapse="")),
                           early=apply(estimate[,"look",,drop=TRUE]<n.look,2,all))
    scenario = bats.res.s2(sample,target=id.target$id)
    res_H0   = list(estimate = estimate,
                    target   = list(par=tar.p,global=tar.g),
                    efficacy = list(par=eff.p,global=eff.g),
                    futility = list(par=fut.p,global=fut.g),
                    sample=sample,scenario=scenario)
    trial_H0    = trial_r 
     
}else{
    H0 = FALSE
}

##
## output
##

cat("\tResults 1\n")    
look        = id.look[,c("pos","id","n","m")]
FE  = data.frame(pos=1:ncol(X),id=colnames(X),target=FALSE,
                 row.names = colnames(X))
FE[id.target$id,"target"] = TRUE
if(H0){FE[,'Beta (H0)'] = beta0}
if(H1){FE[,'Beta (H1)'] = beta}
# 
# browser("browser")
cat("\tResults 2\n")    
par = list(RAR=RAR,group=id.group[,c("pos","id","reference")],
           seed=id.seed,H0=H0,H1=H1)
out = list(beta = FE, look = look, par=par)    
cat("\tResults 3\n")        
if(H0){
    out$H0 = res_H0
    if(extended>0){out$H0$trial = trial_H0}
    }
if(H1){
    out$H1 = res_H1
    if(extended>0){out$H1$trial = trial_H1}
    }
cat("\tResults 4\n") 
#---
out$call <- call
out$type <- "glm"
#---
class(out) = "bats"
out
}

