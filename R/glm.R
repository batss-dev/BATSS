#' @name batss.glm
#' @title BATSS for generalised linear models
#' @description Simulation of Bayesian adaptive trials with GLM endpoint using INtegrated Laplace Approximation (INLA).  
#' @param model an object of class '\link[stats]{formula}' indicating a symbolic description of the model to be fitted (as in the \link[stats]{lm} and \link[stats]{glm} functions).
#' @param var A list. Each entry corresponds to a variable described under '`model`' and indicates the name of a function allowing to generate variates (like \link[stats]{rnorm} and \link[stats]{rexp}, for example). The list names have to match the variable names unded in '`model`' and its first element should correspond to the model outcome. The grouping variable corresponding to the target parameters has to be of class '\link[base]{factor}' with levels corresponding to the names indicated in argument `prob0` (see below).
#' @param var.control An optional list of control parameters for the functions indicated in '`var`'. The names of the list items need to correspond to the names used in '`var`'. Each element is another list with names of the elements corresponding to the parameter names of the functions specified in '`var`'. 
#' @param family A character string indicating the name of the conditional distribution as described in the package INLA (check \link[INLA]{inla.list.models}). Default set to '`gaussian`'.
#' @param link A character string describing the link function to be used in the model to relate the outcome to the set of predictors: 'identity', 'log', 'logit', 'probit', 'robit', 'cauchit', 'loglog' and 'cloglog' are the currently available options. Default set to 'identity'.
#' @param beta A numerical vector of parameter values for the linear predictor. Its length has to match the number of column of the **X** matrix induced by the formula indicated under '`model`' (check \link[stats]{model.matrix}).
#' @param which A numerical vector indicating the position of the target `beta` parameters.
#' @param R a vector of natural numbers to be used as seeds (check \link[base]{set.seed}) for the different Monte Carlo trials (the vector length will thus correspond to the number of Monte Carlo trials). When `R` is a scalar, seeds `1` to `R` are used, where `R` corresponds to the number of Monte Carlo trials. 
#' @param alternative A vector of strings providing the one-sided direction of the alternative hypothesis corresponding to each target parameter indicated under '`which`' (in the same order). Possibilities are 'greater' (default) or 'less'. If the vector is of length 1, the same direction will be used for all target parameter tests.
#' @param RAR A function defining the response-adaptive randomisation probabilities of each group - reference group included - with the same group names and ordering as used in '`prob0`'. Arguments of this function will typically consider BATSS 'ingredients'. Check [RAR.trippa] and [RAR.optimal] for examples. If `RAR = NULL` (default), the probabilities/ratios indicated under `prob0` will be used throughout (fixed allocation probabilities).
#' @param RAR.control An optional list of control parameters for the function provided in '`RAR`'. 
#' @param N A scalar indicating the maximum sample size.
#' @param interim A list of parameters related to interim analyses. Currently, only '`recruited`' is available.  It consists in a vector of integers indicating the number of completed observations at each look, last excluded, in increasing order.
#' @param prob0 A named vector with initial allocation probabilities. Names need to correspond to the levels of the grouping variable. If `RAR = NULL`, these probabilities/ratios will be used throughout (fixed allocation probabilities).
#' @param delta.eff A vector (of length equal to the number of looks (i.e., number of interims + 1)) of clinically meaningful treatment effect values (on the linear predictor scale) to be used to define the efficacy-related posterior probabilities for each target parameter at each look. If a scalar is provided, the same value is used at each look. The default is `delta.eff = 0`. 
#' @param delta.fut A vector (of length equal to the number of looks (i.e., number of interims + 1)) of clinically meaningful treatment effect values (on the linear predictor scale) to be used to define the futility-related posterior probabilities for each target parameter at each look. If a scalar is provided, the same value is used at each look. The default is `delta.fut = delta.eff`. 
#' @param delta.RAR A vector (of length equal to the number of looks (i.e., number of interims + 1)) of clinically meaningful treatment effect values (on the linear predictor scale) to be used to define the RAR-related posterior probabilities for each target parameter at each look. If a scalar is provided, the same value is used at each interim analysis. The default is `delta.RAR = 0`. Note that, when a vector is provided, its last value is ignored as no randomisation is made at the last look.
#' @param eff.arm A function defining if efficacy has been achieved at a given look given the information available at that stage a given target parameter. The output of this function must be a \link[base]{logical} (of length 1). Arguments of this function will typically consider BATSS 'ingredients'. Check [eff.arm.simple] and [eff.arm.infofract] for examples. 
#' @param eff.arm.control An optional list of parameters for the function indicated in '`eff.arm`'.
#' @param eff.trial A function defining if the trial can be stopped for efficacy given the output of the function indicated in '`eff.arm`'. The output of this function must be a \link[base]{logical} of length one. Arguments of this function will typically only consider the BATSS ingredient `eff.target`. Check [eff.trial.all] and [eff.trial.any] for examples. When `eff.trial = NULL` (default), the trial stops for efficacy when *all* target parameters are found to be effective (like in [eff.trial.all]). 
#' @param eff.trial.control An optional list of parameters for the function indicated in '`eff.trial`'.
#' @param fut.arm A function defining if futility has been achieved at a given look given the information available at that stage for each target parameter. The output of this function must be a \link[base]{logical} (of length 1). Arguments of this function will typically consider BATSS 'ingredients'. Check [fut.arm.simple] to see an example of such a function. 
#' @param fut.arm.control An optional list of parameters for the function indicated in '`fut.arm`'.
#' @param fut.trial A function defining if the trial can be stopped for futility given the output of the function indicated in '`fut.arm`'. The output of this function must be a \link[base]{logical} of length one. Arguments of this function will typically only consider the BATSS ingredient `fut.target`. Check [fut.trial.all] for an example of such a function. When `fut.trial = NULL` (default), the trial stops for futility when *all* target parameters are found to be futile (like in [fut.trial.all]).
#' @param fut.trial.control An optional list of parameters for the function indicated in '`fut.trial`'.
#' @param H0 A logical indicating whether the simulation should also consider the case with all target parameters set to 0 to check the probability of rejecting the hypothesis that the target parameter value is equal to 0 individually (pairwise type I error) or globally (family-wise error rate). Default set to `H0=TRUE`.
#' @param computation A character string indicating how the computation should be performed. Possibilities are 'parallel' or 'sequential' with default `computation="parallel"` meaning that the computation is split between `mc.cores`. 
#' @param mc.cores An integer indicating the number of CPUs to be used when `computation="parallel"` (Default to 3 if no global '`mc.cores`' global option is available via \link[base]{getOption}).
#' @param extended an integer indicating the type of results to be returned. 0 (default) provides summary statistics, 1 adds the results of each Monte Carlo trial and 2 additionally returns each Monte Carlo dataset. [batss.combine] requires extended > 0 as the function needs to merge results of different sets of seeds.
#' @param ... Additional arguments to control fitting in \link[INLA]{inla}.
#' @returns The function [batss.glm] returns an S3 object of class 'batss' with available print/summary/plot functions.
#' @export
#' @seealso [summary.batss()] and [plot.batss()] for detailed summaries and plots, and [batss.combine()] to combine different evaluations of [batss.glm] considering the same trial design but different sets of seeds (useful for cluster computation). 
#' @examples
#'\dontrun{
#' # Example: 
#' # * Gaussian conditional distribution with sigma = 5
#' # * 3 groups with group means 'C' = 1 (ref), 'T1' = 2, 'T2' = 3,
#' #     where higher means correspond to better outcomes 
#' # * 5 interim analyses occurring when n = 100, 120, 140, 160, and 180
#' # * fixed and equal allocation probabilities per arm (i.e., no RAR)
#' # * max sample size = 200 
#' # * efficacy stop per arm when the prob of the corresponding parameter 
#' #     being greater than 0 is greater than 0.975 (?eff.arm.simple)
#' # * futility stop per arm when the prob of the corresponding parameter 
#' #     being greater than 0 is smaller than 0.05 (?fut.arm.simple) 
#' # * trial stop once all arms have stopped (?eff.trial.all and ?fut.trial.all)
#' #     or the max sample size was reached 
#' 
#' sim = batss.glm(model            = y ~ group,   
#'                 var              = list(y     = rnorm,
#'                                         group = alloc.balanced),
#'                 var.control      = list(y = list(sd = 5)),
#'                 beta             = c(1, 1, 2),
#'                 which            = c(2:3),
#'                 alternative      = "greater",
#'                 R                = 25,
#'                 N                = 200,
#'                 interim          = list(recruited = seq(100, 180, 20)),
#'                 prob0            = c(C = 1/3, T1 = 1/3, T2 = 1/3),
#'                 eff.arm          = eff.arm.simple,
#'                 eff.arm.control  = list(b = 0.975),
#'                 fut.arm          = fut.arm.simple,
#'                 fut.arm.control  = list(b = 0.05),
#'                 computation      = "parallel",
#'                 H0               = TRUE,
#'                 mc.cores         = parallel::detectCores()-1)
#' }
#' @export
batss.glm = function(
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
if (!(family %in% names(INLA::inla.models()$likelihood))){
  stop("invalid 'family' argument, see help files and inla documentation for available families")
} else {
  if (!(family %in% c("gaussian","binomial","nbinomial","poisson")))
    warning("functionality only tested for gaussian, binomial, negative binomial and poisson distributions")
}
if (!(link %in% INLA::inla.models()$likelihood[[family]]$link) || !(link %in% c("identity","log","logit","probit","robit","cauchit","loglog","cloglog")))
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
    # sequential
    if(computation!="parallel"){
        trial_r = lapply(id.seed,batss.trial,
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
    # parallel
    }else{if(computation=="parallel"){           
        # unix via forking
        if(Sys.info()[[1]]!="Windows"){
        trial_r = parallel::mclapply(id.seed,batss.trial,
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
        trial_r = parallel::parLapply(cl=cl,id.seed,batss.trial,
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
    sample   = batss.res.s1(trial_r,group=id.group$id,
                           type=c(apply(estimate[,"type",,drop=FALSE],2:3,paste0,collapse="")),
                           early=apply(estimate[,"look",,drop=TRUE]<n.look,2,all))
    scenario = batss.res.s2(sample,target=id.target$id)
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
    # sequential
    if(computation!="parallel"){
        trial_r = lapply(id.seed,batss.trial,
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
    # parallel    
    }else{if(computation=="parallel"){           
        # unix via forking
        if(Sys.info()[[1]]!="Windows"){
        trial_r = parallel::mclapply(id.seed,batss.trial,
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
        trial_r = parallel::parLapply(cl=cl,id.seed,batss.trial,
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
    }}     
    
    ##
    ## results
    ##
    #browser()
    estimate = batss.res.e(trial_r,id.target)
    tar.p    = batss.res.tp(estimate,id.target)
    tar.g    = batss.res.tg(estimate,id.target)
    eff.p    = batss.res.ep(estimate,id.target,n.look)
    eff.g    = batss.res.eg(estimate,id.target,n.look)
    fut.p    = batss.res.fp(estimate,id.target,n.look)
    fut.g    = batss.res.fg(estimate,id.target,n.look)
    sample   = batss.res.s1(trial_r,group=id.group$id,
                           type=c(apply(estimate[,"type",,drop=FALSE],2:3,paste0,collapse="")),
                           early=apply(estimate[,"look",,drop=TRUE]<n.look,2,all))
    scenario = batss.res.s2(sample,target=id.target$id)
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
class(out) = "batss"
out
}

