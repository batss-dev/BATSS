# estimate matrix
bats.res.e = function(trial_r,id.target){       
    out = array(unlist(lapply(trial_r,function(x){
        type = rep(NA,nrow(x$target))
        for(i in 1:length(type)){
            eff = x$target$efficacy[i]&!is.na(x$target$efficacy[i])
            fut = x$target$futility[i]&!is.na(x$target$futility[i])
            type[i] = ifelse(eff&fut,3,ifelse(eff,1,ifelse(fut,2,0)))
        }
        c(x$target$look,type,x$target$mid)
        })),
        dim=c(nrow(id.target),3,length(trial_r)),
        dimnames=list(id.target$id,c("look","type","mid")))
    out[,"type",][is.na(out[,"type",])] = 0
    out
    }
# target per parameter
bats.res.tp = function(estimate,id.target){
    out = id.target[,c("pos","id","alternative","group")] 
    out$efficacy = apply(estimate[,"type",,drop=FALSE]==1,1,mean)
    out$futility = apply(estimate[,"type",,drop=FALSE]==2,1,mean)
    colnames(out)[1] = ""    
    out
    }
# target global   
bats.res.tg = function(estimate,id.target){
    out = data.frame(pos=1:2,id=c("At least one","All"),
                     alternative="",group="",efficacy=NA,futility=NA)
    out$efficacy[1] = mean(apply(estimate[,"type",,drop=FALSE]==1,3,sum)>0)
    out$futility[1] = mean(apply(estimate[,"type",,drop=FALSE]==2,3,sum)>0)
    out$efficacy[2] = mean(apply(estimate[,"type",,drop=FALSE]==1,3,all)>0)
    out$futility[2] = mean(apply(estimate[,"type",,drop=FALSE]==2,3,all)>0)
    colnames(out)[1] = ""    
    out
    }  
# efficacy per target parameter
bats.res.ep = function(estimate,id.target,n.look){
    out = id.target[,c("pos","id","alternative","group")]
    out$early   = apply((estimate[,"type",,drop=FALSE]==1)*(estimate[,"look",,drop=FALSE]<n.look),1,mean)
    out$last    = apply((estimate[,"type",,drop=FALSE]==1)*(estimate[,"look",,drop=FALSE]==n.look),1,mean)
    out$overall = apply((estimate[,"type",,drop=FALSE]==1),1,mean)
    colnames(out)[1] = ""    
    out
    }
# efficacy global
bats.res.eg = function(estimate,id.target,n.look){
    out = data.frame(pos=1:2,id=c("At least one","All"),
                     alternative="",group="",early=NA,last=NA,overall=NA)
    # at least one
    out$early[1] = mean(apply(estimate[,"type",,drop=FALSE]==1&(estimate[,"look",,drop=FALSE]<n.look),3,sum)>0)
    out$last[1]  = mean(apply(estimate[,"type",,drop=FALSE]==1&(estimate[,"look",,drop=FALSE]==n.look),3,sum)>0)
    out$overall[1]  = mean(apply(estimate[,"type",,drop=FALSE]==1,3,sum)>0)
    # all
    out$early[2] = mean(apply(estimate[,"type",,drop=FALSE]==1&(estimate[,"look",,drop=FALSE]<n.look),3,all)>0)
    out$last[2]  = mean(apply(estimate[,"type",,drop=FALSE]==1&(estimate[,"look",,drop=FALSE]==n.look),3,all)>0)
    out$overall[2]  = mean(apply(estimate[,"type",,drop=FALSE]==1,3,all)>0)
    # out
    colnames(out)[1] = ""    
    out
    }  
# futility per target parameter
bats.res.fp = function(estimate,id.target,n.look){
    out = id.target[,c("pos","id","alternative","group")]
    out$early   = apply((estimate[,"type",,drop=FALSE]==2)*(estimate[,"look",,drop=FALSE]<n.look),1,mean)
    out$last    = apply((estimate[,"type",,drop=FALSE]==2)*(estimate[,"look",,drop=FALSE]==n.look),1,mean)
    out$overall = apply((estimate[,"type",,drop=FALSE]==2),1,mean)
    colnames(out)[1] = ""    
    out
    }
# efficacy global
bats.res.fg = function(estimate,id.target,n.look){
    out = data.frame(pos=1:2,id=c("At least one","All"),
                     alternative="",group="",early=NA,last=NA,overall=NA)
    # at least one
    out$early[1] = mean(apply(estimate[,"type",,drop=FALSE]==2&(estimate[,"look",,drop=FALSE]<n.look),3,sum)>0)
    out$last[1]  = mean(apply(estimate[,"type",,drop=FALSE]==2&(estimate[,"look",,drop=FALSE]==n.look),3,sum)>0)
    out$overall[1]  = mean(apply(estimate[,"type",,drop=FALSE]==2,3,sum)>0)
    # all
    out$early[2] = mean(apply(estimate[,"type",,drop=FALSE]==2&(estimate[,"look",,drop=FALSE]<n.look),3,all)>0)
    out$last[2]  = mean(apply(estimate[,"type",,drop=FALSE]==2&(estimate[,"look",,drop=FALSE]==n.look),3,all)>0)
    out$overall[2]  = mean(apply(estimate[,"type",,drop=FALSE]==2,3,all)>0)
    # out
    colnames(out)[1] = ""    
    out
    }  
bats.res.s1 = function(trial_r,group,type,early){
    size = as.data.frame(matrix(unlist(lapply(trial_r,function(x,group){
        x$look[max(x$target$look,na.rm=TRUE),paste0("n(",group,")")]
        },group=group)),byrow=TRUE,ncol=length(group),
        dimnames=list(names(trial_r),group)))
    cbind(size,type,early)    
    }
bats.res.s2 = function(sample,target){
    tablew = table(sample$type)
    tablew = tablew[order(tablew,decreasing=TRUE)]
    out    = data.frame(pos=1:length(tablew),id=names(tablew),
                        overall=c(tablew)/sum(tablew),early=NA)
    early  = round(tapply(sample$early,sample$type,mean),2)
    out[names(early),"early"] = early
    for(i in 1:length(target)){
        out = cbind(out,as.numeric(substr(out$id,i,i)))
        colnames(out)[ncol(out)] = target[i]
        }
    out = cbind(out[,-c(3:4)],out[,3:4])    
    colnames(out)[1] = ""
    out
    }
## useful functions (other)
bats.trial = function(int,data,model,link,family,beta,prob0,
                     RAR,RAR.control,
                     eff.arm,eff.trial,
                     eff.arm.control,eff.trial.control,
                     fut.arm,fut.trial,
                     fut.arm.control,fut.trial.control,
                     id.target,n.target,
                     id.look,n.look,
                     id.group,n.group,groupvar,
                     twodelta,delta.eff,delta.fut,delta.RAR,
                     var,var.control,id.var,n.var,
                     #linux.os=linux.os,
                     extended,...){
# int=2
                     
    # cat(paste0("\t start:",int,"\n"))
    set.seed((n.look+1)*int)  

    # generate data for initial panel
    n = m = N = prob = ref = target = ref = active = mu = posterior = NULL
    env = new.env()
    assign("m",id.look[1,"m"], envir = env)
    assign("n",id.look[1,"n"], envir = env)
    assign("prob",prob0 , envir = env)
    assign("var",var , envir = env)
    assign("var.control",var.control , envir = env)
    assign("N", id.look$n[n.look], envir = env)  
    assign("ref",id.group$ref, envir = env)
    
    #call functions from 'var'
    pos.col <- 2                                                                             # initialize the column indicator (starting at two, column one is the response)
    for(vw in 2:length(var)) {                                                               # cycle through all variables in 'var' - starting from 2 because the first variable here is the response
      tmp_nam <- names(var)[vw]                                                              # store current variable name
      args_ <- plyr::.(n=m,m=m,prob=prob)                                                    # set function arguments, these are preset in the trial function loop
      if (tmp_nam %in% names(var.control)) args_ <- c(args_, var.control[[tmp_nam]])         # add additional arguments if specified in 'var.control'
      tmp_var <- R.utils::doCall(var[[tmp_nam]], envir = env, args = args_)                  # call variable generating function
      if (!is.matrix(tmp_var)) {                                                             # check if the generated data is NOT a matrix 
        data[, pos.col] <- tmp_var                                                           # fill column in 'data' according to position indication
        pos.col <- pos.col+1                                                                 # increase position indicator
      } else {
        colnames(tmp_var) <- paste0(tmp_nam,1:dim(tmp_var)[2])                               # name columns of matrix 'name'1,'name'2, etc
        for (jj in 1:dim(tmp_var)[2]) {
          data[, pos.col] <- tmp_var[,jj]                                                    # cycle through columns and fill 'data' accordingly
          pos.col <- pos.col+1                                                               # keep track of position in 'data'
        }
      }
    }
    #X = model.matrix(as.formula(paste0("~",strsplit(model,"~")[[1]][2])),data=data)
    X <- model.matrix(model[-2], data = data)                                            
    #---
    XB = X%*%beta
    assign("mu",switch(link,
                       "identity" = XB,
                       "log" = exp(XB),
                       "logit" = inla.link.logit(XB, inverse=TRUE),
                       "probit" = inla.link.probit(XB, inverse=TRUE),
                       "robit" = inla.link.robit(XB, inverse=TRUE),
                       "cauchit" = inla.link.cauchit(XB, inverse=TRUE),
                       "loglog" = inla.link.loglog(XB, inverse=TRUE),
                       "cloglog" = inla.link.cloglog(XB, inverse=TRUE)),envir=env)
          
    tmp_nam <- names(var)[1] 
    args_ <- plyr::.(n=m,mu=mu)                                                             # create a quoted(!) list of available 'ingredients' 
    if (!(identical(var[[1]],rbinom) || identical(var[[1]],rnbinom))) {
      names(args_)[1:2] <- formalArgs(var[[1]])[1:2]  
    } else {
      if (identical(var[[1]],rbinom)) {
        names(args_)[1:2] <- c("n","prob")
      }
    }                                                                                       # rename the list objects to the names required by specified formula (NOTE: the order of the items is set, if a function requires a different order this will not work, clever rearranging may be needed)
    if (tmp_nam %in% names(var.control)) args_ <- c(args_, var.control[[tmp_nam]])          # add extra arguments if provided
    data[, id.var[1]] = R.utils::doCall(var[[1]], args = args_, envir = env)                # execute in 'env' environment with unused arguments allowed
    #---
    #cat("B")
    # prepare
    posterior.fun = function(inf,fit,delta){
        prob = inla.pmarginal(delta, fit$marginals.fixed[[unlist(inf[1])]])
        ifelse(inf[2]=="greater",1-prob,prob)
    }    
    mx.posterior_eff.lt = mx.posterior_fut.lt = matrix(NA,nrow=n.look,ncol=n.target,dimnames=list(id.look$id,id.target$id))
    if (!is.null(RAR))  mx.posterior_RAR.lt = mx.posterior_eff.lt
    mx.futility.lt = mx.efficacy.lt = matrix(FALSE,nrow=n.look,ncol=n.target,
                             dimnames=list(id.look$id,id.target$id))
    mx.rprob.lt = matrix(NA,nrow=n.look,ncol=n.group,
                         dimnames=list(id.look$id,id.group$id))   
    #cat("C")    
    
    dots <- rlang::dots_list(...,.named=TRUE) # dots=NULL
    
    # loop
    #if(INLA::inla.os.type()=="linux"&!is.na(linux.os)){
    #    INLA::inla.binary.install(os=linux.os,verbose=TRUE,md5.check=FALSE)       
    #    }
    
    for(lw in 1:n.look){# lw=0; lw=lw+1
        # size
        #cat(.p("look:",lw,"\n"))
        temp = table(data[,groupvar])
        id.look[lw,names(temp)] = temp 
        assign("n",temp, envir = env)  
        assign("ref",id.group$ref, envir = env) 
        #cat("D")
        # fit 
        if ("control.family" %in% names(dots)) {
          control.link <- list(control.link=list(model=link))
          dots$control.family <- c(dots$control.family,control.link)
          # fit = inla(formula=model, data=data, family=family,
           #           verbose=FALSE,dots)     does nor  work like this
          fit = do.call(inla,c(list(formula=model, data=data, family=family,
                     verbose=FALSE),dots))
        } else {
          # fit = inla(formula=model, data=data, family=family,
          #            control.family=list(control.link=list(model=link)),
          #            verbose=FALSE,dots)
          fit = do.call(inla,c(list(formula=model, data=data, family=family,
                     control.family=list(control.link=list(model=link)),
                     verbose=FALSE),dots))
        }
        #cat("E")           
        # posteriors, efficacy and futility
        aw = id.target$active
        if (!is.null(eff.arm)) {
          mx.posterior_eff.lt[lw,aw] = apply(id.target[aw,c("id","alternative"),drop=FALSE],1,
                                           posterior.fun,fit=fit,delta=delta.eff[lw]) 
        } else {
          mx.posterior_eff.lt[lw,aw] = NA
        }
        if (twodelta || (is.null(eff.arm) && !is.null(fut.arm))){
            mx.posterior_fut.lt[lw,aw] = apply(id.target[aw,c("id","alternative"),drop=FALSE],1,
                                               posterior.fun,fit=fit,delta=delta.fut[lw])               
        }else{
            if (!is.null(fut.arm)) {
              mx.posterior_fut.lt[lw,aw] = mx.posterior_eff.lt[lw,aw]   
            } else {
              mx.posterior_fut.lt[lw,aw] = NA
            }
        }
        if (!is.null(RAR)) {
          mx.posterior_RAR.lt[lw,aw] = apply(id.target[aw,c("id","alternative"),drop=FALSE],1,
                                             posterior.fun,fit=fit,delta=delta.RAR[lw])
        }
        
        #cat("F")           
        # update mx.futility.lt and mx.efficacy.lt
        for(tw in 1:n.target){
            if(aw[tw]){
                # efficacy
                assign("posterior",mx.posterior_eff.lt[lw,tw], envir = env)
                    # assign("target",names(id.look[lw,names(temp)])==id.target[tw,"group"],envir = env)
                if (is.null(eff.arm) || is.na(delta.eff[lw])) {
                    mx.efficacy.lt[lw,tw] = FALSE
                  } else {
                    mx.efficacy.lt[lw, tw] = R.utils::doCall(eff.arm, args = c(plyr::.(posterior=posterior,n=n,N=N,target=target,ref=ref),eff.arm.control), envir = env)        #call function instead of parsing and evaluating string
                }
                # futility
                if(twodelta || (is.null(eff.arm) && !is.null(fut.arm))){
                    assign("posterior",mx.posterior_fut.lt[lw,tw], envir = env)
                    }
                if (is.null(fut.arm) || is.na(delta.fut[lw])) {
                    mx.futility.lt[lw,tw] = FALSE
                  } else {
                    mx.futility.lt[lw, tw] = R.utils::doCall(fut.arm, args = c(plyr::.(posterior=posterior,n=n,N=N,target=target,ref=ref),fut.arm.control), envir = env)        #call function instead of parsing and evaluating string
                  }
                #---
            }else{
                mx.efficacy.lt[lw,tw] = FALSE
                mx.futility.lt[lw,tw] = FALSE
            }
        }
        #cat("G")           
        eff.target = apply(mx.efficacy.lt[1:lw,,drop=FALSE],2,any)
        fut.target = apply(mx.futility.lt[1:lw,,drop=FALSE],2,any)

        if (!is.null(eff.arm)) eff.stop = eff.trial(eff.target) else eff.stop = FALSE 
        if (!is.null(fut.arm)) fut.stop = fut.trial(fut.target) else fut.stop = FALSE
        #---
        # efficacy       
        if(any(mx.efficacy.lt[lw,aw])){
            # identify arms
            ew = which(mx.efficacy.lt[lw,]&aw)
            # inactive arms according to eff.trial
            id.target$active[ew] = FALSE
            id.group[id.target$group[ew],"active"] = FALSE
            # save estimate and adapt list of target 
            id.target$look[ew]     = lw
            id.target$efficacy[ew] = TRUE                
            id.target[ew,c("low","mid","high")] = fit$summary.fixed[id.target$id[ew],
                                                  c("0.025quant","mean","0.975quant")]
        }
        # futility
        if(any(mx.futility.lt[lw,aw])){
            # identify arms
            fw = which(mx.futility.lt[lw,]&aw)
            # inactive arms according to fut.trial
            id.target$active[fw] = FALSE
            id.group[id.target$group[fw],"active"] = FALSE
            # save estimate and adapt list of target 
            id.target$look[fw]     = lw
            id.target$futility[fw] = TRUE                
            id.target[fw,c("low","mid","high")] = fit$summary.fixed[id.target$id[fw],
                                                  c("0.025quant","mean","0.975quant")]
        }        
        # stop trial due to no active parameters or last look
        all.stop = (eff.stop|fut.stop)|
                   all(!id.target$active)| 
                   lw==n.look        
        if(all.stop){
            if(any(id.target$active)){
                aw = which(id.target$active)
                id.target$look[aw]     = lw
                id.target[aw,c("low","mid","high")] = fit$summary.fixed[id.target$id[aw],
                                                      c("0.025quant","mean","0.975quant")]
            }
            break
        # continue
        }else{
            # prob per group
            if(!is.null(RAR)){
                # prob per group
                assign("posterior",mx.posterior_RAR.lt[lw,id.target$active], envir = env)
                assign("active",id.group$active, envir = env)
                #prob = .eval(RAR,envir=env) 
                assign("n",id.look[lw, id.group$id],envir = env)                                                 #assign ingredients to environment 'env' 
                assign("ref",id.group$ref,envir = env)
                assign("N",id.look$n[n.look],envir = env)
                assign("RAR.control", RAR.control, envir = env)
                prob = R.utils::doCall(RAR, args = c(plyr::.(posterior=posterior,n=n,N=N,ref=ref,active=active), RAR.control) ,envir = env)      #call function RAR in environment 'env'
                #---
            }else{
                prob = prob0[id.group$active]
            }
            names(prob) = id.group$id[id.group$active]
            mx.rprob.lt[lw,names(prob)] = prob/sum(prob)
            
            # predictors
            assign("n",id.look[lw+1,"n"],envir=env)
            assign("m",id.look[lw+1,"m"],envir=env)
            assign("prob",prob,envir=env)
            set.seed(lw+(n.look+1)*int)
 
            assign("var.control",var.control,envir=env)
            covar <- vector("list",length(var[-1]))
            for (ii in 1:length(var[-1])) {
              tmp_nam <- names(var)[ii+1]
              args_ <- plyr::.(n=m,m=m,prob=prob)
              if (tmp_nam %in% names(var.control)) args_ <- c(args_, var.control[[tmp_nam]])
              covar[[ii]] <- R.utils::doCall(var[[ii+1]], envir = env, args = args_)    
              if (is.matrix(covar[[ii]])) colnames(covar[[ii]]) <- paste0(tmp_nam,1:dim(covar[[ii]])[2])
            }
            
            #---
            new = as.data.frame(matrix(NA,id.look[lw+1,"m"],n.var,
                                dimnames=list(paste0(lw+1,"-",1:id.look[lw+1,"m"]),id.var)))
            
            pos.col <- 2
            for (var.count in 1:(length(var)-1)){
              if (!is.matrix(covar[[var.count]])) {
                new[,pos.col] = covar[[var.count]]
                pos.col <- pos.col+1
              } else {
                for (jj in 1:dim(covar[[var.count]])[2]) {
                  new[,pos.col] = covar[[var.count]][,jj]
                  pos.col <- pos.col+1
                }
              }
            }
                
            # response
            X = model.matrix(model[-2], data = new)                                           # call model matrix from formula object 
            #---
            if(ncol(X)!=length(beta)){"ncol(X) != length(beta)"}
            XB = X%*%beta[colnames(X)]
            assign("mu",switch(link,
                               "identity" = XB,
                               "log" = exp(XB),
                               "logit" = inla.link.logit(XB, inverse=TRUE),
                               "probit" = inla.link.probit(XB, inverse=TRUE),
                               "robit" = inla.link.robit(XB, inverse=TRUE),
                               "cauchit" = inla.link.cauchit(XB, inverse=TRUE),
                               "loglog" = inla.link.loglog(XB, inverse=TRUE),
                               "cloglog" = inla.link.cloglog(XB, inverse=TRUE)),envir=env)
            
            tmp_nam <- names(var)[1] 
            args_ <- plyr::.(n=m,mu=mu)                                                             # create a quoted(!) list of available 'ingredients' 
            if (!(identical(var[[1]],rbinom) || identical(var[[1]],rnbinom))) {
              names(args_)[1:2] <- formalArgs(var[[1]])[1:2]  
            } else {
              if (identical(var[[1]],rbinom)) {
                names(args_)[1:2] <- c("n","prob")
              } 
            }
            
            # rename the list objects to the names required by specified formula (NOTE: the order of the items is set, if a function requires a different order this will not work, clever rearranging may be needed)
            if (tmp_nam %in% names(var.control)) args_ <- c(args_, var.control[[tmp_nam]])          # add extra arguments if provided
            new[, id.var[1]] = R.utils::doCall(var[[1]], args = args_, envir = env)
            
            # data
            data = rbind(data,new)
       }# end continue
       #cat(".")
    }# end loop
    
    # output
    # cat(paste0("\t end:",int,"\n"))

    colnames(mx.rprob.lt)      = paste0("r(",colnames(mx.rprob.lt),")")
    colnames(id.look)[-c(1:4)] = paste0("n(",colnames(id.look)[-c(1:4)],")")
    colnames(mx.posterior_eff.lt)  = paste0("pe(",colnames(mx.posterior_eff.lt),")")
    colnames(mx.posterior_fut.lt)  = paste0("pf(",colnames(mx.posterior_fut.lt),")")
    list(target = id.target, look = cbind(id.look,mx.posterior_eff.lt,mx.posterior_fut.lt,mx.rprob.lt),
         data   = if(extended==2){data}else{NULL})
}

