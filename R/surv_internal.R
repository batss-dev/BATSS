batss.surv.res.s1 = function(trial_r,group,type,early){
  size = as.data.frame(matrix(unlist(lapply(trial_r,function(x,group){
    x$look[max(x$target$look,na.rm=TRUE),c(paste0("n(",group,")"),paste0("t(",group,")"),paste0("ev(",group,")"))]
  },group=group)),byrow=TRUE,ncol=length(c(paste0("n(",group,")"),paste0("t(",group,")"),paste0("ev(",group,")"))),
  dimnames=list(names(trial_r),c(paste0("n(",group,")"),paste0("t(",group,")"),paste0("ev(",group,")")))))
  t=sapply(trial_r,function(x) x$tot.time)
  cbind(size,t,type,early)
}

batss.surv.trial = function(int,data,model,family,hr,prob0,n.max,
                     RAR,RAR.control,
                     eff.arm,eff.trial,
                     eff.arm.control,eff.trial.control,
                     fut.arm,fut.trial,
                     fut.arm.control,fut.trial.control,
                     id.target,n.target,
                     id.look,n.look,
                     id.group,n.group,groupvar,
                     twodelta,delta.eff,delta.fut,delta.RAR,
                     surv,surv.control,var,var.control,
                     cens,cens.control,
                     accr,accr.control,accr.type,
                     fup,interim,
                     id.var,n.var,extended,...){

  set.seed((n.look+1)*int)
  
  #generate vector of entry times
  if (!is.null(accr)){
    if (accr.type == "random") {
      entry <- cumsum(c(0,R.utils::doCall(accr, args = c(list(n=id.look$n[n.look]-1), accr.control))))
    } else {
      entry <- c(ifelse(!is.null(accr.control$min),accr.control$min,0),
                 sort(R.utils::doCall(accr, args = c(list(n=id.look$n[n.look]-2), accr.control))), accr.control$max)
    }
    
    #calculate n and m specifically for entry times in case of interims set at time points
    if (!is.null(interim$time) && !is.na(interim$time[1])) {
      #id.look$n[1:(n.look-1)] <- unlist(foreach::foreach(ub=interim$time) %do% {sum(entry<ub)})
      id.look$n[1:(n.look-1)] <- unlist(purrr::map(interim$time, ~sum(entry < .x)))
      id.look$m <- with(id.look,c(n[1],diff(n)))
      
      data <- data[1:id.look$m[1],]
    }
  }
  
  # generate data for initial panel
  env = new.env()
  assign("m",id.look$m[1], envir = env)
  assign("prob",prob0 , envir = env)
  assign("var.control", var.control , envir = env)
  assign("surv.control",surv.control , envir = env)
  
  assign("N", id.look$n[n.look], envir = env)
  assign("ref",id.group$ref, envir = env)
  assign("fup",fup, envir = env)
  pos.col <- 3                                                                             # initialize the column indicator (starting at two, column one is the response)
  for(vw in 1:length(var)) {                                                               # cycle through all variables in 'var' - starting from 2 because the first variable here is the response
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
  
  #cat("A")
  # model matrix (ignored response side of formula object)
  X_tmp <- model.matrix(model[-2], data = data)
  X <- as.matrix(X_tmp[,-1])                         # remove intercept in surv mode
  colnames(X) <- colnames(X_tmp)[-1]
  
  X.rows <- nrow(X)                                  # prepare call to 'simsurv'
  tmp_dat <- data.frame(id=1:X.rows,X)
  tmp_betas <- log(hr)
  names(tmp_betas) <- colnames(X)
  
  assign("tmp_dat", tmp_dat, envir = env)
  assign("tmp_betas", tmp_betas, envir = env)
  
  args_ <- plyr::.(x = tmp_dat,betas = tmp_betas)
  args_ <- c(args_,surv.control)
  
  #simulate time to event data
  data[,1:2] <- R.utils::doCall(surv, alwaysArgs = args_, envir = env)[,2:3]
 
  #add censoring
  if (!is.null(cens)){
    assign("cens.control", cens.control, envir = env)
    cens <-  R.utils::doCall(cens, args = c(plyr::.(n=m), cens.control), envir = env)
    data$status[data$time>cens] <- 0
    data$time[data$time>cens] <- cens[data$time>cens]
  }
  
  #add entry times from vector 'entry'
  if (!is.null(accr)) {
    data$entry <- sample(entry[1:id.look$m[1]])
  } else {
    data$entry <- 0
  }
  
  if (!is.null(n.max)) {
    n.tot <- table(data[,names(var)[1]])
    
    maxed.out <- id.reached <- NULL
    if (any(n.max<n.tot)) {
      count <- 1
      prob_ <- prob0
      
      while (any(n.max<n.tot) && sum(prob_)!=0) {
        
        #identify groups that reached n.max 
        which.ind <- which(n.max<n.tot)
        which.group <- names(prob0)[which.ind]
        which.n <- n.max[which.ind]
        
        #identify the observation first reaching any n.max
        data_ordered <- cbind(data[order(data$entry),],count)
        n.count <- ave(data_ordered$count, data_ordered[,names(var)[1]], FUN=cumsum)
        data_ordered <- cbind(data_ordered,n.count)
        group.cutoffs <- vector()
        for (i in 1:length(which.group)) {
          group.cutoffs[i] <- which(data_ordered[,names(var)[1]]==which.group[i] & data_ordered$n.count==which.n[i])
        }
        cutoff <- min(group.cutoffs)
        id.cutoff <- which.group[which(group.cutoffs==cutoff)[1]]
  
        #drop rest and change prob
        data <- data_ordered[1:cutoff,1:(dim(data_ordered)[2]-2)]
        rownames(data) <- paste0(1, "-", 1:cutoff)
        
        prob_[names(prob_)==id.cutoff] <- 0
        if (sum(prob_)!=0) {
          #recalc probabilities
          prob_ <- prob_/sum(prob_)
          
          # generate data for rest of m
          env = new.env()
          assign("m",id.look$m[1]-cutoff, envir = env)
          assign("prob",prob_, envir = env)
          assign("var.control", var.control , envir = env)
          assign("surv.control",surv.control , envir = env)
          assign("N", id.look$n[n.look], envir = env)
          assign("ref",id.group$ref, envir = env)
          assign("fup",fup, envir = env)
          
          pos.col <- 3                                                                             # initialize the column indicator (starting at two, column one is the response)
          covar <- vector("list",length(var))
          for (ii in 1:length(var)) {
            tmp_nam <- names(var)[ii]
            args_ <- plyr::.(n=m,m=m,prob=prob_)
            if (tmp_nam %in% names(var.control)) args_ <- c(args_, var.control[[tmp_nam]])
            covar[[ii]] <- R.utils::doCall(var[[ii]], envir = env, args = args_)
            if (is.matrix(covar[[ii]])) colnames(covar[[ii]]) <- paste0(tmp_nam,1:dim(covar[[ii]])[2])
          }
          
          new   = as.data.frame(matrix(NA, id.look$m[1]-cutoff, n.var + 2,
                                       dimnames = list(paste0(1, "-", (cutoff+1):id.look$m[1]), c("time", "status", id.var))))
          pos.col <- 3
          for (var.count in 1:length(var)){
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
          
          #model matrix
          X_tmp <- model.matrix(model[-2], data = new)
          X <- matrix(X_tmp[, -1], ncol = ncol(X_tmp) - 1)
          colnames(X) <- colnames(X_tmp)[-1]
          
          #preparing call tom 'simsurv'
          X.rows <- nrow(X)
          tmp_dat <- data.frame(id = 1:X.rows, X)
          
          tmp_betas <- log(hr[colnames(X)])
          names(tmp_betas) <- colnames(X)
          
          assign("tmp_dat", tmp_dat, envir = env)
          assign("tmp_betas", tmp_betas, envir = env)
          
          args_ <- plyr::.(x = tmp_dat, betas = tmp_betas)
          args_ <- c(args_, surv.control)
          
          #generate time to event data
          new[, 1:2] <- R.utils::doCall(surv, alwaysArgs = args_, envir = env)[,2:3]
          
          #apply censoring
          if (!is.null(cens)) {
            assign("cens.control", cens.control, envir = env)
            cens <- R.utils::doCall(cens, args = c(plyr::.(n = m), cens.control), envir = env)
            new$status[new$time > cens] <- 0
            new$time[new$time > cens] <- cens[new$time > cens]
          }
          
          #add delayed entry data
          if (!is.null(accr)) {
            tmp_l <- entry[(cutoff + 1):id.look$m[1]]
            new$entry <- if (length(tmp_l) == 1) tmp_l else sample(tmp_l)
          }
          
          # data
          data = rbind(data, new)
        }
        
        maxed.out <- c(maxed.out,id.cutoff)
        n.tot <- table(data[,names(var)[1]])
      }
    }
  }
  
  if (!is.null(interim$time) && !is.null(interim$event)) {
    data_tmp <- data[order(data$time+data$entry),]
    data_tmp$cumev <- cumsum(data_tmp$status)
    
    if (interim$event.type=="all") {
      event_tmp <- data_tmp$cumev
    } else {
      if (interim$event.type=="control") {
        event_tmp <- cumsum(ifelse(data_tmp[labels(terms(model))[1]]==names(prob0)[1],data_tmp$status,0))
      } else {
        if (interim$event.type=="min") {
          tmp_count <- cumsum(ifelse(data_tmp[labels(terms(model))[1]]==names(prob0)[1],data_tmp$status,0))
          for (k in 2:length(prob0[id.group$active])) {
            tmp_count <- cbind(tmp_count,cumsum(ifelse(data_tmp[labels(terms(model))[1]]==names(prob0[id.group$active])[k],data_tmp$status,0)))
          }
          event_tmp <- apply(tmp_count,1,min)
        }
      }
    }
    
    time_tmp <- ((data_tmp$time+data_tmp$entry)[event_tmp==interim$event[1]])[1]
    data <- data_tmp[data_tmp$entry<=time_tmp,-ncol(data_tmp)]
    
    id.look$n[1] <- nrow(data)
    
    id.look$n[2:(n.look - 1)] <- unlist(purrr::map(interim$time[-1] + time_tmp, ~sum(entry < .x)))
    id.look$m <- with(id.look,c(n[1],diff(n)))
    
    interim$time_tmp <- c(time_tmp,interim$time[-1]+time_tmp)
  }
  
  if (!is.null(interim$event)){
    if (interim$event.type=="cplusone") {
      skipcounter <- 0
      event.list <- rep(list(interim$event),length(prob0)-1)
      names(event.list) <- names(prob0)[-1]
    }
  }
  
  #cat("B")
  # prepare
  posterior.fun = function(inf,fit,delta){
    prob = INLA::inla.pmarginal(delta, fit$marginals.fixed[[unlist(inf[1])]])
    ifelse(inf[2]=="greater",1-prob,prob)
  }
  mx.posterior_eff.lt = mx.posterior_fut.lt = matrix(NA,nrow=n.look,ncol=n.target,dimnames=list(id.look$id,id.target$id))
  if (!is.null(RAR))  mx.posterior_RAR.lt = mx.posterior_eff.lt
  mx.futility.lt = mx.efficacy.lt = matrix(FALSE,nrow=n.look,ncol=n.target,
                                           dimnames=list(id.look$id,id.target$id))
  mx.rprob.lt = matrix(NA,nrow=n.look,ncol=n.group,
                       dimnames=list(id.look$id,id.group$id))
  
  #cat("C")
  dots <- rlang::dots_list(...,.named=TRUE)
  
  #cat("A")
  # errorcounter <- rep(FALSE,n.look)
  
  for(lw in 1:n.look){# lw=0; lw=lw+1
    
    # size
    if (!is.null(interim$event) & is.null(interim$time) & !(lw==n.look)) {
      data_tmp <- data[order(data$time+data$entry),]
      data_tmp$cumev <- cumsum(data_tmp$status)
      
      if (interim$event.type=="all") {
        event_tmp <- data_tmp$cumev
      } else {
        if (interim$event.type=="control") {
          event_tmp <- cumsum(ifelse(data_tmp[labels(terms(model))[1]]==names(prob0)[1],data_tmp$status,0))
        } else {
          if (interim$event.type=="min") {
            tmp_count <- cumsum(ifelse(data_tmp[labels(terms(model))[1]]==names(prob0)[1],data_tmp$status,0))
            for (k in 2:length(prob0[id.group$active])) {
              tmp_count <- cbind(tmp_count,cumsum(ifelse(data_tmp[labels(terms(model))[1]]==names(prob0[id.group$active])[k],data_tmp$status,0)))
            }
            event_tmp <- apply(tmp_count,1,min)
          }
          
          else {
            if (interim$event.type=="cplusone") {
              if (sum(sapply(event.list[id.group$active[-1]],length))==0) {                   #if no combined events in control and treatment i to spend
                skipcounter <- skipcounter+1                                                  #store how many times not evaluated (for row deletion)
                id.look$n[lw] <- NA                                                           #change n,m to prepare row deletion
                id.look$m[lw] <- NA
                id.look$n[lw+1] <- nrow(data)
                id.look$m[lw+1] <- id.look$n[lw+1]-id.look$n[lw-skipcounter]
                next
              }
              tmp_obs <- na.omit(sapply(event.list, function(x) ifelse(is.null(x),NA,x[1])))            #build the vector of current combined events to check for in the data
              time_tmp_vec <- vector()                                                                  #initialize time points when next combination is reached for all treatments
              
              #calculate vector of combined events for first active treatment
              event_tmp <- cumsum(ifelse((data_tmp[labels(terms(model))[1]]==names(prob0)[1] | data_tmp[labels(terms(model))[1]]==names(prob0[id.group$active])[2]),
                                         data_tmp$status,0))
              
              #calculate time of reaching combined events specified by user
              time_tmp_vec[1] <- ifelse(length(event.list[id.group$active[-1]][[1]])==0,NA,min((data_tmp$time+data_tmp$entry)[(event_tmp==event.list[id.group$active[-1]][[1]][1])]))
              
              #repeat if more than one treatment is active
              if (sum(id.group$active)>2) {
                for (k in 3:sum(id.group$active)) {
                  event_tmp <- cumsum(ifelse((data_tmp[labels(terms(model))[1]]==names(prob0)[1] | data_tmp[labels(terms(model))[1]]==names(prob0[id.group$active])[k]),
                                             data_tmp$status,0))
                  
                  time_tmp_vec[k-1] <- ifelse(length(event.list[id.group$active[-1]][[k-1]])==0,NA,min((data_tmp$time+data_tmp$entry)[(event_tmp==event.list[id.group$active[-1]][[k-1]][1])]))
                }
              }
              time_tmp <- min(time_tmp_vec,na.rm=TRUE)                                        #find treatments reaching the set event boundaries first
              del_event <- which(time_tmp_vec==time_tmp)                                      #set vector of indices for these treatments (attention: vector length changes with active treatments)
              tmp_names <- (names(prob0)[id.group$active])[del_event+1]                       #store names of treatments reaching the event boundary
              
              #change names
              id.look$id[lw] <- rownames(id.look)[lw] <- paste0("n(e_",names(prob0)[1],(names(prob0)[id.group$active])[del_event[1]+1],
                                                                ")=",event.list[id.group$active[-1]][[del_event[1]]][1])
              for (k in del_event) {
                event.list[id.group$active[-1]][[k]] <- event.list[id.group$active[-1]][[k]][-1]    #remove used boundary from the treatment group that reached it
              }
            }
          }
        }
      }
      if (interim$event.type!="cplusone") {
        time_tmp <- ((data_tmp$time+data_tmp$entry)[(event_tmp==interim$event[lw])])[1]
      }
      
      data <- data_tmp[data_tmp$entry<=time_tmp,-ncol(data_tmp)]
      id.look$n[lw] <- nrow(data)
      id.look$m[lw]  <-  id.look$n[lw]- ifelse(lw==1,0,id.look$n[lw-1])
      id.look$m[lw+1] <- id.look$n[lw+1]-id.look$n[lw]
    }
    
    if (!is.null(n.max)){
    #calculate cut groups
      n.tot <- table(data[,names(var)[1]])
      
      #identify groups that reached n.max after data being truncated because of interim strategy
      which.ind <- which(n.max==n.tot)
      id.reached <- names(prob0)[which.ind] 
    }
    
    temp = table(data[,groupvar])
    id.look[lw,names(temp)] = temp
    
    assign("n",temp, envir = env)
    assign("ref",id.group$ref[id.group$active], envir = env)
    assign("interim$time",interim$time, envir = env)
    assign("interim$event",interim$event, envir = env)
    
    #cat("D")
    # fit
    if (!(lw==n.look)) {  # change times for looks prior to final
      
      data_calc <- data
      
      #update times and events
      #if (!is.null(accr)) {
      if (!is.null(interim$time)){
        if (!is.null(interim$event)) {
          data_calc$time <- ifelse(data$time+data$entry <= interim$time_tmp[lw],data$time,interim$time_tmp[lw]-data$entry)
          data_calc$status <- ifelse(data$time+data$entry <= interim$time_tmp[lw],data$status,0)
        } else {
          data_calc$time <- ifelse(data$time+data$entry <= interim$time[lw],data$time,interim$time[lw]-data$entry)
          data_calc$status <- ifelse(data$time+data$entry <= interim$time[lw],data$status,0)
        }
      } else {
        if (!is.null(interim$event)){
          data_calc$time <- ifelse(data$time+data$entry <= time_tmp,data$time,time_tmp-data$entry)
          data_calc$status <- ifelse(data$time+data$entry <= time_tmp,data$status,0)
        } else {
          data_calc$time <- ifelse(data$time+data$entry <= max(data$entry),data$time,max(data$entry)-data$entry)
          data_calc$status <- ifelse(data$time+data$entry <= max(data$entry),data$status,0)
        }
      }
      #}

      #calculate observed times
      id.look[lw,"t(n)"] <- sum(data_calc$time)
      id.look[lw,paste0("t(",names(temp),")")] <- aggregate(reformulate(names(var)[1],response="time"),FUN=sum,data=data_calc)[,2]
      
      #calculate observed event
      id.look[lw,"ev(n)"] <- sum(data_calc$status)
      id.look[lw,paste0("ev(",names(temp),")")] <- aggregate(reformulate(names(var)[1],response="status"),FUN=sum,data=data_calc)[,2]
      
      #calculate study duration
      id.look[lw,"t"] <- ifelse(!is.null(interim$time),ifelse(!is.null(interim$event),interim$time_tmp[lw],interim$time[lw]),
                                ifelse(!is.null(interim$event),time_tmp,entry[id.look$n[lw]]))
      
      #fit model
      fit <- do.call(INLA::inla,c(list(formula = model, family = family, data=data_calc, verbose=FALSE),dots))
      
      # tryCatch(
      #   {
      #     fit <- inla(model, family = family, data=data_calc, ...,verbose=FALSE)
      #   },
      #   error = function(e) {
      #     errorcounter[lw] <- TRUE
      #   }
      # )
    } else {    #final look
      
      if (exists("maxt",where=surv.control) && is.null(surv.control$maxt)) {
        data_final <- data
        data$time <- ifelse(data_final$time+data_final$entry <= max(data_final$entry)+fup,data_final$time,max(data_final$entry)+fup-data_final$entry)
        data$status <- ifelse(data_final$time+data_final$entry <= max(data_final$entry)+fup,data_final$status,0)
      }
      
      #calculate observed times
      id.look[lw,"t(n)"] <- sum(data$time)
      id.look[lw,paste0("t(",names(temp),")")] <- aggregate(reformulate(names(var)[1],response="time"),FUN=sum,data=data)[,2]
      
      #calculate observed events
      id.look[lw,"ev(n)"] <- sum(data$status)
      id.look[lw,paste0("ev(",names(temp),")")] <- aggregate(reformulate(names(var)[1],response="status"),FUN=sum,data=data)[,2]
      
      #calculate study duration
      id.look[lw,"t"] <- max(data$entry+data$time)
      
      #fit model
      fit <- do.call(INLA::inla,c(list(formula = model, family = family, data=data, verbose=FALSE),dots))
      
      # tryCatch(
      #   {
      #     fit <- inla(model, family = family, data=data, ...,verbose=FALSE)
      #   },
      #   error = function(e) {
      #     errorcounter[lw] <- TRUE
      #   }
      # )
    }
    
    #cat("E")
    # posteriors, efficacy and futility
    aw = id.target$active
    if (!is.null(interim$event)) {if (interim$event.type=="cplusone" && lw!=n.look) aw[which(names(prob0[-1])!=tmp_names)] = FALSE}       #in case of the 'control plus treatment i'-method, change the vector of evaluation accordingly
    if (all(aw==FALSE)) next                                                                         #for empty rows
    
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
        assign("target",names(id.look[lw,names(temp)])==id.target[tw,"group"],
               envir = env)
        assign("curr.look",lw,envir = env)
        assign("n.look",n.look,envir = env)
        assign("posterior",mx.posterior_eff.lt[lw,tw], envir = env)
        if (is.null(eff.arm) || is.na(delta.eff[lw])) {
          mx.efficacy.lt[lw,tw] = FALSE
        } else {
          mx.efficacy.lt[lw, tw] = R.utils::doCall(eff.arm, args = c(plyr::.(posterior=posterior,n=n,N=N,target=target,ref=ref,curr.look=curr.look,n.look=n.look),eff.arm.control), envir = env)        #call function instead of parsing and evaluating string
        }
        #---
        if (twodelta || (is.null(eff.arm) && !is.null(fut.arm))){
          assign("posterior",mx.posterior_fut.lt[lw,tw], envir = env)
        }
        if (is.null(fut.arm) || is.na(delta.fut[lw])) {
          mx.futility.lt[lw,tw] = FALSE
        } else {
          mx.futility.lt[lw, tw] = R.utils::doCall(fut.arm, args = c(plyr::.(posterior=posterior,n=n,N=N,target=target,ref=ref,curr.look=curr.look,n.look=n.look),fut.arm.control), envir = env)        #call function instead of parsing and evaluating string
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
    
    if (!is.null(eff.arm)) eff.stop = eff.trial(eff.target) else eff.stop = FALSE                 #call function directly (I don't think this would have worked with a user specified function before)
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
    
    # #set maxed out arms to inactive
    # if (!is.null(n.max) && !is.null(id.reached)) {
    #   id.target$active[which(is.element(id.target$group,id.reached))] = FALSE
    #   id.group[id.target$group[which(is.element(id.target$group,id.reached))],"active"] = FALSE
    # }
    
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
      if (id.look$m[lw+1]!=0) {  #if not last look
        # prob per group
        if(!is.null(RAR)){
          # prob per group
          assign("posterior",mx.posterior_RAR.lt[lw,id.target$active], envir = env)
          assign("active",id.group$active, envir = env)
          #prob = .eval(RAR,envir=env)
          assign("n",unlist(id.look[lw, id.group$id]),envir = env)                                         #assign ingredients to environment 'env'
          assign("n.ev",id.look[lw,paste0("ev(",id.group$id,")")],envir = env) 
          assign("ref",id.group$ref,envir = env)
          assign("N",id.look$n[n.look],envir = env)
          assign("RAR.control", RAR.control, envir = env)
          prob = R.utils::doCall(RAR, args = c(plyr::.(posterior=posterior,n=n, n.ev=n.ev,N=N,ref=ref,active=active), RAR.control), envir = env)      #call function RAR in environment 'env'
        }else{
          prob = prob0[id.group$active]
        }
        names(prob) = id.group$id[id.group$active]
        
        if(!is.null(n.max)) prob[is.element(names(prob),id.reached)] <- 0
        
        if (sum(prob,na.rm=TRUE)!=0) prob <- prob/sum(prob)
        mx.rprob.lt[lw,names(prob)] = prob
          
        if (sum(prob,na.rm=TRUE)!=0) {
          
          # predictors
          assign("n", id.look[lw + 1, "n"], envir = env)
          assign("m", id.look[lw + 1, "m"], envir = env)
          assign("prob", prob, envir = env)
          set.seed(lw + (n.look + 1) * int)
          assign("var.control", var.control, envir = env)
          assign("surv.control", surv.control, envir = env)
          assign("fup",fup, envir = env)
          
          covar <- vector("list",length(var))
          for (ii in 1:length(var)) {
            tmp_nam <- names(var)[ii]
            args_ <- plyr::.(n=m,m=m,prob=prob)
            if (tmp_nam %in% names(var.control)) args_ <- c(args_, var.control[[tmp_nam]])
            covar[[ii]] <- R.utils::doCall(var[[ii]], envir = env, args = args_)
            if (is.matrix(covar[[ii]])) colnames(covar[[ii]]) <- paste0(tmp_nam,1:dim(covar[[ii]])[2])
          }
          
          new   = as.data.frame(matrix(NA, id.look[lw + 1, "m"], n.var + 2,
                                       dimnames = list(paste0(lw + 1, "-", 1:id.look[lw + 1, "m"]), c("time", "status", id.var))))
          pos.col <- 3
          for (var.count in 1:length(var)){
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
          
          #model matrix
          X_tmp <- model.matrix(model[-2], data = new)
          X <- matrix(X_tmp[, -1], ncol = ncol(X_tmp) - 1)
          colnames(X) <- colnames(X_tmp)[-1]
          
          #preparing call tom 'simsurv'
          X.rows <- nrow(X)
          tmp_dat <- data.frame(id = 1:X.rows, X)
          
          tmp_betas <- log(hr[colnames(X)])
          names(tmp_betas) <- colnames(X)
          
          assign("tmp_dat", tmp_dat, envir = env)
          assign("tmp_betas", tmp_betas, envir = env)
          
          args_ <- plyr::.(x = tmp_dat, betas = tmp_betas)
          args_ <- c(args_, surv.control)
          
          #generate time to event data
          new[, 1:2] <- R.utils::doCall(surv, alwaysArgs = args_, envir = env)[,2:3]
          
          #apply censoring
          if (!is.null(cens)) {
            assign("cens.control", cens.control, envir = env)
            cens <- R.utils::doCall(cens, args = c(plyr::.(n = m), cens.control), envir = env)
            new$status[new$time > cens] <- 0
            new$time[new$time > cens] <- cens[new$time > cens]
          }
          
          #add delayed entry data
          if (!is.null(accr)) {
            tmp_l <- entry[(id.look$n[lw] + 1):id.look$n[lw + 1]]
            new$entry <- if (length(tmp_l) == 1) tmp_l else sample(tmp_l)
          }
          
          # data
          data = rbind(data, new)
          
          
          #control for max n's
          if (!is.null(n.max)) {
            n.tot <- table(data[,names(var)[1]])
            prob_ <- prob
            
            if (any(n.max<n.tot)) {
              while (any(n.max<n.tot) && sum(prob_)!=0) {
                count <- 1
                
                #identify groups that reached n.max 
                which.ind <- which(n.max<n.tot)
                which.group <- names(prob0)[which.ind]
                which.n <- n.max[which.ind]
                
                #identify the observation first reaching any n.max
                data_ordered <- cbind(data[order(data$entry),],count)
                n.count <- ave(data_ordered$count, data_ordered[,names(var)[1]], FUN=cumsum)
                data_ordered <- cbind(data_ordered,n.count)
                group.cutoffs <- rep(NA,length(which.group))
                for (i in 1:length(which.group)) {
                  group.cutoffs[i] <- which(data_ordered[,names(var)[1]]==which.group[i] & data_ordered$n.count==which.n[i])
                }
                
                #identify group reaching max first
                cutoff <- min(group.cutoffs)
                id.cutoff <- which.group[which(group.cutoffs==cutoff)]
                
                #drop rest and change prob
                data <- data_ordered[1:cutoff,1:(dim(data_ordered)[2]-2)]
  
                prob_[names(prob_)==id.cutoff] <- 0
                if (sum(prob_)!=0) {
                  prob_ <- prob_/sum(prob_)
                
                  #generate data for rest of m
                  env = new.env()
                  assign("n", id.look[lw + 1, "n"], envir = env)
                  assign("m",id.look[lw + 1, "n"]-cutoff, envir = env)
                  assign("prob",prob_, envir = env)
                  set.seed(lw + (n.look + 1) * int)
                  assign("var.control", var.control , envir = env)
                  assign("surv.control",surv.control , envir = env)
                  assign("N", id.look$n[n.look], envir = env)
                  assign("ref",id.group$ref, envir = env)
                  assign("fup",fup, envir = env)
                  
                  pos.col <- 3                                                                             # initialize the column indicator (starting at two, column one is the response)
                  covar <- vector("list",length(var))
                  for (ii in 1:length(var)) {
                    tmp_nam <- names(var)[ii]
                    args_ <- plyr::.(n=m,m=m,prob=prob_)
                    if (tmp_nam %in% names(var.control)) args_ <- c(args_, var.control[[tmp_nam]])
                    
                    covar[[ii]] <- R.utils::doCall(var[[ii]], envir = env, args = args_)
                    if (is.matrix(covar[[ii]])) colnames(covar[[ii]]) <- paste0(tmp_nam,1:dim(covar[[ii]])[2])
                  }
                  
                  new   = as.data.frame(matrix(NA, id.look[lw + 1, "n"]-cutoff, n.var + 2,
                                               dimnames = list(paste0(lw + 1, "-", (cutoff+1):id.look[lw + 1, "n"]), c("time", "status", id.var))))
                  
                  pos.col <- 3
                  for (var.count in 1:length(var)){
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
                  
                  #model matrix
                  X_tmp <- model.matrix(model[-2], data = new)
                  X <- matrix(X_tmp[, -1], ncol = ncol(X_tmp) - 1)
                  colnames(X) <- colnames(X_tmp)[-1]
                  
                  #preparing call tom 'simsurv'
                  X.rows <- nrow(X)
                  tmp_dat <- data.frame(id = 1:X.rows, X)
                  
                  tmp_betas <- log(hr[colnames(X)])
                  names(tmp_betas) <- colnames(X)
                  
                  assign("tmp_dat", tmp_dat, envir = env)
                  assign("tmp_betas", tmp_betas, envir = env)
                  
                  args_ <- plyr::.(x = tmp_dat, betas = tmp_betas)
                  args_ <- c(args_, surv.control)
                  
                  #generate time to event data
                  new[, 1:2] <- R.utils::doCall(surv, alwaysArgs = args_, envir = env)[,2:3]
                  
                  #apply censoring
                  if (!is.null(cens)) {
                    assign("cens.control", cens.control, envir = env)
                    cens <- R.utils::doCall(cens, args = c(plyr::.(n = m), cens.control), envir = env)
                    new$status[new$time > cens] <- 0
                    new$time[new$time > cens] <- cens[new$time > cens]
                  }
                  
                  #add delayed entry data
                  if (!is.null(accr)) {
                    tmp_l <- entry[(cutoff + 1):id.look$n[lw + 1]]
                    new$entry <- if (length(tmp_l) == 1) tmp_l else sample(tmp_l)
                  }
                  
                  # data
                  data = rbind(data, new) 
                }
                
                maxed.out <- c(maxed.out,id.cutoff)
                n.tot <- table(data[,names(var)[1]])
              }
            }
          }
        }

        
      }
    }# end continue
  }# end loop
  
  t.trial <- max(id.look[,"t"],na.rm=TRUE)
  
  colnames(mx.rprob.lt)      = paste0("r(",colnames(mx.rprob.lt),")")
  colnames(id.look)[colnames(id.look)%in%names(temp)] = paste0("n(",colnames(id.look)[colnames(id.look)%in%names(temp)],")")
  colnames(mx.posterior_eff.lt)  = paste0("pe(",colnames(mx.posterior_eff.lt),")")
  colnames(mx.posterior_fut.lt)  = paste0("pf(",colnames(mx.posterior_fut.lt),")")
  
  tmp_look <- cbind(id.look,mx.posterior_eff.lt,mx.posterior_fut.lt,mx.rprob.lt)

  # correct number of interims as they are dependent on the trial progression in case of "event" and "cplusone"
  if (!is.null(interim$event)) {
    if (interim$event.type=="cplusone") {
      if (any(is.na(tmp_look$n))) {
        
        tmp_look <- tmp_look[!is.na(tmp_look$n),]
        tmp_look$pos[tmp_look$id=="final"] <- dim(tmp_look)[1]
        id.target$look[id.target$look==n.look] <- dim(tmp_look)[1]
      }
    }
  }
  
  list(target = id.target, look = tmp_look, tot.time = t.trial, last.rec = max(data$entry), data = if(extended==2){data}else{NULL})
}


utils::globalVariables(c("m", "posterior", "n", "N", "target", "ref", "curr.look", "n.ev", "active"))
