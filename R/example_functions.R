# ----- ALLOCATION -----

#' @name alloc.simple
#' @title Simple allocation function
#' @description [alloc.simple] independently randomises each unit to a group (i.e., flips a coin for each unit) so that the observed allocation probabilities may be far from the target ones.
#' @param m scalar of number of participants to be allocated
#' @param prob named vector of allocation ratio or probabilities
#' @return [alloc.simple] returns an object of class \link[base]{factor} with levels as provided in the input vector 'prob' and length 'm'.
#' @examples
#' alloc.simple(100, prob = c(A=.4,B=.6))
#' table(alloc.simple(100, prob = c(A=.4,B=.6)))
#' table(alloc.simple(100, prob = c(A=.4,B=.6)))
#' @export
alloc.simple = function(m,prob){
  prob = abs(prob)/sum(abs(prob)) # normalise probabilities
  factor(sample(names(prob), m, replace=TRUE, prob = prob), levels=names(prob))
}

#' @name alloc.balanced
#' @title Balanced allocation function
#' @description [alloc.balanced] first allocates the largest possible number of units to the different groups given their exact target probabilities and then assigns randomly the remaining units to the different groups according to multinomial draws. This method leads to observed allocation probabilities matching the target ones when m*prob is an integer for each group and to observed allocation probabilities (on average) closer to the target ones compared to [alloc.simple].
#' @param m scalar of number of participants to be allocated
#' @param prob named vector of allocation ratio or probabilities
#' @return [alloc.balanced] returns an object of class \link[base]{factor} with levels as provided in the input vector 'prob' and length 'm'.
#' @examples
#' alloc.balanced(100, prob = c(A=.4,B=.6))
#' table(alloc.balanced(100, prob = c(A=.4,B=.6)))
#' table(alloc.balanced(100, prob = c(A=.4,B=.6)))
#' @export
alloc.balanced = function(m,prob){
  prob = abs(prob)/sum(abs(prob))
  m0.g = floor(prob*m)
  m0   = sum(m0.g)
  factor(rep(names(prob),m0.g+rmultinom(1,m-m0,prob)),
         levels=names(prob))
}


# ----- EFFICACY, ARMS -----

#' @name eff.arm.simple
#' @title Simple arm efficacy stop 
#' @description allows stopping an arm for efficacy when the probability of the corresponding target parameter being greater than delta is greater than a fixed value b
#' @param posterior posterior probability of P(theta>delta) > b
#' @param b probability boundary for efficacy
#' @return [eff.arm.simple] returns a logical constant.
#' @export
eff.arm.simple = function(posterior,b){
  posterior > b
}

#' @name eff.arm.infofract
#' @title information-fraction based arm efficacy stop 
#' @description allows stopping an arm for efficacy at a given look when the probability of the corresponding target parameter being greater than delta is greater than a function of the information fraction at that look.
#' @param posterior posterior probability of P(theta>delta) > b
#' @param b probability boundary for efficacy
#' @param n number of recruited participants
#' @param N total (planned) sample size
#' @param p tuning parameter
#' @return [eff.arm.infofract] returns a logical constant.
#' @export
eff.arm.infofract = function(posterior,b,n,N,p){
  posterior > (1-(b*(sum(n)/N)^p))
}


# ----- EFFICACY, TRIAL -----

#' @name eff.trial.all
#' @title trial efficacy stop 
#' @description allows stopping the trial for efficacy at a given look if all active treatments reached efficacy
#' @param eff.target logical vector with length of active intervention arms
#' @return [eff.trial.all] returns a logical constant.
#' @export
eff.trial.all = function(eff.target){all(eff.target)}

#' @name eff.trial.any
#' @title trial efficacy stop 
#' @description allows stopping the trial for efficacy at a given look if at least one active treatment reached efficacy
#' @param eff.target logical vector with length of active intervention arms
#' @return [eff.trial.any] returns a logical constant.
#' @export
eff.trial.any = function(eff.target){any(eff.target)}


# ----- FUTILITY, ARMS -----

#' @name fut.arm.simple
#' @title arm futility stop 
#' @description allows stopping an arm for futility when the probability of the corresponding target parameter being greater than delta is smaller than a fixed value 'b'
#' @param posterior posterior probability of P(theta>delta) < b
#' @param b probability boundary for futility
#' @return [fut.arm.simple] returns a logical constant.
#' @export
fut.arm.simple = function(posterior,b){
  posterior < b
}


# ----- FUTILITY, TRIAL -----

#' @name fut.trial.all
#' @title trial futility stop 
#' @description allows stopping the trial for futility at a given look if all active treatments reached efficacy
#' @param fut.target logical vector with length of active intervention arms
#' @return [fut.trial.all] returns a logical constant.
#' @export
fut.trial.all = function(fut.target){all(fut.target)}


# ----- RAR -----

#' @name RAR.trippa
#' @title RAR of Trippa et al. (2012)
#' @description define the group allocation probabilities based on the response adaptive randomisation rule of Trippa et al. (2012)
#' @param posterior posterior probability of P(theta>delta) < b
#' @param n current sample size
#' @param N total sample size
#' @param ref vector indicating the reference group
#' @param active vector indicating which treatments are active
#' @param gamma scaling factor
#' @param eta scaling factor
#' @param nu scaling factor
#' @return [RAR.trippa] returns a vector of probabilities with length of active.
#' @export
RAR.trippa = function(posterior,n,N,ref,active,gamma,eta,nu){
  g = sum(active)
  h = gamma*(sum(n)/N)^eta
  p = rep(NA,g)

  # reference
  p[1] = (exp(max(n[!ref])-n[ref])^nu)/(g-1)
  # targets
  p[2:g] = (posterior^h)/(sum(posterior^h))
  #
  unlist(p)
}

#' @name RAR.optimal
#' @title 'Optimal' control allocation
#' @description technically not response adaptive but keeps allocation ratio to control at the square root of active intervention arms
#' @param active vector indicating which treatments are active
#' @return [RAR.optimal] returns a vector of probabilities with length of active.
#' @export
RAR.optimal = function(active){
  K <-  sum(active)-1
  tot <- K+sqrt(K)

  p = c(sqrt(K), rep(1,K))/tot
}
