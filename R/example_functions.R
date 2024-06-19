# ----- ALLOCATION -----

#' @name alloc.simple
#' @title Simple allocation function
#' @description [alloc.simple] independently randomises each unit to a group (i.e., flips a coin for each unit) so that the observed allocation probabilities may be far from the target ones. This strategy is often considered to be a poor choice.
#' @param m a scalar of number of participants to be allocated.
#' @param prob a named vector of allocation ratio or probabilities.
#' @return [alloc.simple] returns an object of class \link[base]{factor} of length '`m`' with levels matching the names of the vector '`prob`'.
#' @seealso [alloc.balanced()], another group allocation function.
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
#' @param m a scalar of number of participants to be allocated.
#' @param prob a named vector of allocation ratio or probabilities.
#' @return [alloc.balanced] returns an object of class \link[base]{factor} of length '`m`' with levels matching the names of the vector '`prob`'.
#' @seealso [alloc.simple()], another group allocation function.
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
#' @description allows stopping an arm for efficacy when the probability of the corresponding target parameter being greater than `delta.eff` is greater than a fixed value `b`.
#' @param posterior the BATS ingredient '`posterior`' corresponding, in this context, to the (posterior) probability P(beta>delta.eff) for a given target parameter.
#' @param b the cut-off value to declare efficacy.
#' @return [eff.arm.simple] returns a logical constant.
#' @export
eff.arm.simple = function(posterior,b){
  posterior > b
}

#' @name eff.arm.infofract
#' @title information-fraction based arm efficacy stop 
#' @description allows stopping an arm for efficacy at a given look when the probability of the corresponding target parameter being greater than delta is greater than a function of the information fraction at that look.
#' @param posterior the BATS ingredient '`posterior`' corresponding, in this context, to the (posterior) probability P(beta>delta.eff) for a given target parameter.
#' @param b a tuning parameter (to be defined in `eff.arm.control`).
#' @param n the BATS ingredient '`n`' corresponding to the vector of number of recruited participants per arm including the control group.
#' @param N the BATS ingredient '`N' corresponding to the maximum (planned) sample size.
#' @param p a tuning parameter (to be defined in `eff.arm.control`).
#' @return [eff.arm.infofract] returns a logical constant.
#' @export
eff.arm.infofract = function(posterior,b,n,N,p){
  posterior > (1-(b*(sum(n)/N)^p))
}


# ----- EFFICACY, TRIAL -----

#' @name eff.trial.all
#' @title trial efficacy stop 
#' @description allows stopping the trial for efficacy if *all* target parameters reached efficacy at the look of interest or before.
#' @param eff.target the BATS ingredient '`eff.target`' corresponding to a \link[base]{logical} vector of the same length as argument `which` (i.e., the number of target parameters) indicating if efficacy was reached for each target parameter at that stage or at a previous stage.
#' @return [eff.trial.all] returns a logical constant.
#' @export
eff.trial.all = function(eff.target){all(eff.target)}

#' @name eff.trial.any
#' @title trial efficacy stop 
#' @description allows stopping the trial for efficacy if *at least one* target parameter reached efficacy at the look of interest.
#' @param eff.target the BATS ingredient '`eff.target`' corresponding to a \link[base]{logical} vector of the same length as argument `which` (i.e., the number of target parameters) indicating if efficacy was reached for each target parameter at that stage or at a previous stage.
#' @return [eff.trial.any] returns a logical constant.
#' @export
eff.trial.any = function(eff.target){any(eff.target)}


# ----- FUTILITY, ARMS -----

#' @name fut.arm.simple
#' @title arm futility stop 
#' @description allows stopping an arm for futility when the probability of the corresponding target parameter being greater than delta is smaller than a fixed value 'b'
#' @param posterior the BATS ingredient '`posterior`' corresponding, in this context, to the (posterior) probability P(beta>delta.fut) for a given target parameter.
#' @param b the cut-off value to declare futility (to be defined in `fut.arm.control`).
#' @return [fut.arm.simple] returns a logical constant.
#' @export
fut.arm.simple = function(posterior,b){
  posterior < b
}


# ----- FUTILITY, TRIAL -----

#' @name fut.trial.all
#' @title trial futility stop 
#' @description allows stopping the trial for efficacy if *all* active treatment reached futility at the look of interest or before.
#' @param fut.target the BATS ingredient '`fut.target`' corresponding to a \link[base]{logical} vector of the same length as argument `which` (i.e., the number of target parameters) indicating if futility was declared for each target parameter at that stage or at a previous stage.
#' @return [fut.trial.all] returns a logical constant.
#' @export
fut.trial.all = function(fut.target){all(fut.target)}


# ----- RAR -----

#' @name RAR.trippa
#' @title RAR of Trippa et al. (2012)
#' @description define the group allocation probabilities based on the response adaptive randomisation rule of Trippa et al. (2012)
#' @param posterior the BATS ingredient '`posterior`' corresponding, in this context, to the (posterior) probability P(beta>delta.RAR) for a given target parameter.
#' @param n the BATS ingredient '`n`' corresponding to the vector of number of recruited participants per arm including the control group at the look of interest.
#' @param N the BATS ingredient '`N' corresponding to the maximum (planned) sample size.
#' @param ref the BATS ingredient '`ref`' corresponding to a \link[base]{logical} vector of length "number of arms plus 1 (for the reference group)" and indicating which group is the reference one.
#' @param active the BATS ingredient '`active`' corresponding to a \link[base]{logical} vector of the same length as argument `which` (i.e., the number of target parameters) indicating if each arm is active at the look of interest.
#' @param gamma a scaling factor (to be defined in `fut.arm.control`).
#' @param eta a scaling factor (to be defined in `fut.arm.control`).
#' @param nu a scaling factor (to be defined in `fut.arm.control`).
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
#' @param active the BATS ingredient '`active`' corresponding to a \link[base]{logical} vector of the same length as argument `which` (i.e., the number of target parameters) indicating if each arm is active at the look of interest.
#' @return [RAR.optimal] returns a vector of probabilities with length of active.
#' @export
RAR.optimal = function(active){
  K <-  sum(active)-1
  tot <- K+sqrt(K)

  p = c(sqrt(K), rep(1,K))/tot
}
