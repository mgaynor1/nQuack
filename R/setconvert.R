#' Calculate Variance from Mean, Tau, and Sequencing Error
#'
#' @description This function is used to replace variance in mixture model sets.
#'
#' @param set A list of lists, each of the lists must contain avec, mvec, and svec.
#' @param tau Sequence overdispersion parameter for read counts.
#' @param error Sequencing error rate.
#'
#' @returns Mean and variance for the associated tau and error.
setconvert <- function(set, tau, error){
  setout <- list()
  for(i in 1:length(set)){
    varout <- muvarcalcvec(set[[i]]$mvec, tau, error)
    setout[[i]] <- list(avec = c(set[[i]]$avec),
                        mvec = c(set[[i]]$mvec),
                        svec = c(varout$var))

  }
  return(setout)
}
