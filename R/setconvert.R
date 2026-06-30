#' Calculate Variance from Mean, Tau, and Sequencing Error
#'
#' @description This function is used to replace variance in mixture model sets.
#'
#' @param set A list of lists, each of the lists must contain avec, mvec, and svec.
#' @param tau Sequence overdispersion parameter for read counts.
#' @param error Sequencing error rate.
#'
#' @examples
#'  set <- c()
#'  set[[1]] =  list(avec = c(1.00), mvec = c(0.50), svec = c(0.01));
#'  set[[2]] =  list(avec = c(0.50, 0.50), mvec = c(0.67, 0.33), svec = c(0.01, 0.01));
#'  exset <- setconvert(set, tau = 0.01, error = 0.001)
#' @returns Mean and variance for the associated tau and error.
#' @export
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
