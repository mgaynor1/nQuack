#' Variance calculation from Mean, Tau, and Sequencing Error
#'
#' @description This function is used to calculate variance.
#'
#' @param mu Vector of means.
#' @param tau Sequence overdispersion parameter for read counts.
#' @param error Sequencing error rate.
#'
#' @returns Mean and variance for the associated tau and error.

muvarcalcvec <- function(mu, tau, error){
  var.out <- c()
  for(z in 1:length(mu)){
    m <- ((mu[z]*(1-error))+((1-mu[z])*error))
   var.out[z] <- (m*(1-m)*tau)
  }
  parm.list <- list(var = var.out)
  return(parm.list)
}
