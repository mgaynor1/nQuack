#' Bootstrapping - Expected Maximization - Choose your distribution and type
#'
#' @description This function was made to assist with bootstrap replication for a set of models run a subset of models based on a selected distribution and type.
#' There are many limitations to this function to make this tractable, as there are 128 models that could be run with our package.
#' Here we do not include models or comparisons we found unhelpful, this includes the nQuire implementation and log-likelihood ratio tests.
#'
#' @inheritParams bestquack
#' @param nboots Number of bootstrap replicates to examine.
#'
#' @return BIC scores and log-likelihood (LL) for included mixture models.
#'  For both, the smallest score is the most likely model.


quackNboots <- function(xm, nboots = 100, distribution, type, uniform,
                      mixtures = c("diploid", "triploid", "tetraploid", "hexaploid", "pentaploid"),
                      samplename,
                      trunc = c(0.0,0.0),  lowvar = FALSE,
                      tau = NA, error = NA){

  og <- bestquack(xm, distribution, type, uniform, mixtures, samplename, trunc, lowvar, tau, error)
  ogpick <- (quackit(og))$winnerBIC
  message(paste0("The best pick for the original data is ", ogpick))

  bootits <- c()
  for(i in 1:nboots){
    bootits[[i]] <- resample_xm(xm, nrow(xm))
  }

  bootscheck <- c()
  for(i in 1:nboots){
  check <- bestquack(bootits[[i]], distribution, type, uniform, mixtures, samplename, trunc, lowvar, tau, error)
  bootscheck[i] <- (quackit(check))$winnerBIC
  }

 results <- data.frame(matrix(nrow = 2, ncol = length(mixtures)))
 colnames(results) <- mixtures
 rownames(results) <- c("original", "bootstrap.replicates")
 results[1, ogpick] <- 1

  out <- data.frame(table(bootscheck))
  for(j in 1:nrow(out)){
    results[2, as.character(out[j, 1])] <- as.numeric(out[j, 2])
}
results$sample <- samplename
  return(results)

}
