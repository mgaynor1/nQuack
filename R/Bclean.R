#' Remove noise with the beta distribution
#'
#' @description Here we filter allele frequencies with a beta mixture model that contains 5 mixtures:
#' three mixtures representing cytotypes included in nQuack and two mixtures representing a U-shaped distribution.
#' We constrained the first three mixtures to have shape and scale parameters above 1,
#' while the last two mixtures shape and scale are constrained to be less than 1.
#' With this implementation of expectation-maximization, we utilizes the scaled probability of each data point belonging
#' to each mixture model to remove site where the probability of belonging to a U-shaped mixture is higher than the
#' probability of belonging to any other mixture. Due to the computational time needed
#' to run the expectation-maximization algorithm, by default, we simple calculate this probability matrix with the E-step and do not
#' run the complete algorithm.
#'
#' @param xm  Matrix with total coverage and coverage for a randomly sampled allele.
#' @param plot Default to TRUE. The plots do not share the same y-axis, so careful interpretation is key.
#'  Warning, if nothing is removed, the plot of removed data will be missing.
#' @param quick Default to TRUE. If set as FALSE, the expectation-maximization algorithm will be run in full.
#'
#' @return Numeric matrix with total coverage and coverage for a randomly sampled allele.
#'
#' @importFrom graphics hist par
#' @importFrom stats runif
#'

Bclean <- function(xm, plot = TRUE, quick = TRUE){
  # Setup data
  xm <- as.matrix(xm)
  xi <- xm[,2]/xm[,1]

  cL <- min(xi)
  cU <- max(xi)
  # Model

  #tcalcfull <- alphabetacalcvec(mu = c(0.20, 0.33, 0.50, 0.67, 0.80), var = c(0.01, 0.01, 0.01, 0.01, 0.01))
  #set <-  list(avec = c(0.11, 0.22, 0.33, 0.22, 0.11, 0.01), t1vec = c(tcalc[,1], 0.5), t2vec = c(tcalc[,2], 0.33))

  tcalc <- alphabetacalcvec(mu = c(0.287, 0.50, 0.713), var = c(0.01, 0.01, 0.01))
  set <-  list(avec = c(0.25, 0.25, 0.25, 0.125, 0.125), t1vec = c(tcalc[,1], 0.5, 0.33), t2vec = c(tcalc[,2], 0.33, 0.5))

  if(quick == FALSE){
    message("Stay tune: Starting expectation-maximization algorithm now! nQuack calculations may take a long time to complete.")
    checkB <- emstepB3(set, xi, 1000, 0.1, c(cL,cU))
    pir <- checkB$pir
  }else if(quick == TRUE){
    checkB <- estepB3(set, xi, c(cL,cU))
    pir <- checkB$zprob
  }

  # Make Remove Set
  removelist <- c()

  for(i in 1:length(xi)){
    if(which(pir[i,] == max(pir[i,])) == 4 |which(pir[i,] == max(pir[i,])) == 5 ){
      removelist[i] <- 0
    }else{
      removelist[i] <- 1
    }
  }

  if(min(xi) > 0){
    os = round(min(xi)*100)
    size = 100-(2*os)
  }else if(min(xi) == 0 | is.na(min(xi))){
    os = 1
    size = 100
  }

  if(plot == TRUE){
    xikeep <- xi[which(removelist == 1)]
    xiremove <- xi[which(removelist == 0)]

    par(mfrow = c(1, 3))

    hist(xi,xlim = c(0, 1), main = "Original Data", xlab = "Allele Frequency", breaks = size )

    if(length(xikeep) == 0){
      plot(0, type="n", axes=FALSE, xlab="", ylab="")
    }else{
      hist(xikeep, xlim = c(0, 1), main = "Data Retained - nQuack", xlab = "Allele Frequency", breaks = size)
    }

    if(length(xiremove) == 0){
      plot(0, type="n", axes=FALSE, xlab="", ylab="")
    } else{
      hist(xiremove,  xlim = c(0, 1), main = "Data Removed - nQuack",  xlab = "Allele Frequency", breaks = size)
    }
  }

  xmout <- xm[which(removelist == 1),]
  return(xmout)
}
