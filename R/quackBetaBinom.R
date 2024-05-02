#' Model Selection - Expectation Maximization - Beta-Binomial Mixture
#'
#' @description This function uses the expectation maximization of both
#' the beta-binomial and beta-binomial-uniform mixture models for model selection.
#' Here we can run up to 32 mixture models.
#'
#' @param xm Matrix with two columns with total coverage and coverage
#'  for a randomly sampled allele.
#' @param samplename Name of sample to be included in output.
#' @param cores Threads available to run process in parallel.
#' @param parallel default = FALSE, set to true if cores > 1.
#' @param trunc List of two values representing the lower and upper bounds for
#'  allele frequency truncation , \eqn{c_{L}} and \eqn{c_{U}}. If allele frequency
#'  truncation was done to remove error, then you do not need to truncate
#'  the expected. If no truncation has been done, this should be set to c(0,0),
#'  which is the default.
#' @param tau Sequencing overdispersion parameter. If tau and error are provided,
#'  the variance of each mixture will be inferred from these values.
#'  If not, the variance by default is equal to 0.01 or 0.001.
#' @param error Sequencing error rate. If tau and error are provided,
#' the variance of each mixture will be inferred from these values. If not, the variance by default is equal to 0.01 or 0.001.
#' @param lowvar Default to FALSE. When false, variance is equal to 0.01.
#'  If set to TRUE and tau and error are not provided,
#'  the variance will be set as 0.001.
#' @param free default = FALSE, skip the free model calculation and does not
#'   calculate delta log-likelihood.
#'
#' @return BIC scores and log-likelihood (LL) mixture models including diploid,
#'  triploid, tetraploid, pentaploid, and hexaploid. When free = TRUE,
#'  the delta log-likelihood (dLL) is calculated based on
#'  the associated free model (without or with a uniform mixture).
#'  For BIC or delta-log likelihood, the smallest score is the most likely model.
#'  For LL, the largest score is the most likely model.
#'
#' @importFrom foreach foreach %dopar% %:%
#' @importFrom future plan availableCores multisession
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#'

quackBetaBinom <- function(xm, samplename, cores, parallel= FALSE, trunc = c(0.0,0.0),
                           lowvar = FALSE, tau = NA, error = NA, free = FALSE){
  cat(paste0(" \t\t <(.)__ <(.)__ <(-)__",
             "\n \t\t  (___/  (___/  (___/  nQuack-in-progress", "\n"))

  parallel <- ifelse(cores > 1, TRUE, FALSE)
  message(paste0("parallel set to ", parallel))

  # Input data setup
  xm <- as.matrix(xm)
  n <- nrow(xm)
  lnn <- log(n)

  # Input parameters
  set <- list()
  setU <- list()
  if(lowvar == FALSE){
    ## Model Params
    pu =  list(avec = c(0.11, 0.22, 0.33, 0.22, 0.11, 0.01), mvec = c(0.20, 0.33, 0.50, 0.67, 0.80), svec = c(0.01, 0.01, 0.01, 0.01, 0.01));
    p = list(avec = c(0.11, 0.22, 0.34, 0.22, 0.11), mvec = c(0.20, 0.33, 0.50, 0.67, 0.80), svec = c(0.01, 0.01, 0.01, 0.01, 0.01));

    set[[1]] =  list(avec = c(1.00), mvec = c(0.50), svec = c(0.01));
    set[[2]] =  list(avec = c(0.50, 0.50), mvec = c(0.67, 0.33), svec = c(0.01, 0.01));
    set[[3]] =  list(avec = c(0.33, 0.33, 0.33), mvec = c(0.25, 0.50, 0.75), svec = c(0.01, 0.01, 0.01));
    set[[4]] =  list(avec = c(0.25, 0.25, 0.25, 0.25), mvec = c(0.20, 0.40, 0.60, 0.80), svec = c(0.01, 0.01, 0.01, 0.01));
    set[[5]] =  list(avec = c(0.20, 0.20, 0.20, 0.20, 0.20), mvec = c(0.17, 0.33, 0.50, 0.67, 0.83), svec = c(0.01, 0.01, 0.01, 0.01, 0.01));

    setU[[1]] =  list(avec = c(0.9, 0.1), mvec = c(0.50), svec = c(0.01));
    setU[[2]] =  list(avec = c(0.45, 0.45, 0.1), mvec = c(0.67, 0.33), svec = c(0.01, 0.01));
    setU[[3]] =  list(avec = c(0.3, 0.3, 0.3, 0.1), mvec = c(0.25, 0.50, 0.75), svec = c(0.01, 0.01, 0.01));
    setU[[4]] =  list(avec = c(0.225, 0.225, 0.225, 0.225, 0.1), mvec = c(0.20, 0.40, 0.60, 0.80), svec = c(0.01,0.01, 0.01, 0.01));
    setU[[5]] =  list(avec = c(0.18, 0.18, 0.18, 0.18, 0.18, 0.1), mvec = c(0.17, 0.33, 0.50, 0.67, 0.83), svec = c(0.01, 0.01, 0.01, 0.01, 0.01));
    vv <- 0.01
  }else{
    pu =  list(avec = c(0.11, 0.22, 0.33, 0.22, 0.11, 0.01), mvec = c(0.20, 0.33, 0.50, 0.67, 0.80), svec = c(0.001, 0.001, 0.001, 0.001, 0.001));
    p = list(avec = c(0.11, 0.22, 0.34, 0.22, 0.11), mvec = c(0.20, 0.33, 0.50, 0.67, 0.80), svec = c(0.001, 0.001, 0.001, 0.001, 0.001));

    set[[1]] =  list(avec = c(1.00), mvec = c(0.50), svec = c(0.001));
    set[[2]] =  list(avec = c(0.50, 0.50), mvec = c(0.67, 0.33), svec = c(0.001, 0.001));
    set[[3]] =  list(avec = c(0.33, 0.33, 0.33), mvec = c(0.25, 0.50, 0.75), svec = c(0.001, 0.001,0.001));
    set[[4]] =  list(avec = c(0.25, 0.25, 0.25, 0.25), mvec = c(0.20, 0.40, 0.60, 0.80), svec = c(0.001,0.001, 0.001, 0.001));
    set[[5]] =  list(avec = c(0.20, 0.20, 0.20, 0.20, 0.20), mvec = c(0.17, 0.33, 0.50, 0.67, 0.83), svec = c(0.001, 0.001, 0.001, 0.001, 0.001));

    setU[[1]] =  list(avec = c(0.9, 0.1), mvec = c(0.50), svec = c(0.001));
    setU[[2]] =  list(avec = c(0.45, 0.45, 0.1), mvec = c(0.67, 0.33), svec = c(0.001, 0.001));
    setU[[3]] =  list(avec = c(0.3, 0.3, 0.3, 0.1), mvec = c(0.25, 0.50, 0.75), svec = c(0.001, 0.001, 0.001));
    setU[[4]] =  list(avec = c(0.225, 0.225, 0.225, 0.225, 0.1), mvec = c(0.20, 0.40, 0.60, 0.80), svec = c(0.001,0.001, 0.001, 0.001));
    setU[[5]] =  list(avec = c(0.18, 0.18, 0.18, 0.18, 0.18, 0.1), mvec = c(0.17, 0.33, 0.50, 0.67, 0.83), svec = c(0.001, 0.001, 0.001, 0.001, 0.001));
    vv <- 0.001
     }


  if((is.na(tau) + is.na(error)) == 2 ){
    set <- set
    setU <- setU
  }else{
    set <- setconvert(set, tau, error)
    setU <- setconvert(setU, tau, error)
  }


  # Set up parm list for BIC calculations
  parmlistN <- c()
  parmlistN[[1]] <- c(0, 2, 3, 4, 5)
  parmlistN[[2]] <- c(2, 4, 6, 8, 10) #fixed_2
  parmlistN[[3]] <- c(1, 2, 3, 4, 5)  #fixed_3

  parmlistU <- c()
  parmlistU[[1]] <- c(2, 3, 4, 5, 6)
  parmlistU[[2]] <- c(3, 5, 7, 9, 11)
  parmlistU[[3]] <- c(1, 2, 3, 4, 5)

   # EM setup
  niter = 1000
  epsilon = 0.1

  if(free == TRUE){
    type = "free"

    message("Free Model Calculations Started")
    # Free models
    LL1 <- data.frame(LL = (emstepBB(p, xm, niter, epsilon, trunc, "free"))$loglikelihood, type = "free", mixture = "free", distribution = "beta-binomial")
    LL2 <- data.frame(LL = (emstepBBU(pu, xm, niter, epsilon, trunc, "free"))$loglikelihood, type = "free", mixture = "free", distribution = "beta-binomial-uniform")

    LL1$BIC <- 0
    LL2$BIC <-  0

  }else if(free == FALSE){
    message("Free Model Skipped. Log-likelihood ratio will not be included")
  }
  # Fixed models
  if(is.na(cores)){
    cores <- (future::availableCores()-1)
  }

  if(parallel == TRUE){
    cl <- parallel::makeCluster(cores)

    ## Register parallel
    doParallel::registerDoParallel(cl)

    bvec <- c("fixed", "fixed_2", "fixed_3")
    mixturekey <- c("diploid", "triploid", "tetraploid", "pentaploid", "hexaploid")
    tt <- ifelse(sum(trunc) > 0, TRUE, FALSE)
    opts <- list(chunkSize=2)

    message("Calculating likelihood of each mixture with a beta-binomial distibution.")
    d <- foreach(type=1:3, .combine='rbind', .options.nws=opts) %:%
      foreach::foreach(iter = 1:5, .combine='rbind', .packages = "nQuack") %dopar% {
        if(iter == 1 & type == 1){
          ab <- nQuack::alphabetacalc(set[[1]]$mvec, set[[1]]$svec)
          if(tt == TRUE){
            obj <- list(loglikelihood = sum(log(extraDistr::dbbinom(xm[,2], xm[,1], ab[1,1], ab[2,1], log = FALSE)/(extraDistr::pbbinom((xm[,1]*trunc[2]), xm[,1], ab[1,1], ab[2,1]) - extraDistr::pbbinom((xm[,1]*trunc[1]), xm[,1], ab[1,1], ab[2,1])))))
          }else{
            obj <- list(loglikelihood = sum(extraDistr::dbbinom(xm[,2], xm[,1], ab[1,1], ab[2,1], log = TRUE)))

          }
        } else{
          obj <-  nQuack::emstepBB(set[[iter]], xm, niter, epsilon, trunc,  type = bvec[type])
        }

        BICout <- ((-2*obj$loglikelihood) + (lnn*parmlistN[[type]][iter]))
        data.frame(LL=obj$loglikelihood, type=bvec[type], mixture=mixturekey[iter], distribution = "beta-binomial", BIC = BICout)
      }

   message("Calculating likelihood of each mixture with a beta-binomial + uniform distibution.")
    opts <- list(chunkSize=2)
    du <- foreach(type=1:3, .combine='rbind', .options.nws=opts) %:%
      foreach::foreach(iter = 1:5, .combine='rbind', .packages = "nQuack") %dopar% {
        obj <-  nQuack::emstepBBU(setU[[iter]], xm, niter, epsilon, trunc,  type = bvec[type])
        BICout <- ((-2*obj$loglikelihood) + (lnn*parmlistU[[type]][iter]))
        data.frame(LL=obj$loglikelihood, type=bvec[type], mixture=mixturekey[iter], distribution = "beta-binomial-uniform", BIC = BICout)
      }

    # Stop cluster
    parallel::stopCluster(cl)

  }else if(parallel == FALSE){
    bvec <- c("fixed", "fixed_2", "fixed_3")
    mixturekey <- c("diploid", "triploid", "tetraploid", "pentaploid", "hexaploid")
    tt <- ifelse(sum(trunc) > 0, TRUE, FALSE)

    d <- as.data.frame(matrix(ncol = 5, nrow = 0))
    colnames(d) <- c("LL", "type", "mixture", "distribution", "BIC")
    message("Calculating likelihood of each mixture with a beta-binomial distibution.")
    for(type in 1:3){
      for(iter in 1:5){
        if(iter == 1 & type == 1){
          ab <- alphabetacalc(0.5, vv)
          if(tt == TRUE){
            obj <- list(loglikelihood = sum(log(extraDistr::dbbinom(xm[,2], xm[,1], ab[1,1], ab[2,1], log = FALSE)/(extraDistr::pbbinom((xm[,1]*trunc[2]), xm[,1], ab[1,1], ab[2,1]) - extraDistr::pbbinom((xm[,1]*trunc[1]), xm[,1], ab[1,1], ab[2,1])))))
          }else{
            obj <- list(loglikelihood = sum(extraDistr::dbbinom(xm[,2], xm[,1], ab[1,1], ab[2,1], log = TRUE)))
          }
        } else{
          obj <-  emstepBB(set[[iter]], xm, niter, epsilon, trunc,  type = bvec[type])
        }

        BICout <- ((-2*obj$loglikelihood) + (lnn*parmlistN[[type]][iter]))
        d <- rbind(d, data.frame(LL=obj$loglikelihood, type=bvec[type], mixture=mixturekey[iter], distribution = "beta-binomial", BIC = BICout))
      }
    }

    du <- as.data.frame(matrix(nrow = 0, ncol = 5))
    colnames(du) <- c("LL", "type", "mixture", "distribution", "BIC")
    message("Calculating likelihood of each mixture with a beta-binomial + uniform distibution.")

    for(type in 1:3){
      for(iter in 1:5){
        obj <-  emstepBBU(setU[[iter]], xm, niter, epsilon, trunc,  type = bvec[type])
        BICout <- ((-2*obj$loglikelihood) + (lnn*parmlistU[[type]][iter]))
        du <- rbind(du, data.frame(LL=obj$loglikelihood, type=bvec[type], mixture=mixturekey[iter], distribution = "beta-binomial-uniform", BIC = BICout))
      }
    }


  }
  if(free == TRUE){
    dLL <- c()
    for(kk in 1:nrow(d)){
     dLL[kk] <-(LL1$LL - d$LL[kk])
    }
    d$dLL <- dLL

    dULL <- c()
    for(kk in 1:nrow(du)){
      dULL[kk] <-(LL2$LL - du$LL[kk])
    }
    du$dLL <- dULL

    LL1$dLL <- 0
    LL2$dLL <- 0


    out <- rbind(LL1, LL2, d, du)

  }else if(free == FALSE){
    out <- rbind( d, du)
  }
  out$sample <- samplename
  colnames(out)
  return(out)
}


