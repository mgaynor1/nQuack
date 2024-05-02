#' Model Selection - Expectation Maximization - Choose your distribution and type
#'
#' @description This function was made to run a subset of models based on a selected distribution and type.
#' There are many limitations to this function to make this tractable, as there are 128 models that could be run with our package.
#' Here we do not include models or comparisons we found unhelpful, this includes the nQuire implementation and log-likelihood ratio tests.
#'
#' @param xm Matrix with two columns with total coverage and coverage
#'  for a randomly sampled allele.
#' @param distribution May be set to normal, beta, or beta-binomial. We do not include the implementation with nQuire.
#' @param type May be equal to fixed, fixed_2, or fixed_3.
#' @param uniform If equal to 1, a uniform mixture is included. If equal to 0, no uniform mixture is included.
#' @param mixtures Defaults to `c("diploid", "triploid", "tetraploid", "hexaploid", "pentaploid")`.
#' @param samplename Name of sample to be included in output.
#' @param trunc List of two values representing the lower and upper bounds for
#'  allele frequency truncation ,\eqn{c_{L}} and \eqn{c_{U}}. If allele frequency
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
#'
#' @return BIC scores and log-likelihood (LL) for the included mixture models.
#'  For BIC, the smallest score is the most likely model.
#'  For LL, the largest score is the most likely model.
#'
#' @importFrom stats dbeta pbeta


bestquack <- function(xm, distribution, type, uniform,
                      mixtures = c("diploid", "triploid", "tetraploid", "hexaploid", "pentaploid"),
                      samplename,
           trunc = c(0.0,0.0),  lowvar = FALSE,
           tau = NA, error = NA){

    cat(paste0(" \t\t <(.)__ <(.)__ <(-)__",
               "\n \t\t  (___/  (___/  (___/  nQuack-in-progress", "\n"))

    # Input data setup
    xm <- as.matrix(xm)
    xi <- (xm[,2]/xm[,1])
    n <- length(xi)
    lnn <- log(n)

     # Input parameters
    set <- list()
    setU <- list()

    if(lowvar == FALSE){

      set[[1]] =  list(avec = c(1.00), mvec = c(0.50), svec = c(0.01));
      set[[2]] =  list(avec = c(0.50, 0.50), mvec = c(0.67, 0.33), svec = c(0.01, 0.01));
      set[[3]] =  list(avec = c(0.33, 0.33, 0.33), mvec = c(0.25, 0.50, 0.75), svec = c(0.01, 0.01,0.01));
      set[[4]] =  list(avec = c(0.25, 0.25, 0.25, 0.25), mvec = c(0.20, 0.40, 0.60, 0.80), svec = c(0.01,0.01, 0.01, 0.01));
      set[[5]] =  list(avec = c(0.20, 0.20, 0.20, 0.20, 0.20), mvec = c(0.17, 0.33, 0.50, 0.67, 0.83), svec = c(0.01, 0.01, 0.01, 0.01, 0.01));

      setU[[1]] =  list(avec = c(0.9, 0.1), mvec = c(0.50), svec = c(0.01));
      setU[[2]] =  list(avec = c(0.45, 0.45, 0.1), mvec = c(0.67, 0.33), svec = c(0.01, 0.01));
      setU[[3]] =  list(avec = c(0.3, 0.3, 0.3, 0.1), mvec = c(0.25, 0.50, 0.75), svec = c(0.01, 0.01, 0.01));
      setU[[4]] =  list(avec = c(0.225, 0.225, 0.225, 0.225, 0.1), mvec = c(0.20, 0.40, 0.60, 0.80), svec = c(0.01,0.01, 0.01, 0.01));
      setU[[5]] =  list(avec = c(0.18, 0.18, 0.18, 0.18, 0.18, 0.1), mvec = c(0.17, 0.33, 0.50, 0.67, 0.83), svec = c(0.01, 0.01, 0.01, 0.01, 0.01));
      vv <- 0.01
      }else if(lowvar == TRUE){

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

    # Mixture subset
    mixturekey <- c("diploid", "triploid", "tetraploid", "hexaploid", "pentaploid")
    if(length(mixtures) != 5){
      subme <- which(mixturekey %in% mixtures)

      set <- set[subme]
      setU <- setU[subme]
    }else{
      set <- set
      setU <- setU
    }


    if(distribution == "normal"){
      # Set up parm list for BIC calculations
      parmlistN <- c()
      parmlistN[[1]] <- c(1,2, 3, 4, 5)
      parmlistN[[2]] <- c(2,4, 6, 8, 10)
      parmlistN[[3]] <- c(1,2, 3, 4, 5)

      parmlistU <- c()
      parmlistU[[1]] <- c(2,3, 4, 5, 6)
      parmlistU[[2]] <- c(3,5, 7, 9, 11)
      parmlistU[[3]] <- c(1,2, 3, 4, 5)

      mixturekey <- c("diploid", "triploid", "tetraploid", "hexaploid", "pentaploid")

      bvec <- c("fixed", "fixed_2", "fixed_3")
      subme2 <- which(bvec %in% type)
      parmlistN <- parmlistN[[subme2]]
      parmlistU <- parmlistU[[subme2]]
      if(length(mixtures) != 5){
        subme <- which(mixturekey %in% mixtures)

        parmlistN <- parmlistN[subme]
        parmlistU <- parmlistU[subme]
      }else{
        parmlistN <- parmlistN
        parmlistU <- parmlistU
      }

      niter = 1000
      epsilon = 0.1

      if(uniform == 0){
        d <- as.data.frame(matrix(nrow = 0, ncol = 5))
        colnames(d) <- c("LL", "type", "mixture", "distribution", "BIC")
        for(i in 1:length(mixtures)){
          obj <-  nQuack::emstepNA(set[[i]], xi, niter, epsilon, trunc,  type = type)
          BICout <- ((-2*obj$loglikelihood) + (lnn*parmlistN[i]))
          d <- rbind(d, data.frame(LL=obj$loglikelihood, type=type, mixture=mixtures[i], distribution = "normal", BIC = BICout))
        }
      }else if(uniform == 1){
        d <- as.data.frame(matrix(nrow = 0, ncol = 5))
        colnames(d) <- c("LL", "type", "mixture", "distribution", "BIC")


        for(i in 1:length(mixtures)){
            obj <-  nQuack::emstepNUA(setU[[i]], xi, niter, epsilon, trunc,  type = type)
            BICout <- ((-2*obj$loglikelihood) + (lnn*parmlistU[i]))
            d <- rbind(d, data.frame(LL=obj$loglikelihood, type=type, mixture=mixtures[i], distribution = "normal-uniform", BIC = BICout))
          }
      }


    }else if(distribution == "beta"){
      # Set up parm list for BIC calculations
      parmlistN <- c()
      parmlistN[[1]] <- c(0, 2, 3, 4, 5)
      parmlistN[[2]] <- c(2, 4, 6, 8, 10) #fixed_2
      parmlistN[[3]] <- c(1, 2, 3, 4, 5)  #fixed_3

      parmlistU <- c()
      parmlistU[[1]] <- c(2, 3, 4, 5, 6)
      parmlistU[[2]] <- c(3, 5, 7, 9, 11)
      parmlistU[[3]] <- c(1, 2, 3, 4, 5)
      mixturekey <- c("diploid", "triploid", "tetraploid", "hexaploid", "pentaploid")

      bvec <- c("fixed", "fixed_2", "fixed_3")
      subme2 <- which(bvec %in% type)
      parmlistN <- parmlistN[[subme2]]
      parmlistU <- parmlistU[[subme2]]
      if(length(mixtures) != 5){
        subme <- which(mixturekey %in% mixtures)

        parmlistN <- parmlistN[subme]
        parmlistU <- parmlistU[subme]
      }else{
        parmlistN <- parmlistN
        parmlistU <- parmlistU
      }

      # EM setup
      niter = 1000
      epsilon = 0.1
      tt <- ifelse(sum(trunc) > 0, TRUE, FALSE)

      if(uniform == 0){
        d <- as.data.frame(matrix(nrow = 0, ncol = 5))
        colnames(d) <- c("LL", "type", "mixture", "distribution", "BIC")

        for(i in 1:length(mixtures)){
            if(mixtures[i] == "diploid" & type == "fixed"){
              ab <- alphabetacalc(set[[1]]$mvec, set[[1]]$svec)
              if(tt == TRUE){
                obj <- list(loglikelihood = c(sum(log(dbeta(xi, ab[1,1], ab[2,1], log = FALSE)/(pbeta(trunc[2], ab[1,1], ab[2,1], 1, 0) - pbeta(trunc[1], ab[1,1], ab[2,1], 1, 0))))))
              }else{
                obj <- list(loglikelihood = c(sum(dbeta(xi, ab[1,1], ab[2,1], log = TRUE))))
              }
            } else{
              obj <-  emstepB(set[[i]], xi, niter, epsilon, trunc,  type = type)
            }

            BICout <- ((-2*obj$loglikelihood) + (lnn*parmlistN[i]))
            d <- rbind(d, data.frame(LL=obj$loglikelihood, type=type,  mixture=mixtures[i], distribution = "beta", BIC = BICout))
        }

      }else if (uniform == 1){

        d <- as.data.frame(matrix(nrow = 0, ncol = 5))
        colnames(d) <- c("LL", "type", "mixture", "distribution", "BIC")

        for(i in 1:length(mixtures)){
            obj <-  emstepBU(setU[[i]], xi, niter, epsilon, trunc,  type = type)
            BICout <- ((-2*obj$loglikelihood) + (lnn*parmlistU[i]))
            d <- rbind(d, data.frame(LL=obj$loglikelihood, type=type, mixture=mixtures[i], distribution = "beta-uniform", BIC = BICout))
          }
      }


    }else if(distribution == "beta-binomial"){
      n <- nrow(xm)
      lnn <- log(n)
      # Set up parm list for BIC calculations
      parmlistN <- c()
      parmlistN[[1]] <- c(0, 2, 3, 4, 5)
      parmlistN[[2]] <- c(2, 4, 6, 8, 10) #fixed_2
      parmlistN[[3]] <- c(1, 2, 3, 4, 5)  #fixed_3

      parmlistU <- c()
      parmlistU[[1]] <- c(2, 3, 4, 5, 6)
      parmlistU[[2]] <- c(3, 5, 7, 9, 11)
      parmlistU[[3]] <- c(1, 2, 3, 4, 5)
      mixturekey <- c("diploid", "triploid", "tetraploid", "hexaploid", "pentaploid")

      bvec <- c("fixed", "fixed_2", "fixed_3")
      subme2 <- which(bvec %in% type)
      parmlistN <- parmlistN[[subme2]]
      parmlistU <- parmlistU[[subme2]]
      if(length(mixtures) != 5){
        subme <- which(mixturekey %in% mixtures)

        parmlistN <- parmlistN[subme]
        parmlistU <- parmlistU[subme]
      }else{
        parmlistN <- parmlistN
        parmlistU <- parmlistU
      }

      # EM setup
      niter = 1000
      epsilon = 0.1
      tt <- ifelse(sum(trunc) > 0, TRUE, FALSE)

      if(uniform == 0){
        d <- as.data.frame(matrix(ncol = 5, nrow = 0))
        colnames(d) <- c("LL", "type", "mixture", "distribution", "BIC")

        for(i in 1:length(mixtures)){
          if(mixtures[i] == "diploid" & type == "fixed"){
              ab <- alphabetacalc(0.5, vv)
              if(tt == TRUE){
                obj <- list(loglikelihood = sum(log(extraDistr::dbbinom(xm[,2], xm[,1], ab[1,1], ab[2,1], log = FALSE)/(extraDistr::pbbinom((xm[,1]*trunc[2]), xm[,1], ab[1,1], ab[2,1]) - extraDistr::pbbinom((xm[,1]*trunc[1]), xm[,1], ab[1,1], ab[2,1])))))
              }else{
                obj <- list(loglikelihood = sum(extraDistr::dbbinom(xm[,2], xm[,1], ab[1,1], ab[2,1], log = TRUE)))
              }
             }else{
              obj <-  emstepBB(set[[i]], xm, niter, epsilon, trunc,  type = type)
            }

            BICout <- ((-2*obj$loglikelihood) + (lnn*parmlistN[i]))
            d <- rbind(d, data.frame(LL=obj$loglikelihood, type= type, mixture=mixtures[i], distribution = "beta-binomial", BIC = BICout))
          }





      }else if (uniform == 1){
        d <- as.data.frame(matrix(ncol = 5, nrow = 0))
        colnames(d) <- c("LL", "type", "mixture", "distribution", "BIC")
        for(i in 1:length(mixtures)){
          obj <-  emstepBBU(setU[[i]], xm, niter, epsilon, trunc,  type = type)
          BICout <- ((-2*obj$loglikelihood) + (lnn*parmlistU[i]))
          d <- rbind(d, data.frame(LL=obj$loglikelihood, type=type, mixture=mixtures[i], distribution = "beta-binomial-uniform", BIC = BICout))

        }
      }

    }else{
      message("distribution may equal normal, beta, or beta-binomial")
    }

  out <- d
  out$sample <- samplename
  return(out)
}
