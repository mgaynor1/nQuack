#' Model Selection - Based on BIC or Log-Likelihood
#'
#' @description This function is for model interpretation.
#'
#' @param model_out Data frame containing, at minimum, columns labeled LL, type, mixture, distribution, and BIC.
#' @param summary_statistic May be equal to BIC or LL.
#' @param mixtures Defaults to `c("diploid", "triploid", "tetraploid", "hexaploid", "pentaploid")`.
#'
#' @examples
#' out <- quackNormal(xm[1:100,], samplename = "sample1", cores = 1)
#' goose <- quackit(out)
#' @returns Returns data frame with the most likely model for each set of mixtures.
#' Includes the best and second best mixtures, as well as the difference between the two.
#' We only use BIC or LL to compare within each distribution and type.
#' To identify the most accurate model, you will need to compare accuracy across distributions
#'  and types using a set of known samples. The distributions include
#'  Normal, Beta, and Beta-Binomial - each with and without a uniform mixture.
#'  The type indicates which parameters are estimated for the mixtures:
#'  all parameters (`type = 'free'`, only used to calculate delta log-likelihood),
#'  only alpha  (`type = 'fixed'`), only alpha and variance (`type = 'fixed_2'`),
#'   and only variance (`type ='fixed_3`) to be estimated for each mixture.
#'
#' @export
quackit <- function(model_out, summary_statistic = "BIC",
                    mixtures = c("diploid", "triploid", "tetraploid", "hexaploid", "pentaploid")){
  if(summary_statistic == "dLL"){
    summary_statistic <- "LL"
  }else{
    summary_statistic <- summary_statistic
  }
  # Set up
  if(length(mixtures) != 5){
    model_out <- model_out[which(model_out$mixture %in% mixtures), ]
  }else{
    model_out <- model_out
  }
  if(length(which(model_out$mixture== "free")) != 0){
    model_out <-  model_out[which(model_out$mixture != "free"), ]
  } else{
    model_out <- model_out
  }

  # Subsample by distribution type and fixed type
 results <- c()

 dl <- unique(model_out$distribution)
 for(j in 1:length(dl)){
   dset <- model_out[which(model_out$distribution == dl[j]), ]
   fl <- unique(dset$type)
   winners <- as.data.frame(matrix(nrow = length(fl), ncol = 6 ))
   for(k in 1:length(fl)){
     dfsub <- dset[which(dset$type == fl[k]), ]
     winnerBIC <- (dfsub[which(dfsub[[summary_statistic]] == min(dfsub[[summary_statistic]] )), ])[[summary_statistic]]
     subsub <- dfsub[-c(which(dfsub[[summary_statistic]]  == min(dfsub[[summary_statistic]] ))), ]
     secondBIC <- min(subsub[[summary_statistic]] )
     BICdif <- (winnerBIC - secondBIC)
     winners[k, 1] <- (winnerBIC - secondBIC)
     winners[k, 2] <- (dfsub[which(dfsub[[summary_statistic]]  == winnerBIC), ])$mixture
     winners[k, 3] <- (dfsub[which(dfsub[[summary_statistic]]  == secondBIC), ])$mixture
     winners[k, 4] <- unique(dfsub$distribution)
     winners[k, 5] <- unique(dfsub$type)
     winners[k, 6] <- unique(dfsub$sample)

   }
   results[[j]] <- winners
 }

 out2go <-  do.call(rbind, results)

 if(summary_statistic == "BIC"){
   colnames( out2go) <- c("BICdif", "winnerBIC", "secondBIC", "Distribution", "Type", "sample")
 } else if(summary_statistic == "LL"){
   colnames( out2go) <- c("LLdif", "winnerLL", "secondLL", "Distribution", "Type", "sample")
 }
 return(out2go)

}
