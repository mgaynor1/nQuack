#' Model Selection - Based on BIC or Log-Likelihood
#'
#' @description This function is for model interpretation.
#'
#' @param model_out Data frame containing, at minimum, columns labeled LL, type, mixture, distribution, and BIC.
#' @param summary_statistic May be equal to BIC or dLL.
#' @param mixtures Defaults to `c("diploid", "triploid", "tetraploid", "hexaploid", "pentaploid")`.
#'
#' @return Returns data frame with the most likely model for each set of mixtures. Includes the best and second best mixtures, as well as the difference between the two.
#'
#'
quackit <- function(model_out, summary_statistic = "BIC",
                    mixtures = c("diploid", "triploid", "tetraploid", "hexaploid", "pentaploid")){
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
 } else if(summary_statistic == "dLL"){
   colnames( out2go) <- c("dLLdif", "winnerdLL", "seconddLL", "Distribution", "Type", "sample")
 }
 return(out2go)

}
