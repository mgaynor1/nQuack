#' Simulate Allele Counts for Single Individual - Simple Approach
#'
#' @description This function is used to simulate coverage of each allele at
#'   biallelic heterozygous sites assuming a binomial distribution.
#'
#' @param mvec Vector of means.
#' @param s.size Number of biallelic sites to generate. Defaults to 50000. Warning,
#'    the number of sites generated will not be the number of sites returned
#'    due to filtering steps.
#' @param cover Coverage of sites.
#' @param sampled Default as TRUE. Will randomly sample allele A or allele B, then return a data
#'    frame with total coverage and coverage of a randomly sampled allele will be returned.
#'
#' @return If sampled = FALSE, a data frame with total coverage, coverage of allele A,
#'    and coverage of allele B will be returned.  If sampled = TRUE, a data frame with total coverage
#'    and coverage of a randomly sampled allele will be returned.




sim.ind.simple <- function(mvec, cover = 100, s.size = 50000, sampled = TRUE){

  all.sites <- data.frame(matrix(nrow = s.size, ncol = 3))
  all.sites[, 1] <- cover
  all.sites[, 2] <- rbinom(n = s.size, size = cover, prob = mvec)

  for(i in 1:s.size) {
    all.sites[i, 3] <- all.sites[i, 1] - all.sites[i, 2]
  }

  # Remove homozygotes
  all.sites <- all.sites[which(all.sites[,1] >= 2 & all.sites[,2] != 0 & all.sites[,3] != 0),]

  # Reduce data duplication
  if(sampled == TRUE){
    all.sites.sample <- data.frame(matrix(nrow = nrow(all.sites), ncol = 2))
    for(i in 1:nrow(all.sites)){
      all.sites.sample[i, 1] <- all.sites[i, 1]
      all.sites.sample[i, 2] <- sample(x = all.sites[i, 2:3], size = 1, prob = c(0.5, 0.5))
    }
    return(all.sites.sample)
  } else{
    colnames(all.sites) <- c("Total.Coverage", "Allele.A", "Allele.B")
    return(all.sites)
  }
}
