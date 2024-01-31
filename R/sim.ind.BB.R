#' Simulate Allele Counts for Single Individual - Beta-Binomial Distribution
#'
#' @description This function is used to simulate coverage of each allele at
#'   biallelic heterozygous sites assuming a beta binomial distribution. Here we
#'   sample sequence depth from a truncated poisson distribution between a set minimum, maximum,
#'   and lambda. Only heterozygous sites are returned.
#'   Based on input variables, the sites may be filtered based on the total coverage (`filter.coverage`), allele
#'   sequencing coverage (`filter.error`), or allele frequency (`filter.freq`).
#'
#' @param mvec Vector of mean values of allele frequency.
#' @param avec Vector of alpha values representing the proportion expected of each mean.
#' @param svec Vector of variance values.
#' @param error Sequencing error rate. Default as 0.001, or very low error.
#' @param s.size Number of biallelic sites to generate. Defaults to 50000. Warning,
#'    the number of sites generated will not be the number of sites returned
#'    due to filtering steps.
#' @param max.coverage Maximum sequencing depth per site. Defaults to 20.
#' @param min.coverage Minimum sequencing depth per site. Defaults to 2.
#' @param lambda Set lambda for the truncated poisson distrubtion. Defaults to 11.
#' @param filter.coverage Default as TRUE. Filters to only retain sites where total sequencing depth
#'    is greater than the provided minimum coverage and less than the max quantile depth (set with the max.depth.quantile.prob).
#' @param max.depth.quantile.prob Maximum depth quantile probability. Defaults to 0.9.
#' @param filter.error Default as TRUE. Filter to only retain sites where allele coverage
#'    is greater than the sequencing error rate times the total coverage, but
#'    less than one minus the sequencing error rate times the total coverage.
#' @param filter.freq Default as FALSE. When set to true, sites are filtered based on provided `trunc`.
#' @param trunc List of two values representing the lower and upper bounds,\eqn{c_{L}} and \eqn{c_{U}}.
#'    Defaults as c(0,0) to represent no truncation.
#' @param sampled Default as TRUE. Will randomly sample allele A or allele B, then return a data
#'    frame with total coverage and coverage of a randomly sampled allele will be returned.
#'
#' @return If sampled = FALSE, a data frame with total coverage, coverage of allele A,
#'    and coverage of allele B will be returned. If sampled = TRUE, a data frame with total coverage
#'    and coverage of a randomly sampled allele will be returned.
#'
#'
#' @importFrom truncdist rtrunc
#' @importFrom stats rbeta rbinom quantile



sim.ind.BB <- function(mvec, avec, svec,  error = 0.001,
                           s.size = 50000, lambda = 11,
                           max.coverage = 20, min.coverage = 2,
                           filter.coverage = TRUE, max.depth.quantile.prob = 0.9,
                           filter.error = TRUE,
                           filter.freq = FALSE, trunc = c(0, 0), sampled = TRUE){

  alpha.beta <- alphabetacalcvec(mvec, svec)
  smax <- (length(avec))
  all.sites <- data.frame(matrix(nrow = s.size, ncol = 3))
  coverage.dist <- truncdist::rtrunc(s.size, spec="pois", a = min.coverage, b=max.coverage, lambda=lambda)

  for(i in 1:s.size) {
    ## Which heterozygote
    which.type <- sample(1:smax, size = 1, replace = TRUE, prob = avec)
    ## Sequencing Depth
    cover <- sample(x= coverage.dist, size=1)
    ## Fill out the matrix
    P <- rbeta(n = 1, alpha.beta[which.type, 1], alpha.beta[which.type, 2])
    all.sites[i, 2] <-  rbinom(n = 1, size = cover, prob = P )
    all.sites[i, 1] <- cover
    all.sites[i, 3] <- all.sites[i, 1] - all.sites[i, 2]
  }


  # Filter simulated data
  ## Remove homozygote sites
  all.sites <- all.sites[which(all.sites[,1] >= 2 & all.sites[,2] != 0 & all.sites[,3] != 0),]

  if(filter.coverage == TRUE){
    max <- ceiling(quantile(all.sites[,1], max.depth.quantile.prob))
    all.sites <- all.sites[which(all.sites[,1] > min.coverage & (all.sites[,1] < max)), ]
  }

  if(filter.error == TRUE){
    ## If Coverage(A) > e*Coverage(A+B)
    ## If Covergae(A) < (1-e)*Coverage(A+B)
    all.sites <- all.sites[which(all.sites[,2] > (error*all.sites[,1]) & (all.sites[,2] < ((1-error)*all.sites[,1]) )), ]
    all.sites <- all.sites[which(all.sites[,3] > (error*all.sites[,1]) & (all.sites[,3] < ((1-error)*all.sites[,1]) )), ]
  }

  if(filter.freq == TRUE){
    if(sum(trunc) > 0){
      for(i in 1:nrow(all.sites)){
        all.sites[i,4] <-  (all.sites[i,2]/all.sites[i,1])
        all.sites[i,5] <- (all.sites[i,3]/all.sites[i,1])
      }

      all.sites <-  all.sites[which(all.sites[,4] >= trunc[1] & all.sites[,4] <= trunc[2] & all.sites[,5] >= trunc[1] & all.sites[,5] <= trunc[2]), ]

    } else{
      print("Set trunc")
    }
  }

  all.sites <- all.sites[,1:3]

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



