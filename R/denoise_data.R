#' Denoise Data
#'
#' @description Here we filter allele frequencies with a normal + uniform
#' mixture model. nQuack utilizes the scaled probability of each data point belonging
#' to each mixture model, which is inferred in the expectation maximization algorithm.
#' We remove allele frequencies where probability of belonging to uniform mixture is higher than their
#' probability of belonging to any other mixture. We also implement nQuire's denoise method here, which
#' utilizes the inferred alpha parameter and a histogram of base frequencies to filter the data.
#'
#' @param xm  Matrix with total coverage and coverage for a randomly sampled allele.
#' @param plot Default to TRUE. The plots do not share the same y-axis, so careful interpretation is key.
#'  Warning, if nothing is removed, the plot of removed data will be missing.
#' @param filter Indicates which method to remove data based upon. Options: 'both', 'nquire', or 'nquack'. nQuack
#'  utilizes the scaled probability of each data point belonging to each mixture model, removing sites where the probability of belonging to uniform mixture is higher than their
#' probability of belonging to any other mixture. nQuire utilizes the inferred alpha parameter and a histogram of base frequencies to filter the data.
#'
#' @return Numeric matrix with total coverage and coverage for a randomly sampled allele.
#'
#' @importFrom graphics hist par
#' @importFrom stats runif
#'

denoise_data <- function(xm, plot = TRUE, filter = "both"){

  # Setup data
  xm <- as.matrix(xm)
  xi <- xm[,2]/xm[,1]

  # Free Model Run
  set <-  list(avec = c(0.11, 0.22, 0.33, 0.22, 0.11, 0.01), mvec = c(0.20, 0.33, 0.50, 0.67, 0.80), svec = c(0.01, 0.01, 0.01, 0.01, 0.01));
  nout <- emstepNUA(set, xi, 1000, 0.1, c(0,0), "free" )
  pir <- nout$pir

  # nQuire's way:
  if(min(xi) > 0){
    os = round(min(xi)*100)
    size = 100-(2*os)
  }else if(min(xi) == 0 | is.na(min(xi))){
    os = 1
    size = 100
  }

  min = ((length(xi) * nout$parm.list$avec[6])/size)
  details <- hist(xi, breaks = size, plot = FALSE)
  y <- details$counts
  for(i in 1:size){
   f = ((y[i] - min)/y[i])

    if(f < 0.01 | is.na(f)){
      y[i] = 0.01
    }else{
      y[i] = f
    }
  }
  rl <- c()
  for(i in 1:length(xi)){
    val = round(xi[i]*100)
    ys = val-os
    if(ys == 0 | is.na(ys)){ys <- 1}else if(ys > size){ ys <- size-1}
    if((val >= os) & (val <= (100-os)) & (runif(1) < y[ys]) ){
      rl[i] <- 1
    } else{
      rl[i] <- 0
    }
  }

  # nQuack's way
  removelist <- c()
  for(i in 1:length(xi)){
  if(which(pir[i,] == max(pir[i,])) == 6){
      removelist[i] <- 0
    }else{
      removelist[i] <- 1
    }
  }

  if(filter == "both"){
    if(plot == TRUE){
      xikeep <- xi[which(removelist == 1)]
      xiremove <- xi[which(removelist == 0)]
      xiQkeep <- xi[which(rl == 1)]
      xiQremove <- xi[which(rl == 0)]
      bothkeep <- xi[-which(removelist == 0 | rl == 0)]
      par(mfrow = c(3, 3))

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

      hist(xi,xlim = c(0, 1), main = "Original Data", xlab = "Allele Frequency", breaks = size )
      if(length(xiQkeep) == 0){
        plot(0, type="n", axes=FALSE, xlab="", ylab="")
      }else{
        hist(xiQkeep, xlim = c(0, 1), main = "Data Retained - nQuire", xlab = "Allele Frequency", breaks = size)
      }
      if(length(xiQremove) == 0){
        plot(0, type="n", axes=FALSE, xlab="", ylab="")
      }else{
        hist(xiQremove,  xlim = c(0, 1), main = "Data Removed - nQuire",  xlab = "Allele Frequency", breaks = size)
      }

      plot(0, type="n", axes=FALSE, xlab="", ylab="")

      if(length(bothkeep) == 0){
          plot(0, type="n", axes=FALSE, xlab="", ylab="")
        }else{
          hist(bothkeep,xlim = c(0, 1), main = "Combined - Data Returned",  xlab = "Allele Frequency", breaks = size )
        }
    }

    xmout  <- xm[-which(removelist == 0 | rl == 0), ]
    return(xmout)
  }else if(filter == "nquire"){
    if(plot == TRUE){

      xiQkeep <- xi[which(rl == 1)]
      xiQremove <- xi[which(rl == 0)]
      par(mfrow = c(1, 3))

      hist(xi,xlim = c(0, 1), main = "Original Data", xlab = "Allele Frequency", breaks = size )

      if(length(xiQkeep) == 0){
        plot(0, type="n", axes=FALSE, xlab="", ylab="")
      }else{
        hist(xiQkeep, xlim = c(0, 1), main = "Data Retained - nQuire", xlab = "Allele Frequency", breaks = size)
      }
      if(length(xiQremove) == 0){
        plot(0, type="n", axes=FALSE, xlab="", ylab="")
      }else{
        hist(xiQremove,  xlim = c(0, 1), main = "Data Removed - nQuire",  xlab = "Allele Frequency", breaks = size)
      }
    }
    xmout <- xm[-which(rl == 0),]
    return(xmout)
  } else if(filter == "nquack"){
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
      }else{
        hist(xiremove,  xlim = c(0, 1), main = "Data Removed - nQuack",  xlab = "Allele Frequency", breaks = size)
      }
    }

    xmout <- xm[-which(removelist == 0),]
    return(xmout)
  }


}
