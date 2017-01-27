#' Predict values for the plateau envelope models produced by \code{fit.glm.env} and
#' \code{fit.bugs.env}.
#'
#' \code{envelope.predict} predicts values for the plateau envelope models produced by the functions
#' \code{fit.glm.env} and \code{fit.bugs.env}.
#'
#' @param envelope.fit The fitted object from either \code{fit.glm.env} or
#' \code{fit.bugs.env}.
#' 
#' @param newdata A matrix or data frame of new data for prediction.
#' 
#' @return A set of predictions from the supplied fitted model evaluated at the covariate values in
#' object \code{newdata}.
#' @seealso \code{link{map.plot}}, \code{link{fit.bugs.env}}
#' @export
envelope.predict <- function(envelope.fit,newdata){
  x.clim <- as.matrix(envelope.fit$x.clim)
  x.clim.names <- colnames(x.clim)
  n.x.clim <- ncol(x.clim)
  
  # modified standized part  
  x.min <- apply(x.clim,2,min)
  x.max <- apply(x.clim,2,max)
  x.diff <- x.max-x.min
  
  scaleTo <- function(anydata){
    for (ss in 1: ncol(x.clim) ) {
      anydata[,ss] <- (anydata[,ss]-x.min[ss])/x.diff[ss]
    }
    return(anydata)
  }
  
  scaleBack <- function(anydata){
    for (ss in 1: ncol(x.clim) ) {
      anydata[,ss] <- anydata[,ss] * x.diff[ss]+ x.min[ss]
    }
    return(anydata)
  }
  
  newdata.std  <- scaleTo(newdata)
  x.clim.std <-   scaleTo(x.clim)
  
  beta.top <- 2*n.x.clim
  beta.mat <- matrix(exp(envelope.fit$par[1:beta.top]),byrow=TRUE,ncol=2)
  ax <- envelope.fit$par[(beta.top+1):(beta.top+n.x.clim)]
  az <- envelope.fit$par[beta.top+n.x.clim+2]
  beta0 <- az-exp(envelope.fit$par[beta.top+n.x.clim+1])
  if(!n.x.clim==1){ # if more than 1 variable, then consider interaction
    gamma.mat <- array(NA,dim=c(n.x.clim,n.x.clim))
    gamma.mat[t(upper.tri(gamma.mat))] <- envelope.fit$par[(beta.top+n.x.clim+3):length(envelope.fit$par)]
    gamma.mat <- t(gamma.mat)
    gamma.mat[gamma.mat>20] <- 20 # To stop Inf in the line below (gives same result - 1)
    gamma.mat <- exp(gamma.mat)/(1+exp(gamma.mat))
    # Combinations of row, column with values
    gamma.rc <- which(!is.na(gamma.mat),arr.ind=TRUE)
    gamma.rc.row <- gamma.rc[,"row"]
    gamma.rc.col <- gamma.rc[,"col"]
    gamma.min.1 <- beta.mat[gamma.rc.row,1]*beta.mat[gamma.rc.col,1]
    gamma.min.2 <- beta.mat[gamma.rc.row,2]*beta.mat[gamma.rc.col,2]
    gamma.min <- -(2/(n.x.clim-1))*sqrt(pmin(gamma.min.1,gamma.min.2))
    gamma.max.1 <- beta.mat[gamma.rc.row,1]*beta.mat[gamma.rc.col,2]
    gamma.max.2 <- beta.mat[gamma.rc.row,2]*beta.mat[gamma.rc.col,1]
    gamma.max <- (2/(n.x.clim-1))*sqrt(pmin(gamma.max.1,gamma.max.2))
    for(i in 1:(n.x.clim*(n.x.clim-1)/2)){
      gamma.mat[gamma.rc.row[i],gamma.rc.col[i]] <- gamma.min[i]+gamma.mat[gamma.rc.row[i],gamma.rc.col[i]]*(gamma.max[i]-gamma.min[i])
    }
  }
  
  # if(missing(x.labels)){
  #   x.labels <- paste("Variable",1:n.x.clim)
  # }
  # x.plot.seq <- seq(x.plot.lims[1],x.plot.lims[2],length=len)
  # x.plot <- array(NA,dim=c(len^n.x.clim,n.x.clim))
  # for(i in 1:n.x.clim){
  #   x.plot[,i] <- rep(x.plot.seq,times=len^(i-1),each=len^(n.x.clim-i))
  # }
  
  x.plot <- newdata.std
  
  ax.mat <- matrix(rep(ax,each=nrow(x.plot)),ncol=n.x.clim)
  x.plot.a.centre <- x.plot-ax.mat
  x.mat.ind <- (x.plot.a.centre > 0)+1  # decide left or right of the peak
  beta.mat.prod <- array(dim=c(nrow(x.plot),n.x.clim))
  for(i in 1:n.x.clim){
    beta.mat.prod[,i] <- beta.mat[i,x.mat.ind[,i]]
  }
  #colnames(x.plot) <- x.clim.names
  x.sum.main <- rowSums(beta.mat.prod*x.plot.a.centre^2)
  if(n.x.clim!=1){
    # Combinations of row, column with values
    gamma.rc <- which(!is.na(gamma.mat),arr.ind=TRUE)
    n.interactions <- nrow(gamma.rc)
    int.mat <- x.plot.a.centre[,gamma.rc]
    int.mat.1 <- as.matrix(int.mat[,1:n.interactions])
    int.mat.2 <- as.matrix(int.mat[,(n.interactions+1):(2*n.interactions)])
    x.sum.int <- rowSums(int.mat.1*int.mat.2*matrix(rep(gamma.mat[gamma.rc],nrow(x.plot)),ncol=n.interactions,byrow=TRUE))
  }else{
    x.sum.int <- 0
  }
  x.envelope <- pmin(az-(sqrt(x.sum.main+x.sum.int)),beta0)
  logityrange <- 1/(1+exp(-x.envelope))
  
  finalprediction <- logityrange
  
  return(finalprediction)
}
