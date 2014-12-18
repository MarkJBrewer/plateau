#' Evaluate the plateau envelope function
#'
#' \code{env.fn} evaluates the plateau envelope for a given set of parameters
#' and a given set of climate covariates. Called by the internal fitting
#' function \code{glm.env.fn} and the map plotting function \code{plot.map}.
#'
#' @param pars The vector of envelope parameters, length 2p+p+2+p(p-1)/2.
#' @param x.clim The n by p matrix of climate covariates.
#' @param slope.limit A scalar putting an upper bound on the envelope slopes;
#' limit is approximately \code{exp(slope.limit)}.
#' @return A list of several named objects: the evaluations of the envelope
#' in a vector \code{x.envelope}; and two objects (\code{beta.mat} and
#' \code{ax.vec}) used by the \code{glm.env.fn} function when penalising
#' the deviance.
#' @export
env.fn <- function(pars,x.clim,slope.limit=25){
    n.x.clim <- ncol(x.clim)
    n.data <- nrow(x.clim)
    beta.top <- 2*n.x.clim
    # Extract the parameters into their respective objects
    beta.mat <- matrix(exp(pars[1:beta.top]),byrow=TRUE,ncol=2)
    ax.vec <- pars[(beta.top+1):(beta.top+n.x.clim)]
    az <- pars[beta.top+n.x.clim+2]
    beta0 <- az-exp(pars[beta.top+n.x.clim+1])
    ax.mat <- matrix(rep(ax.vec,each=n.data),ncol=n.x.clim)
    x.clim.a.centre <- x.clim-ax.mat
    x.mat.ind <- (x.clim.a.centre > 0)+1
    beta.mat.prod <- array(dim=c(n.data,n.x.clim))
    for(i in 1:n.x.clim){
        beta.mat.prod[,i] <- beta.mat[i,x.mat.ind[,i]]
    }
    # the contributions from the main effects of climate
    x.sum.main <- rowSums(beta.mat.prod*x.clim.a.centre*x.clim.a.centre)
    if(n.x.clim!=1){ # no interactions if p==1
        gamma.mat <- array(NA,dim=c(n.x.clim,n.x.clim))
        gamma.mat[t(upper.tri(gamma.mat))] <- pars[(beta.top+n.x.clim+3):length(pars)]
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
        n.interactions <- nrow(gamma.rc)
        for(i in 1:n.interactions){
            gamma.mat[gamma.rc.row[i],gamma.rc.col[i]] <- gamma.min[i]+gamma.mat[gamma.rc.row[i],gamma.rc.col[i]]*(gamma.max[i]-gamma.min[i])
        }
        int.mat <- x.clim.a.centre[,gamma.rc]
        int.mat.1 <- as.matrix(int.mat[,1:n.interactions])
        int.mat.2 <- as.matrix(int.mat[,(n.interactions+1):(2*n.interactions)])
        # the contributions from the interactions of climate variables
        x.sum.int <- rowSums(int.mat.1*int.mat.2*matrix(rep(gamma.mat[gamma.rc],n.data),ncol=n.interactions,byrow=TRUE))
    }else{
        x.sum.int <- 0
    }
    if(any(!is.finite(x.sum.int))){
        return(list(x.envelope=0))
    }
    if(any(x.sum.main+x.sum.int < 0)){
        if(any(x.sum.main+x.sum.int < -1e-8)){
            stop("Negative square root in env.fn.")
        }else{
            x.sum.int[x.sum.main+x.sum.int < 0] <- x.sum.main[x.sum.main+x.sum.int < 0]+1e16
        }
    }
    env.fn.object <- list(x.envelope=pmin(az-(sqrt(x.sum.main+x.sum.int)),beta0))
    env.fn.object$beta.mat <- beta.mat
    env.fn.object$ax.vec <- ax.vec
    return(env.fn.object)
}
