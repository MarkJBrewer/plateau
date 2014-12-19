#' Fit a plateau envelope via a (non-spatial) GLM and non-linear optimisation
#'
#' \code{fit.glm.env} fits a plateau envelope via a (non-spatial) GLM and
#' non-linear optimisation. This function exists to provide scope for a quick,
#' exploratory analysis and for obtaining starting values for a Bayesian
#' implementation of the plateau envelope fitted via MCMC.
#'
#' @param y The binary response variable (taking values 0 or 1 for absence and
#' presence respectively).
#' @param x.clim The n by p matrix of climate covariates.
#' @param initial.pars.input A vector of length 2p+p+2+p(p-1)/2 containing
#' starting values for each parameter; if missing, the code will suggest
#' starting values.
#' @param random.search A logical variable, which when set to \code{TRUE}
#' will cause the generation of \code{n.iter} random sets of starting values.
#' @param n.iter If \code{random.search} is \code{TRUE}, then \code{n.iter}
#' is number of sets of random starting values to use.
#' @param constrain.beta A p by 2 matrix of logicals in order to indicate
#' which beta parameters should be constrained to not vary too
#' much from its pair - set at most only one of these to \code{TRUE}
#' in each row. A value of \code{TRUE} in column 1 suggests that
#' there are not enough data to estimate the "up" (left-hand)
#' part of the envelope, and in column 2 correspondingly for
#' the "down" (right-hand) part of the envelope. If
#' betas should be constrained.
#' @param slope.limit Scalar putting an upper bound on the envelope slopes;
#' limit is approximately \code{exp(slope.limit)}.
#' @return A list containing the following elements:
#' \describe{
#' \item{\code{par}}{Vector of estimates of the parameters in the order: (2*p) *log*
#'                    of slopes beta_11,beta_12,beta21,beta22,...; (p) apex
#'                    coordinates of the climate variables ax; (1) *log* distance
#'                    beneath the y apex of the "slice" parameter; (1) the apex
#'                    coordinate in the y direction; and if p>1, (p*(p-1)/2)
#'                    pairwise interaction terms.}
#' \item{\code{value}}{The (scalar) value of the "optimised" parameters.}
#' \item{\code{counts}}{The number of \code{glm.env.fn} function calls
#' (\code{NA} for gradient).}
#' \item{\code{convergence}}{Report from \code{optim()} as to convergence
#' status; beware that 0 does not reliably mean "good" convergence.}
#' \item{\code{message}}{Any message from \code{optim()}.}
#' \item{\code{y}}{The response data, returned mainly to aid the plotting
#' functions.}
#' \item{\code{x.clim}}{The climate data, also for the plotting functions.}
#' }
#' @export
fit.glm.env <- function(y,x.clim,initial.pars.input,random.search=FALSE,n.iter=100,
    constrain.beta=FALSE,slope.limit=7){
    x.clim <- as.matrix(x.clim)
    n.x.clim <- ncol(x.clim)
    # Now standardise the climate variables by mapping onto [0,1]
    x.clim.std <- apply(x.clim,2,function(x){
        x.min <- min(x)
        x.max <- max(x)
        return((x-x.min)/(x.max-x.min))
    })
    # If no initial values supplied, we'll work out our own
    if(missing(initial.pars.input)){
        initial.object <- generate.initial.values(y=y,x.clim=x.clim.std,constrain.beta=constrain.beta)
        initial.pars <- initial.object$initial.pars
        constrain.beta <- initial.object$constrain.beta
    }else{ # If initial values supplied, use those
        initial.pars <- initial.pars.input
    }
    # Obtain our first set of estimate parameters from optim()
    fit1 <- optim(par=initial.pars,fn=glm.env.fn,hessian=FALSE,
        control=list(maxit=10000),y=y,x.clim=x.clim.std,constrain.beta=constrain.beta,
        slope.limit=slope.limit)
    current.pars <- fit1$par
    current.best <- fit1$value
    current.best.fit <- fit1
    if(random.search){ # run through n.iter random starts
        for(i in 1:n.iter){
            new.value <- -9e90
            current.best.i <- 9e90
            old.current.best.i <- current.best.i
            while(old.current.best.i>new.value){ # Keep running optim until actual
                                                 # convergence (possibly local)
                old.current.best.i <- current.best.i
                fit1 <- optim(par=current.pars,fn=glm.env.fn,hessian=FALSE,
                    control=list(maxit=10000),y=y,x.clim=x.clim.std,
                    constrain.beta=constrain.beta,slope.limit=slope.limit)
                new.value <- fit1$value
                if(current.best.i>new.value){
                    current.best.i <- new.value
                    current.pars <- fit1$par
                }
            }
            if(current.best.i<current.best){
                current.best <- current.best.i
                current.best.fit <- fit1
                current.pars <- fit1$par
            }
            current.pars <- generate.initial.values(y=y,x.clim=x.clim.std,constrain.beta=constrain.beta,random=TRUE,pars=current.best.fit$par)
            print("Current parameter values:")
            print(current.pars)
            print("Current best deviance:")
            print(current.best)
            if(i==1){
                cat(paste(i," iteration of ",n.iter,".\n",sep=""))
            }else{
                cat(paste(i," iterations of ",n.iter,".\n",sep=""))
            }
        }
    }else{
        new.value <- -9e90
        current.best.fit <- fit1
        old.current.best <- current.best
        while(old.current.best>new.value){ # Keep running optim until actual
                                           # convergence (possibly local
            old.current.best <- current.best
            current.pars <- fit1$par
            fit1 <- optim(par=current.pars,fn=glm.env.fn,hessian=FALSE,
                control=list(maxit=10000),y=y,x.clim=x.clim.std,
                constrain.beta=constrain.beta,slope.limit=slope.limit)
            new.value <- fit1$value
            if(current.best>new.value){
                current.best <- new.value
                current.best.fit <- fit1
            }
        }
    }
    # Add the data, mainly for the plotting functions
    current.best.fit$y <- y
    current.best.fit$x.clim <- x.clim

    return(current.best.fit)
}
