#' Fit a plateau envelope via a (non-spatial) GLM and non-linear optimisation
#'
#' \code{fit.glm.env} fits a plateau envelope via a (non-spatial) GLM and
#' non-linear optimisation. This function exists to provide scope for a quick,
#' exploratory analysis and for obtaining starting values for a Bayesian
#' implementation of the plateau envelope fitted via MCMC.
#'
#' @param data The data frame (with n rows) containing all the variables for analysis.
#' @param y A string denoting the binary response variable (taking values 0 or 1 for absence and
#' presence respectively); must correspond to a column name in the data frame
#' specied at \code{data}.
#' @param x.clim A vector (length p) of strings denoting which columns in the supplied data
#' frame correspond to the climate covariates; must correspond to column names in the data frame
#' specied at \code{data}.
#' @param x.nonclim A vector (length p2) of strings denoting which columns in the supplied data
#' frame correspond to the  non-climate covariates; must correspond to column names in the data frame
#' specied at \code{data}.
#' @param x.factor A vector (length p3) of strings denoting which columns in the supplied data
#' frame correspond to the non-climate factors; must correspond to column names in the data frame
#' specied at \code{data}.
#' @param initial.pars.input A vector of length 2p+p+2+p(p-1)/2 containing
#' starting values for each parameter; if missing, the code will suggest
#' starting values. If \code{x.factor} is non-null, there will be a total of
#' the number of total levels across all the factors minus p3 extra parameters.
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
#' @param silent Logical flag denoting whether the function runs silently or
#' not. Default is \code{TRUE}.
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
fit.glm.env <- function(data,y,x.clim,x.nonclim=NULL,x.factor=NULL,
    initial.pars.input,random.search=FALSE,n.iter=100,constrain.beta=FALSE,
    slope.limit=7,silent=TRUE){
    contr.env <- function(n){
        return(diag(n)[,1:(n-1)]-diag(n)[,2:n])
    }
    y.name <- y
    y <- data[,y]
    n.x.clim <- length(x.clim)
    x.clim <- as.matrix(data[,x.clim])
    x.nonclim.names <- x.nonclim
    x.factor.names <- x.factor
    if(!is.null(x.nonclim)){
        if(!is.null(x.factor)){
            x.nonclim <- as.matrix(data[,x.nonclim])
            x.factor <- data[,x.factor]
            x.nonclim.formula <- as.formula(paste(y.name,"~offset(x.envelope)-1+",
                paste(x.nonclim.names,collapse="+"),
                paste(x.factor.names,collapse="+"),sep=""))
        }else{
            x.nonclim <- as.matrix(data[,x.nonclim])
            x.nonclim.formula <- as.formula(paste(y.name,"~offset(x.envelope)-1+",
                paste(x.nonclim.names,collapse="+"),sep=""))
        }
    }else{
        if(!is.null(x.factor)){
            x.factor <- data[,x.factor]
            x.nonclim.formula <- as.formula(paste(y.name,"~offset(x.envelope)-1+",
                paste(x.factor.names,collapse="+"),sep=""))
        }else{
            x.nonclim.formula <- NULL
        }
    }
    # Now standardise the climate variables by mapping onto [0,1]
    x.clim.std <- apply(x.clim,2,function(x){
        x.min <- min(x)
        x.max <- max(x)
        return((x-x.min)/(x.max-x.min))
    })
    # If no initial values supplied, we'll work out our own
    if(missing(initial.pars.input)){
        initial.object <- generate.initial.values(y=y,x.clim=x.clim.std,
            x.nonclim=x.nonclim,x.factor=x.factor,
            constrain.beta=constrain.beta)
        initial.pars <- initial.object$initial.pars
        constrain.beta <- initial.object$constrain.beta
    }else{ # If initial values supplied, use those
        initial.pars <- initial.pars.input
    }
    # Obtain our first set of estimate parameters from optim()
    fit1 <- optim(par=initial.pars,fn=glm.env.fn,hessian=FALSE,
        control=list(maxit=10000),data=data,y=y.name,x.clim=x.clim.std,x.nonclim=x.nonclim.names,
        x.nonclim.formula=x.nonclim.formula,x.factor=x.factor.names,
        constrain.beta=constrain.beta,slope.limit=slope.limit)
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
                    control=list(maxit=10000),data=data,y=y.name,x.clim=x.clim.std,x.nonclim=x.nonclim.names,
                    x.nonclim.formula=x.nonclim.formula,x.factor=x.factor.names,
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
            current.pars <- generate.initial.values(y=y,x.clim=x.clim.std,
                x.nonclim=x.nonclim,x.factor=x.factor,
                constrain.beta=constrain.beta,random=TRUE,
                pars=current.best.fit$par)
            if(!silent){
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
                control=list(maxit=10000),data=data,y=y.name,x.clim=x.clim.std,x.nonclim=x.nonclim.names,
                x.nonclim.formula=x.nonclim.formula,x.factor=x.factor.names,
                constrain.beta=constrain.beta,slope.limit=slope.limit)
            new.value <- fit1$value
            if(current.best>new.value){
                current.best <- new.value
                current.best.fit <- fit1
            }
        }
    }
    if(!is.null(x.nonclim)){
      if(!is.null(x.factor)){
        x.envelope <- env.fn(pars=current.best.fit$par,x.clim=x.clim.std,
            slope.limit=slope.limit)$x.envelope
        tempdata <- data[,c(y.name,x.nonclim.names,x.factor.names)]
        tempdata$x.envelope <- x.envelope
        if(length(x.factor.names)==1){
          temp.nlevels <- nlevels(tempdata[,x.factor.names])
          tempdata[,x.factor.names] <- (contr.env(temp.nlevels))[tempdata[,x.factor.names],]
        }else{
          for(i in 1:length(x.factor.names)){
            temp.nlevels <- nlevels(tempdata[,x.factor.names[i]])
            tempdata[,x.factor.names[i]] <- (contr.env(temp.nlevels))[tempdata[,x.factor.names[i]],]
          }
        }
        current.best.fit$nonclim.glm <- glm(x.nonclim.formula,data=tempdata,
            family=binomial)
      }else{
        x.envelope <- env.fn(pars=current.best.fit$par,x.clim=x.clim.std,
                             slope.limit=slope.limit)$x.envelope
        tempdata <- data[,c(y.name,x.nonclim.names)]
        tempdata$x.envelope <- x.envelope
        current.best.fit$nonclim.glm <- glm(x.nonclim.formula,data=tempdata,
              family=binomial)
      }
    }else{
      if(!is.null(x.factor)){
        x.envelope <- env.fn(pars=current.best.fit$par,x.clim=x.clim.std,
                             slope.limit=slope.limit)$x.envelope
        tempdata <- data[,c(y.name,x.factor.names)]
        tempdata$x.envelope <- x.envelope
        if(length(x.factor.names)==1){
          temp.nlevels <- nlevels(tempdata[,x.factor.names])
          tempdata[,x.factor.names] <- (contr.env(temp.nlevels))[tempdata[,x.factor.names],]
        }else{
          for(i in 1:length(x.factor.names)){
            temp.nlevels <- nlevels(tempdata[,x.factor.names[i]])
            tempdata[,x.factor.names[i]] <- (contr.env(temp.nlevels))[tempdata[,x.factor.names[i]],]
          }
        }
        current.best.fit$nonclim.glm <- glm(x.nonclim.formula,data=tempdata,
                                            family=binomial)        
      }else{
        current.best.fit$nonclim.glm <- NULL
      }
    }
    # Add the data, mainly for the plotting functions
    current.best.fit$y <- y
    current.best.fit$x.clim <- x.clim
    current.best.fit$x.nonclim <- x.nonclim
    current.best.fit$x.factor <- x.factor
    
    return(current.best.fit)
}
