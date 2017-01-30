#' Fit a spatial, Bayesian, plateau envelope model via MCMC in WinBUGS
#'
#' \code{fit.bugs.env} fits a spatial GLM with plateau envelope on the climate
#' covariates via MCMC in WinBUGS.
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
#' @param not.spatial Logical scalar specifying whether the model should include
#' a spatial intrinsic CAR component. Defaults to FALSE, i.e. there \emph{is} a
#' spatial component.
#' @param car.sigma Standard deviation of the WinBUGS \code{car.normal}
#' process; needs to be a constant in the case of a binary response.
#' @param num A vector (length n) of the numbers of neighbours of each cell, in
#' order (but see later comment on cliques).
#' @param adj A vector of adjacencies, length 2*number of neighbourhood pairs.
#' @param u Vector of spatial random effects to be supplied as data, of length
#' n; mostly \code{NA}s, meaning the u's should be estimated in
#' \code{car.normal}, but you may wish to fix some values (e.g. to 0
#' for disconnected cells). If u is missing and there is only one "clique" (i.e.
#' connected neighbourhood) then u will be created as a vector of \code{NA}s of
#' the right length. u should be specified otherwise; this is deliberate.
#'
#' The next four inputs are used to set identifiability constraints on the
#' modelling.
#' @param prior.ax A list of up to two p-vector objects, \code{mean} and
#' \code{var}, containing the means and variances of (informative) priors
#' for the apex x co-ordinates ax; if either \code{mean} or \code{var} are
#' omitted, the means and variances are set to 0.5 and 1 respectively.
#' @param prior.beta A list of two p by 2 matrix objects, \code{mean} and
#' \code{var}, containing the means and variances of (informative) priors
#' for the slopes beta; if omitted, the means and variances are
#' set to 100 and 10000 respectively.
#' @param prior.beta0.difference A list of two scalar objects, \code{mean} and
#' \code{var}, containing the mean and variance for the log-difference
#' between the apex az and plateau level beta0; if omitted,
#' the mean and variance are set to 1 and 10 respectively.
#' @param constrain.beta A p by 2 matrix of logicals in order to indicate
#' which beta parameters should be constrained to not vary too
#' much from its pair - set at most only one of these to \code{TRUE}
#' in each row. A value of \code{TRUE} in column 1 suggests that
#' there are not enough data to estimate the "up" (left-hand)
#' part of the envelope, and in column 2 correspondingly for
#' the "down" (right-hand) part of the envelope. If
#' \code{initial.pars.input} is not set, the code will work out which
#' betas should be constrained.
#' @param initial.pars.input Vector of length 2p+p+2+p(p-1)/2 containing
#' starting values for each parameter; if missing, the code
#' works out its own starting value. If \code{x.factor} is non-null, there will be a total of
#' the number of total levels across all the factors minus p3 extra parameters.
#' @param informative.priors List of logical scalars for which informative
#' priors should be used (from: beta, beta0, ax);
#' this option uses the function \code{generate.initial.values} to
#' generate the means, and is an alternative to specifying the
#' priors yourself using \code{prior.ax} above etc. Default has all
#' options set to \code{FALSE}.
#' @param burnin Scalar specifying the number of "burn-in" iterations to be
#' discarded (default 5000).
#' @param post.burnin Scalar specifying the number of subsequent iterations to
#' be retained (default 1000).
#' @param chains Scalar, number of parallel chains to run (default 2).
#' @param thin Scalar, the thinning to apply to the MCMC iterations (default 1).
#' @param working.directory String containing the location of the WinBUGS code
#' file; default NULL, in which case a temporary folder is used.
#' @param bugs.directory String containing location of WinBUGS installation;
#' default "C:/Program Files (x86)/WinBUGS14/" is for Windows
#' 64-bit machines.
#' @param WinBUGS.debug Logical flag as to whether to close WinBUGS after
#' running (default \code{FALSE} implies WinBUGS is closed).
#' @param WinBUGS.code You can supply your own code file, especially useful if
#' you want to use informative priors for external information. The default
#' value is \code{NULL}, and the code file must be in the folder as specified by
#' \code{WinBUGS.code.location}. The effect of a \code{NULL} value that the
#' function \code{write.bugs.model} is called, and the WinBUGS code file is
#' generated automatically.
#' @param WinBUGS.code.location If \code{WinBUGS.code} is not \code{NULL} then
#' this folder is examined for the code file specified.  The default value is
#' \code{NULL}, whereby the code looks for the file in the current directory.
#' @param no.starting.value A list of strings denoting the objects we should
#' not initialise, usually because we set them or calculate them in a bespoke
#' WinBUGS code file. These objects will be set to \code{NA} by this function;
#' to use, create objects such as
#'                           \code{list("beta[2,2]","ax[1]")}
#' for example.
#' @param estimate.p Logical flag specifying whether or not to retain samples
#' for the posterior probabilities of presence for each cell,
#' default \code{FALSE}; be aware using \code{TRUE} can result in slow
#' interaction between R and WinBUGS.
#' @param estimate.u Logical flag specifying whether or not to retain samples
#' for the spatial random effects of presence for each cell,
#' default \code{FALSE}; be aware using \code{TRUE} can result in slow
#' interaction between R and WinBUGS.
#' @param silent Logical flag denoting whether the function runs silently or
#' not. Default is \code{TRUE}.
#'
#' The final inputs are optional and will be ignored otherwise. Specifically,
#' you have the option of setting up more than one \code{car.normal} in case you
#' have disconnected subsets of your neighbourhood structure. The vectors
#' \code{num}, \code{adj} and \code{u} are still used, and are then composed of
#' \code{NClique} parts; the start/end vectors
#' below are used to draw out the relevant data for each clique. Note that if
#' using separate cliques, \code{num} will have length n as before, but the
#' order of corresponding cells will be different (not just 1 to n); instead,
#' \code{num} specifies the numbers of neighbours in each cell in clique 1 (in
#' the same order as in the data set), and then the numbers of neighbours of
#' cells in clique 2, and so on, until the final clique. Vectors \code{adj}
#' and \code{u} are structured similarly. The start and end points of each
#' clique are specified using the following inputs:
#' @param u.clique.start Vector of length \code{NClique}; where each clique
#' starts in \code{u} (and \code{num}) [Note \code{NClique} is determined in
#' the code as the length of this input vector]
#' @param u.clique.end Vector of length \code{NClique}; where each clique ends
#' in \code{u} (and \code{num}).
#' @param adj.clique.start Where each clique starts in \code{adj}.
#' @param adj.clique.end Where each clique ends in \code{adj}.
#' @param clique Vector to say which clique each observation is in, length n.
#' @param clique.i Vector to say which element each observation is of its
#' clique, length n.
#'
#' If none of these six are supplied, they will be set in the code to enable a
#' "standard" (single clique) analysis.
#'
#' @return A list object as returned by the \code{bugs()} function in
#' \pkg{R2WinBUGS}; see the help file for that function for further details.
#'
#' The list is augmented by the response \code{y} and the climate variables
#' \code{x.clim}, primarily to aid the plotting functions, and by
#' the vector which.beta, of length p, which contains 1's and
#' 2's to signify which set of beta parameters to use; see
#' the description of beta below.
#'
#' The variables in WinBUGS are (i indexes the climate variable):
#' \describe{
#' \item{\code{beta[which.beta[i],i,1/2]}}{\code{which.beta} (see above)
#'                    specifies which of two sets of \code{beta} parameters to use.
#'                    There are two \code{betas} in order to facilitate constraining
#'                    either the \code{beta[i,1]} (up-slope) or \code{beta[i,2]} (down-slope)
#'                    is to be constrained to be close to its opposite. If
#'                    \code{beta[i,2]} is constrained, we'll use the first set (i.e.
#'                    \code{which.beta[i]==1}) as this models \code{beta[i,1]} directly and
#'                    the multiplicative difference between that and \code{beta[i,2]};
#'                    if \code{beta[i,1]} is constrained (\code{which.beta[i]==2}) then we
#'                    model \code{beta[i,2]} directly.}
#' \item{\code{ax}}{The x-coordinate of the apex for each climate variable.}
#' \item{\code{beta0}}{Scalar suggesting where the top-slicing should be applied.}
#' \item{\code{az}}{Scalar coordinate on the logit scale of the response axis apex.}
#' \item{\code{gamma}}{(Upper triangle) matrix of pairwise interaction parameters
#'                    between the climate variables; these are constrained
#'                    relative to the \code{betas} in order not to break the geometric
#'                    formula for a cone.}
#' }
#' @export
fit.bugs.env <- function(data,y,x.clim,x.nonclim=NULL,x.factor=NULL,
    not.spatial=FALSE,car.sigma=0.1,num,adj,u,
    prior.ax,prior.beta,prior.beta0.difference,constrain.beta,initial.pars.input,
    informative.priors=list(beta=FALSE, beta0=FALSE, ax=FALSE),
    burnin=5000,post.burnin=1000,chains=2,thin=1,
    working.directory=NULL,silent=TRUE,
    bugs.directory="C:/Program Files (x86)/WinBUGS14/",
    WinBUGS.debug=FALSE,WinBUGS.code=NULL,WinBUGS.code.location=NULL,
    no.starting.value=NULL,estimate.p=FALSE,estimate.u=FALSE,
    u.clique.start,u.clique.end,adj.clique.start,adj.clique.end,clique,clique.i){

    # Move to working directory or temporary directory when NULL
    inTempDir <- FALSE
    if(is.null(working.directory)) {
      # create a temp folder for WINBUGS in the current working directory
      working.directory <- getwd()
      working.directory <- paste(working.directory,"/temp",format(Sys.time(), "%Y%b%d%H%M"),sep="")
      if( !file.exists(working.directory) ) dir.create(working.directory)
      
      savedWD <- getwd()
      setwd(working.directory)
      on.exit(setwd(savedWD), add = TRUE)
      inTempDir <- TRUE
    }
    #y.name <- y
    y <- data[,y]
    n.x.clim <- length(x.clim)
    x.clim <- as.matrix(data[,x.clim])
    beta.top <- 2*n.x.clim
    # Now standardise the climate variables by mapping onto [0,1]
    x.clim.std <- apply(x.clim,2,function(x){
      x.min <- min(x)
      x.max <- max(x)
      return((x-x.min)/(x.max-x.min))
    })
    x.nonclim.names <- x.nonclim
    x.factor.names <- x.factor
    if(!is.null(x.nonclim)){
      if(!is.null(x.factor)){
        x.nonclim <- as.matrix(data[,x.nonclim])
        x.factor <- as.matrix(data[,x.factor])
        
      }else{
        x.nonclim <- as.matrix(data[,x.nonclim])
      }
      n.x.nonclim <- ncol(x.nonclim)
    }else{
      if(!is.null(x.factor)){
        x.factor <- as.matrix(data[,x.factor])
      }
    }
    if(!is.null(x.factor)){
      n.x.factor <- ncol(x.factor)
    }
    # Either use user specified WINBUGS code, or automatically generate it.
    # if WINBUGS code is provided manually, then use it; if not, automatically generate it.
    if(!is.null(WinBUGS.code)){
      WinBUGS.model <- WinBUGS.code
      file.copy(paste(WinBUGS.code.location,WinBUGS.code,sep="/"),
                working.directory)
    }else{
      write.bugs.model(n.x.clim,n.x.nonclim,n.x.factor,not.spatial,
                       working.directory=working.directory)
      WinBUGS.model <- "WINBUGS_code.txt"
    }
    if(!not.spatial){
      # Set the clique parameters if not supplied, to in effect run with only one
      # clique
      if(missing(u.clique.start)){
        u.clique.start <- c(1,NA)
      }
      if(missing(u.clique.end)){
        u.clique.end <- c(length(y),NA)
      }
      if(missing(adj.clique.start)){
        adj.clique.start <- c(1,NA)
      }
      if(missing(adj.clique.end)){
        adj.clique.end <- c(length(adj),NA)
      }
      NClique <- length(u.clique.start)-1
      if(missing(clique)){
        clique <- rep(1,length(y))
      }
      if(missing(clique.i)){
        clique.i <- rep(1,length(y))
      }
      if (missing(u) & NClique==1){
        u <- as.numeric(rep(NA,length(y)))
      }
      clique.length <- u.clique.end[-(NClique+1)]-u.clique.start[-(NClique+1)]+1
      adj.length <- adj.clique.end[-(NClique+1)]-adj.clique.start[-(NClique+1)]+1
    }
    if(is.null(informative.priors$beta)){
      informative.priors$beta <- FALSE
    }
    if(is.null(informative.priors$beta0)){
      informative.priors$beta0 <- FALSE
    }
    if(is.null(informative.priors$ax)){
      informative.priors$ax <- FALSE
    }
    if(!not.spatial){
      # Can't give a clique with one observation to car.normal, so ignore these
      NonSingletonClique <- sum(clique.length>1.5)
      nonsingleton.clique.list <- c(which(clique.length>1.5),NA)
    }
    if(missing(initial.pars.input)){
        if(missing(constrain.beta)){
            initial.object <- generate.initial.values(y=y,x.clim=x.clim.std,
                x.nonclim=x.nonclim,x.factor=x.factor)
        }else{
            initial.object <- generate.initial.values(y=y,x.clim=x.clim.std,
                x.nonclim=x.nonclim,x.factor=x.factor,constrain.beta=constrain.beta)
        }
        initial.pars <- initial.object$initial.pars
        constrain.beta <- initial.object$constrain.beta
    }else{ # If initial values supplied, use those
        initial.pars <- initial.pars.input
    }
    if(missing(constrain.beta)){
        constrain.beta <- matrix(rep(FALSE,2*n.x.clim),ncol=2)
    }
    if(missing(prior.ax)){
        if(informative.priors$ax){
            ax.mu <- initial.pars[(beta.top+1):(beta.top+n.x.clim)]
            ax.var <- rep(0.001,n.x.clim)
        }else{
            ax.mu <- rep(0.5,n.x.clim)
            ax.var <- rep(1,n.x.clim)
        }
    }else{
        if(is.null(prior.ax$mean)){
            if(informative.priors$ax){
                ax.mu <- initial.pars[(beta.top+1):(beta.top+n.x.clim)]
            }else{
                ax.mu <- rep(0.5,n.x.clim)
            }
        }else{
            ax.mu <- prior.ax$mean
        }
        if(is.null(prior.ax$var)){
            if(informative.priors$ax){
                ax.var <- rep(0.001,n.x.clim)
            }else{
                ax.var <- rep(1,n.x.clim)
            }
        }else{
            ax.var <- prior.ax$var
        }
        if(any(is.na(ax.mu))){
            ax.mu[is.na(ax.mu)] <- 0.5
        }
        if(any(is.na(ax.var))){
            ax.var[is.na(ax.var)] <- 1
        }
    }
    alpha1 <- matrix(rep(1,2*n.x.clim),ncol=2)
    alpha2 <- matrix(rep(1,2*n.x.clim),ncol=2)
    if(missing(prior.beta)){
        if(informative.priors$beta){
            beta.mu <- matrix(exp(initial.pars[1:beta.top]),byrow=TRUE,ncol=2)
            beta.var <- matrix(rep(10,2*n.x.clim),ncol=2)
        }else{
            beta.mu <- matrix(rep(100,2*n.x.clim),ncol=2)
            beta.var <- matrix(rep(10000,2*n.x.clim),ncol=2)
        }
    }else{
        beta.mu <- matrix(exp(prior.beta$mean),ncol=2,byrow=TRUE)
        beta.var <- matrix(prior.beta$var,ncol=2,byrow=TRUE)
        beta.mu[is.na(beta.mu)] <- 100
        beta.var[is.na(beta.var)] <- 10000
    }
    if(missing(prior.beta0.difference)){
        if(informative.priors$beta0){
            beta0.mu <- initial.pars[beta.top+n.x.clim+1]
            beta0.var <- 0.1
        }else{
            beta0.mu <- 1
            beta0.var <- 10
        }
    }else{
        beta0.mu <- prior.beta0.difference$mean
        beta0.var <- prior.beta0.difference$var
    }
    if(any(ax.mu< -0.9)){
        for(i in n.x.clim){
            if(ax.mu[i]< -0.9){
                ax.mu[i] <- -0.9
            }
        }
    }
    if(any(ax.mu> 1.9)){
        for(i in n.x.clim){
            if(ax.mu[i]> 1.9){
                ax.mu[i] <- 1.9
            }
        }
    }
    if(!not.spatial){
      car.tau <- 1/car.sigma^2
    }
    WinBUGS.data <- list(y=y,
                     x.clim=x.clim.std,
                     N=length(y),
                     NEnv=ncol(x.clim.std),
                     ax.mu=ax.mu,
                     ax.var=ax.var,
                     beta.mu=beta.mu,
                     beta.var=beta.var,
                     beta0.mu=beta0.mu,
                     beta0.var=beta0.var)
    if(!not.spatial){
      WinBUGS.data$NClique=NClique
      WinBUGS.data$NonSingletonClique=NonSingletonClique
      WinBUGS.data$nonsingleton.clique.list=nonsingleton.clique.list
      WinBUGS.data$u.clique.start=u.clique.start
      WinBUGS.data$u.clique.end=u.clique.end
      WinBUGS.data$adj.clique.start=adj.clique.start
      WinBUGS.data$adj.clique.end=adj.clique.end
      WinBUGS.data$clique=clique
      WinBUGS.data$clique.i=clique.i
      WinBUGS.data$car.tau=car.tau
      WinBUGS.data$adj=adj
      WinBUGS.data$num=num
      WinBUGS.data$u=u
    }
    if(n.x.clim==1){
        WinBUGS.monitor <- c("beta","beta0","ax","az","deviance")
    }else{
        WinBUGS.monitor <- c("beta","beta0","gamma","ax","az","deviance")
    }
    if(estimate.p){
        WinBUGS.monitor <- c(WinBUGS.monitor,"p")
    }
    if(!not.spatial){
      if(estimate.u){
        WinBUGS.monitor <- c(WinBUGS.monitor,"u")
      }
    }
    if(!is.null(x.nonclim)){
        WinBUGS.data$x.nonclim <- as.matrix(x.nonclim)
        WinBUGS.data$NnonEnv <- ncol(x.nonclim)
        WinBUGS.monitor <- c(WinBUGS.monitor,"nonbeta")
    }
    if(!is.null(x.factor)){
        WinBUGS.data$x.factor <- as.matrix(x.factor)
        WinBUGS.data$x.factor.lengths <- numeric(ncol(WinBUGS.data$x.factor))
        for(i in 1:ncol(WinBUGS.data$x.factor)){
            WinBUGS.data$x.factor.lengths[i] <- nlevels(x.factor[,i])
            WinBUGS.data$x.factor[,i] <- as.numeric(x.factor[,i])
        }
        WinBUGS.data$Nfac <- ncol(WinBUGS.data$x.factor)
        WinBUGS.monitor <- c(WinBUGS.monitor,"nonfac")
    }
    WinBUGS.inits <- list()
    n.data <- length(y)
    beta.mat <- matrix(exp(initial.pars[1:beta.top]),byrow=TRUE,ncol=2)
    beta.mat[beta.mat>1000] <- 1000
    ax.vec <- initial.pars[(beta.top+1):(beta.top+n.x.clim)]
    if(n.x.clim==1){
        if(ax.vec[1]>1.1){
            ax.vec[1] <- 1.1
        }
        if(ax.vec[1]< -0.1){
            ax.vec[1] <- -0.1
        }
    }else{
        if(any(ax.vec< -0.8)){
            for(i in n.x.clim){
                if(ax.vec[i]< -0.8){
                    ax.vec[i] <- -0.8
                }
            }
        }
        if(any(ax.vec> 1.8)){
            for(i in n.x.clim){
                if(ax.vec[i]> 1.8){
                    ax.vec[i] <- 1.8
                }
            }
        }
    }
    az <- initial.pars[beta.top+n.x.clim+2]
    beta0 <- az-exp(initial.pars[beta.top+n.x.clim+1])
    ax.mat <- matrix(rep(ax.vec,each=n.data),ncol=n.x.clim)
    if(n.x.clim!=1){
        gamma.mat <- array(NA,dim=c(n.x.clim,n.x.clim))
        gamma.mat[t(upper.tri(gamma.mat))] <- initial.pars[(beta.top+n.x.clim+3):length(initial.pars)]
        gamma.mat <- t(gamma.mat)
        gamma.mat[gamma.mat>20] <- 20 # To stop Inf in the line below (gives same result - 1)
        gamma.mat <- exp(gamma.mat)/(1+exp(gamma.mat))
    }
    x.clim.a.centre <- x.clim.std-ax.mat
    x.mat.ind <- (x.clim.a.centre > 0)+1
    beta.mat.prod <- array(dim=c(n.data,n.x.clim))
    for(i in 1:n.x.clim){
        beta.mat.prod[,i] <- beta.mat[i,x.mat.ind[,i]]
    }
    # the contributions from the main effects of climate
    x.sum.main <- rowSums(beta.mat.prod*x.clim.a.centre*x.clim.a.centre)
    if(n.x.clim!=1){ # no interactions if p==1
        # Combinations of row, column with values
        gamma.rc <- which(!is.na(gamma.mat),arr.ind=TRUE)
        gamma.min.1 <- beta.mat[gamma.rc[,"row"],1]*beta.mat[gamma.rc[,"col"],1]
        gamma.min.2 <- beta.mat[gamma.rc[,"row"],2]*beta.mat[gamma.rc[,"col"],2]
        gamma.min <- -(2/(n.x.clim-1))*sqrt(pmin(gamma.min.1,gamma.min.2))
        gamma.max.1 <- beta.mat[gamma.rc[,"row"],1]*beta.mat[gamma.rc[,"col"],2]
        gamma.max.2 <- beta.mat[gamma.rc[,"row"],2]*beta.mat[gamma.rc[,"col"],1]
        gamma.max <- (2/(n.x.clim-1))*sqrt(pmin(gamma.max.1,gamma.max.2))
        n.interactions <- nrow(gamma.rc)
        gamma.rc.row <- gamma.rc[,"row"]
        gamma.rc.col <- gamma.rc[,"col"]
        gamma.mat.1 <- gamma.mat
        for(i in 1:n.interactions){
            gamma.mat.1[gamma.rc.row[i],gamma.rc.col[i]] <- gamma.min[i]+gamma.mat.1[gamma.rc.row[i],gamma.rc.col[i]]*(gamma.max[i]-gamma.min[i])
        }
        int.mat <- x.clim.a.centre[,gamma.rc]
        int.mat.1 <- as.matrix(int.mat[,1:n.interactions])
        int.mat.2 <- as.matrix(int.mat[,(n.interactions+1):(2*n.interactions)])
        # the contributions from the interactions of climate variables
        x.sum.int <- rowSums(int.mat.1*int.mat.2*matrix(rep(gamma.mat.1[gamma.rc],n.data),ncol=n.interactions,byrow=TRUE))
    }else{
        x.sum.int <- 0
    }
    x.envelope <- pmin(az-(sqrt(x.sum.main+x.sum.int)),beta0)
    data.temp <- data.frame(y=y,x.envelope=x.envelope)
    glmfit <- glm(y~offset(x.envelope)-1,family=binomial,data=data.temp)
    pred.u <- predict(glmfit,type="link")
    mean.presence <- mean(data.temp$y)
    logit.mean.presence <- log(mean.presence/(1-mean.presence))
    if(!not.spatial){
      # Generate initial values for the spatial random effect u's
      u.inits <- rep(NA,length(u))
      # Dividing by 20 produces sensible values here
      u.inits[data.temp$y==1] <- (max(pred.u)-pred.u[data.temp$y==1]-logit.mean.presence)/20
      u.inits[data.temp$y==0] <- (pred.u[data.temp$y==0]-logit.mean.presence-max(pred.u))/20
      u.inits[!is.na(u)] <- NA
    }
    for(i in 1:chains){
        if(i==1){
            if(n.x.clim==1){
                WinBUGS.inits[[i]] <- list(beta0.diff=initial.pars[beta.top+n.x.clim+1],
                    beta=matrix(beta.mat,ncol=2),
                    ax=ax.vec,az=az)
            }else{
                gamma.mat[gamma.mat < 1e-20] <- 1e-20
                WinBUGS.inits[[i]] <- list(beta0.diff=initial.pars[beta.top+n.x.clim+1],
                    beta=matrix(beta.mat,ncol=2),
                    gamma.temp=gamma.mat,
                    ax=ax.vec,az=az)
            }
        }else{
            if(n.x.clim==1){
                WinBUGS.inits[[i]] <- list(beta0.diff=initial.pars[beta.top+n.x.clim+1]*runif(1,0.9,1.1),
                    beta=matrix(beta.mat*matrix(runif(prod(dim(beta.mat)),0.9,1),ncol=2),ncol=2),
                    ax=ax.vec*runif(length(ax.vec),0.9,1.1),az=az*runif(1,0.9,1.1))
            }else{
                WinBUGS.inits[[i]] <- list(beta0.diff=initial.pars[beta.top+n.x.clim+1]*runif(1,0.9,1.1),
                    beta=matrix(beta.mat*matrix(runif(prod(dim(beta.mat)),0.9,1),ncol=2),ncol=2),
                    gamma.temp=gamma.mat*matrix(runif(prod(dim(gamma.mat)),0.9,1),ncol=n.x.clim),
                    ax=ax.vec*runif(length(ax.vec),0.9,1.1),az=az*runif(1,0.9,1.1))
            }
        }
        if(!silent){
            print(paste("Initial values for chain",i,":"))
            print(WinBUGS.inits[[i]])
        }
        if(!not.spatial){
          WinBUGS.inits[[i]]$u <- u.inits
        }
        if(!is.null(no.starting.value)){
            for(j in 1:length(no.starting.value)){
                eval(parse(text=paste("WinBUGS.inits[[",i,"]]$",no.starting.value[[j]]," <- NA",sep="")))
            }
        }
    }
    # Which is the good one? Choose either if both, but ensure large variance above
    which.beta <- rep(2,n.x.clim)
    which.beta[rowSums(constrain.beta)!=0] <- 1
    if(!silent){
        print(which.beta)
    }
    WinBUGS.data$which.beta <- which.beta

    #require(R2WinBUGS)
    if(!silent){
        cat("Starting WinBUGS run - opening WinBUGS now...\n")
    }
    results <- bugs(WinBUGS.data,WinBUGS.inits,WinBUGS.monitor,WinBUGS.model,
        chains,post.burnin+burnin,burnin,thin,
        bugs.directory=bugs.directory,
        working.directory=working.directory,codaPkg=FALSE,debug=WinBUGS.debug,DIC=FALSE)
    if(!silent){
        cat("WinBUGS run completed.\n")
    }  
    parvec <- c(log(t(results$mean$beta)),results$mean$ax,log(results$mean$az-results$mean$beta0),results$mean$az)
    beta.mat <- results$mean$beta
    if(n.x.clim!=1){
        gamma.mat <- array(NA,dim=c(n.x.clim,n.x.clim))
        gamma.mat[t(upper.tri(gamma.mat))] <- c(t(results$mean$gamma))[1:(n.x.clim*(n.x.clim-1)/2)]
        gamma.mat <- t(gamma.mat)
        gamma.rc <- which(!is.na(gamma.mat),arr.ind=TRUE)
        gamma.min.1 <- beta.mat[gamma.rc[,"row"],1]*beta.mat[gamma.rc[,"col"],1]
        if(n.x.clim > 2.5){
            gamma.min.2 <- diag(beta.mat[gamma.rc[,"row"],which.beta[gamma.rc[,"row"]]])*diag(beta.mat[gamma.rc[,"col"],which.beta[gamma.rc[,"col"]]])
        }else{
            gamma.min.2 <- beta.mat[gamma.rc[,"row"],which.beta[gamma.rc[,"row"]]]*beta.mat[gamma.rc[,"col"],which.beta[gamma.rc[,"col"]]]
        }
        gamma.min <- -(2/(n.x.clim-1))*sqrt(pmin(gamma.min.1,gamma.min.2))
        if(n.x.clim > 2.5){
            gamma.max.1 <- beta.mat[gamma.rc[,"row"],1]*diag(beta.mat[gamma.rc[,"col"],which.beta[gamma.rc[,"col"]]])
            gamma.max.2 <- diag(beta.mat[gamma.rc[,"row"],which.beta[gamma.rc[,"row"]]])*beta.mat[gamma.rc[,"col"],1]
        }else{
            gamma.max.1 <- beta.mat[gamma.rc[,"row"],1]*beta.mat[gamma.rc[,"col"],which.beta[gamma.rc[,"col"]]]
            gamma.max.2 <- beta.mat[gamma.rc[,"row"],which.beta[gamma.rc[,"row"]]]*beta.mat[gamma.rc[,"col"],1]
        }
        gamma.max <- (2/(n.x.clim-1))*sqrt(pmin(gamma.max.1,gamma.max.2))
        for(i in 1:(n.x.clim*(n.x.clim-1)/2)){
            gamma.mat[gamma.rc[i,"row"],gamma.rc[i,"col"]] <- (gamma.mat[gamma.rc[i,"row"],gamma.rc[i,"col"]]-gamma.min[i])/(gamma.max[i]-gamma.min[i])
        }
        gamma.mat <- log(gamma.mat/(1-gamma.mat))
        gamma.vec <- t(gamma.mat)[lower.tri(t(gamma.mat))]
        parvec <- c(parvec,gamma.vec)
    }
    results$par <- parvec
    results$y <- y
    results$x.clim <- x.clim
    results$x.nonclim <- x.nonclim
    results$x.factor <- x.factor
    results$which.beta <- which.beta

    return(results)
}
