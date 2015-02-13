#' Generate initial values for plateau envelopes
#'
#' \code{generate.initial.values} generates initial values for the plateau
#' envelopes. It is called by \code{fit.glm.env} and \code{fit.bugs.env}, and is
#' not intended for use interactively.
#'
#' @param y The binary response variable (taking values 0 or 1 for absence and
#' presence respectively).
#' @param x.clim The n by p matrix of climate covariates.
#' @param x.nonclim The n by p2 matrix of non-climate covariates.
#' @param constrain.beta Should ridge penalty be imposed on the betas (slopes)?
#' @param random Logical flag on whether or not to generate random starting values
#' @param pars The vector of envelope parameters, length 2p+p+2+p(p-1)/2
#' @return If \code{random==TRUE}, the output will be a single vector of
#' parameter values in the order required for \code{glm.env.fn}; if
#' \code{random==FALSE}, the output will be a list with two elements:
#' \code{initial.pars}, a single vector of parameter values in the order
#' required for \code{glm.env.fn}, and \code{constrain.beta}, also as for
#' \code{glm.env.fn}.
#' @export
generate.initial.values <- function(y,x.clim,x.nonclim=NULL,
    constrain.beta=FALSE,random=FALSE,pars){

    x.clim <- as.matrix(x.clim)
    n.x.clim <- ncol(x.clim)
    beta.init <- numeric(2*n.x.clim)
    ax.init <- numeric(n.x.clim)
    beta0.init <- numeric(n.x.clim)
    az.init <- numeric(n.x.clim)
    len <- 101
    gamdata <- data.frame(y=y)
    for(i in 1:n.x.clim){
        gamdata[paste("x",i,sep="")] <- x.clim[,i]
    }
    gamdata$x.nonclim <- x.nonclim
    x <- seq(0,1,length=len)
    x.plot <- array(NA,dim=c(len^n.x.clim,n.x.clim))
    for(i in 1:n.x.clim){
        x.plot[,i] <- rep(x,times=len^(i-1),each=len^(n.x.clim-i))
    }
    newgamdat <- list()
    for(i in 1:n.x.clim){
        newgamdat[[paste("x",i,sep="")]] <- x.plot[,i]
    }
    if(!is.null(x.nonclim)){
        newgamdat$x.nonclim <- matrix(rep(colMeans(x.nonclim),len^n.x.clim),
            ncol=ncol(x.nonclim),byrow=TRUE)
    }
    newgamdatmat <- array(NA,dim=c(len^n.x.clim,n.x.clim))
    for(i in 1:n.x.clim){
        newgamdatmat[,i] <- newgamdat[[i]]
    }
    if(random==TRUE){
        if(missing(pars)){
            pars.object <- generate.initial.values(y=y,x.clim=x.clim,
                constrain.beta=constrain.beta)
            pars <- pars.object$initial.pars
        }
        beta.init <- pars[1:(2*n.x.clim)]*c(abs(rnorm(2*n.x.clim,1,0.25)))
        #beta.init[beta.init>10] <- 10
        ax.init <- pars[(2*n.x.clim)+(1:n.x.clim)]*c(abs(rnorm(n.x.clim,1,0.1)))
        beta0.init <- pars[(3*n.x.clim)+1]*c(abs(rnorm(1,1,0.1)))
        az.init <- pars[(3*n.x.clim)+2]*c(abs(rnorm(1,1,0.1)))
        if(n.x.clim!=1){
            gamma.init <- pars[(3*n.x.clim)+2+1:(n.x.clim*(n.x.clim-1)/2)]*c(rnorm(n.x.clim*(n.x.clim-1)/2,1,5))
            pars <- c(beta.init,ax.init,beta0.init,az.init,gamma.init)
        }else{
            pars <- c(beta.init,ax.init,beta0.init,az.init)
        }
        return(pars)
    }
    if(n.x.clim!=1){
        # Try multivariate GAM fit to set starting values for x.clim>=2
        if(is.null(x.nonclim)){
            gam.formula <- paste("y~te(",paste("x",1:(n.x.clim-1),",",sep="",collapse=""),"x",
                n.x.clim,",sp=rep(0.01,",n.x.clim,"))",sep="",collapse="")
        }else{
            gam.formula <- paste("y~te(",paste("x",1:(n.x.clim-1),",",sep="",collapse=""),"x",
                n.x.clim,",sp=rep(0.01,",n.x.clim,"))+x.nonclim",sep="",collapse="")
        }
        gam.formula <- as.formula(gam.formula)
        gamfit3 <- gam(gam.formula,family=binomial,data=gamdata)
        predgam3a <- predict(gamfit3,newdata=newgamdat,type="link")
        maxgam <- array(NA,dim=c(n.x.clim,len^(n.x.clim-1)))
        fitintercept <- array(NA,dim=c(n.x.clim,len^(n.x.clim-1)))
        fitslope <- array(NA,dim=c(n.x.clim,len^(n.x.clim-1)))
        multiplier.mat <- matrix(rep((n.x.clim-1):1,each=len^n.x.clim),ncol=n.x.clim-1)
        multiplier.mat <- (len+1)^(multiplier.mat-1)
        for(j in 1:n.x.clim){
            newgamdatmatnoj <- newgamdatmat[,-j]
            newgamdatmatnoj.rowscores <- rowSums(multiplier.mat*newgamdatmatnoj)
            uniquenoj <- unique(newgamdatmatnoj.rowscores)
            for(i in 1:(len^(n.x.clim-1))){
                the.index <- which(newgamdatmatnoj.rowscores==uniquenoj[i])
                the.y <- predgam3a[the.index]
                the.x <- newgamdatmat[the.index,j]
                maxgam[j,i] <- which.max(the.y)
                the.lm <- lm(the.y~the.x)
                fitintercept[j,i] <- coef(the.lm)[1]
                fitslope[j,i] <- coef(the.lm)[2]
            }
            # Slopes often too shallow, so arbitrarily increase
            beta.init[1+(j-1)*2] <- log(2*abs(max(fitslope[j,])))
            beta.init[2+(j-1)*2] <- log(2*abs(min(fitslope[j,])))
            ax.init[j] <- mean(x[maxgam[j,]])
        }
        # Try univariate fits to set starting values az and beta0
        for(i in 1:n.x.clim){
            # Calling fit.glm.env recursively
        print("Before glm call in generate.initial.values")
            fitted.envelope.glm <- fit.glm.env(y=y,x.clim=x.clim[,i])
        print("Before glm call in generate.initial.values")
            if(fitted.envelope.glm$par[1]>beta.init[1+(i-1)*2]){
                beta.init[1+(i-1)*2] <- fitted.envelope.glm$par[1]
            }
            if(fitted.envelope.glm$par[2]>beta.init[2+(i-1)*2]){
                beta.init[2+(i-1)*2] <- fitted.envelope.glm$par[2]
            }
            beta0.init[i] <- fitted.envelope.glm$par[4]
            az.init[i] <- fitted.envelope.glm$par[5]
        }
        beta0.init <- log(mean(exp(beta0.init)))
        az.init <- mean(az.init)
        beta.mat <- matrix(beta.init,ncol=2,byrow=TRUE)
        # Try setting initial values for the interaction parameters by
        # looking at the correlations between the climate variables when y
        # is 1 (i.e. where there are presences observed)
        cor.mat <- cor(x.clim[y==1,])
        # Map these correlations on [-1,1] to [0,1]
        cor.mat <- (1+cor.mat)/2
        cor.vec <- -cor.mat[t(upper.tri(cor.mat))]
        cor.mat[!upper.tri(cor.mat)] <- NA
        cor.rc <- which(!is.na(cor.mat),arr.ind=TRUE)
        cor.mat <- log(cor.mat/(1-cor.mat))
        initial.pars <- c(beta.init,ax.init,beta0.init,az.init,cor.mat[t(lower.tri(cor.mat))])
        constrain.beta <- matrix(rep(FALSE,n.x.clim*2),ncol=2)
        for(i in 1:n.x.clim){
            if(ax.init[i]+exp(beta0.init)/beta.init[2+(i-1)*2]>=0.99){
                constrain.beta[i,2] <- TRUE
            }else{
                if(ax.init[i]-exp(beta0.init)/beta.init[1+(i-1)*2]<=0.01){
                    constrain.beta[i,1] <- TRUE
                }
            }
        }
    }else{ # when n.x.clim==1, use GAM
        #require(mgcv)
        if(is.null(x.nonclim)){
            gam.formula <- as.formula("y~s(x1)")
        }else{
            gam.formula <- as.formula("y~s(x1)+x.nonclim")
        }
        unigamfit <- gam(gam.formula,family=binomial,data=gamdata)
        plot(unigamfit)
        predunigamfit <- predict(unigamfit,newdata=newgamdat,type="link")
        maxfun <- which.max(predunigamfit)
        ax.init <- x[maxfun]
        if(maxfun==1){
            beta.init[2] <- log((predunigamfit[maxfun]-predunigamfit[length(x)])/(x[length(x)]-x[maxfun]))
            beta.init[1] <- beta.init[2]
        }else{
            if(maxfun==len){
                beta.init[1] <- log((predunigamfit[maxfun]-predunigamfit[1])/(x[maxfun]-x[1]))
                beta.init[2] <- beta.init[1]
            }else{
                beta.init[1] <- log((predunigamfit[maxfun]-predunigamfit[1])/(x[maxfun]-x[1]))
                beta.init[2] <- log((predunigamfit[maxfun]-predunigamfit[length(x)])/(x[length(x)]-x[maxfun]))
            }
        }
        beta.init[beta.init>10] <- 10
        initial.pars <- c(beta.init,ax.init,0,2)
        constrain.beta <- matrix(rep(FALSE,2),ncol=2)
        if(ax.init+exp(beta0.init)/beta.init[2]>=0.99){
            constrain.beta[i,2] <- TRUE
        }else{
            if(ax.init-exp(beta0.init)/beta.init[1]<=0.01){
                constrain.beta[i,1] <- TRUE
            }
        }
    }

    initial.pars.and.constraints <- list()
    initial.pars.and.constraints$initial.pars <- initial.pars
    initial.pars.and.constraints$constrain.beta <- constrain.beta
    print("Initial values selected :")
    print(initial.pars)
    print("Constrain beta? :")
    print(constrain.beta)
    return(initial.pars.and.constraints)
}
