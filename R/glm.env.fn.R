#' Calculate deviance of plateau envelope via GLM
#'
#' \code{glm.env.fn} fits the plateau envelope as part of a GLM. It is called
#' by the function \code{fit.glm.env} and is not intended for use interactively.
#' At present, only binary (logistic) response models for presence/absence
#' data are implemented.
#'
#' @param pars The vector of envelope parameters, length 2p+p+2+p(p-1)/2.
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
#' @param x.nonclim.formula The pre-specified formula if non-climate variables
#' OR factor variables are included.
#' @param constrain.beta Should a ridge penalty be imposed on the betas
#' (slopes)? Takes the form either of: a matrix of dimension p by 2
#' (specifying \code{TRUE} for a ridge penalty for a given beta, or
#' \code{FALSE} otherwise); or of a single value \code{FALSE} to imply no
#' ridge penalties at all.
#' @param slope.limit Scalar putting an upper bound on the envelope slopes;
#' limit is approximately \code{exp(slope.limit)}.
#' @return The (scalar) deviance for fitted GLM using specified envelope
#' (possibly including a ridge penalty if set).
#' @export
glm.env.fn <- function(pars,data,y,x.clim,x.nonclim=NULL,x.nonclim.formula=NULL,
    x.factor=NULL,constrain.beta=FALSE,slope.limit=7){
    n.x.clim <- ncol(x.clim)
    env.fn.object <- env.fn(pars=pars,x.clim=x.clim,slope.limit=slope.limit)
    x.envelope <- env.fn.object$x.envelope
    y.name <- y
    y <- data[,y]
    if(mean(env.fn.object$x.envelope) < -2000){
        return(9e100)
    }else{
        beta.mat <- env.fn.object$beta.mat
        ax.vec <- env.fn.object$ax.vec
        if(is.null(x.nonclim) && is.null(x.factor)){
            glm.temp <- glm(y~offset(x.envelope)-1,family=binomial)
        }else{
            if(!is.null(x.nonclim)){
              if(!is.null(x.factor)){
                tempdata <- data[,c(y.name,x.nonclim,x.factor)]
                tempdata$x.envelope <- x.envelope
              }else{
                tempdata <- data[,c(y.name,x.nonclim)]
                tempdata$x.envelope <- x.envelope
              }
            }else{
              if(!is.null(x.factor)){
                tempdata <- data[,c(y.name,x.factor)]
                tempdata$x.envelope <- x.envelope              }
            }
            if(!is.null(x.factor)){
              if(length(x.factor)==1){
                temp.nlevels <- nlevels(tempdata[,x.factor])
                tempdata[,x.factor] <- (contr.sum(temp.nlevels)+contr.helmert(temp.nlevels))[tempdata[,x.factor],]
              }else{
                for(i in 1:length(x.factor)){
                  temp.nlevels <- nlevels(tempdata[,x.factor[i]])
                  tempdata[,x.factor[i]] <- (contr.sum(temp.nlevels)+contr.helmert(temp.nlevels))[tempdata[,x.factor[i]],]
                }
              }
            }
          #print(names(tempdata))
          #print(str(tempdata))
          #print("glm1")
          #print(summary(glm(x.nonclim.formula,data=tempdata,family=binomial)))
          #print((glm(x.nonclim.formula,data=tempdata,family=binomial)$contrasts))
          #print("glm2")
          #print(glm(x.nonclim.formula,data=tempdata,family=binomial,contrasts=list(FAC1="contr.poly")))
          #print((glm(x.nonclim.formula,data=tempdata,family=binomial)$contrasts))
          glm.temp <- try(glm(x.nonclim.formula,data=tempdata,family=binomial))
        }
        glm.deviance <- glm.temp$deviance
        # Impose penalty on the difference between the two betas for the same
        # if there is a TRUE in the corresponding row of constrain.beta
        if(any(constrain.beta)){
            for(i in 1:n.x.clim){
                if(any(constrain.beta[i,])){
                    #glm.deviance <- glm.deviance+abs(diff(beta.mat[i,]))#^2
                    glm.deviance <- glm.deviance+(diff(beta.mat[i,]))^2
                }
            }
        }
        # Also impose penalty to stop the ax values from straying beyond [-1,2]
        if(any(ax.vec< -0.9)){
            for(i in 1:n.x.clim){
                if(ax.vec[i]< -0.9){
                    glm.deviance <- glm.deviance+abs(ax.vec[i]+0.9)#^2
                }
            }
        }
        if(any(ax.vec>1.9)){
            for(i in 1:n.x.clim){
                if(ax.vec[i]>1.9){
                    glm.deviance <- glm.deviance+abs(ax.vec[i]-1.9)#^2
                }
            }
        }
        # Stop beta becoming too large (causes overflow errors)
        if(any(pars[1:(2*n.x.clim)]>slope.limit)){
            for(i in 1:(2*n.x.clim)){
                if(pars[i]>slope.limit){
                    glm.deviance <- glm.deviance+0.2*(pars[i]-slope.limit)^2 # or abs
                }
            }
        }
        return(glm.deviance)
    }
}
