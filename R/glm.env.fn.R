#' Calculate deviance of plateau envelope via GLM
#'
#' \code{glm.env.fn} fits the plateau envelope as part of a GLM. It is called
#' by the function \code{fit.glm.env} and is not intended for use interactively.
#' At present, only binary (logistic) response models for presence/absence
#' data are implemented.
#'
#' @param pars The vector of envelope parameters, length 2p+p+2+p(p-1)/2.
#' @param y The binary response variable (taking values 0 or 1 for absence and
#' presence respectively).
#' @param x.clim The n by p matrix of climate covariates.
#' @param x.nonclim The n by p2 matrix of climate covariates.
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
glm.env.fn <- function(pars,y,x.clim,x.nonclim=NULL,constrain.beta=FALSE,
    slope.limit=7){
    n.x.clim <- ncol(x.clim)
    env.fn.object <- env.fn(pars=pars,x.clim=x.clim,slope.limit=slope.limit)
    x.envelope <- env.fn.object$x.envelope
    if(sum(env.fn.object$x.envelope)==0){
        return(9e100)
    }else{
        beta.mat <- env.fn.object$beta.mat
        ax.vec <- env.fn.object$ax.vec
        # exp.offset <- exp(offset(x.envelope))
        if(is.null(x.nonclim)){
            glm.temp <- glm(y~offset(x.envelope)-1,family=binomial)
        }else{
            glm.temp <- glm(y~offset(x.envelope)-1+x.nonclim,family=binomial)
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
