#' Plot the plateau envelopes produced by \code{fit.glm.env} and
#' \code{fit.bugs.env}.
#'
#' \code{envelope.plot} plots the plateau envelopes produced by the functions
#' \code{fit.glm.env} and \code{fit.bugs.env}.
#'
#' @param envelope.fit The fitted object from either \code{fit.glm.env} or
#' \code{fit.bugs.env}.
#' @param type String taking values \code{"persp"} or \code{"contour"}
#' to specify which type of plot is required; ignored if fitted object has
#' only one climate variable (when a univariate plot with rug is plotted).
#' @param x.plot.lims Vector of two numerics containing the endpoints for the
#' plotting regions in the (scaled) climate variables.
#' @param len The number of points to use for plotting in each dimension.
#' @param plot.vars If p>2 (currentl implying p==3), need to know which
#' two variables to plot; this vector of two climate column names makes
#' the selection.
#' @param x.labels Labels to plot on the x-axes, notionally the variable
#' names; if omitted, the code will use "Variable 1" etc.
#' @param fix.values Vector of numeric length p with named elements
#' corresponding to the values (and names) of the non-plotted variables;
#' ignored if p<3.
#' @param close.points String indicating how to plot points in contour or
#' perspective plots: the value "shade" plots filled circles for points
#' which are "close" to the chosen values of the held out variables, and
#' empty circles otherwise; the value "plot" is as for shade but without
#' the empty circles; and for any other value of \code{close.points},
#' the full data set will be plotted as empty circles.
#' @param contour.levels If specified, the probability values to be plotted
#' in the contour plots; this is especially important when using
#' \code{plot.envelope} to make a "video".
#' @return For the case of a single climate variable, the plot produced is a
#' line graph, possibly with rugplot and points. For two climate covariates,
#' either a contour plot or a perspective plot is produced. For three climate
#' covariates, the same choice of contour or perspective plot is made for two
#' of the variables for a fixed value of the third covariate.
#' @seealso \code{link{map.plot}}, \code{link{fit.bugs.env}}
#' @export
envelope.plot <- function(envelope.fit,type="persp",x.plot.lims=c(0,1),len=100,
    plot.vars,x.labels,fix.values,close.points="",contour.levels,...){
    x.clim <- as.matrix(envelope.fit$x.clim)
    x.clim.names <- colnames(x.clim)
    n.x.clim <- ncol(x.clim)
    # Now standardise the climate variables by mapping onto [0,1]
    x.clim.std <- apply(x.clim,2,function(x){
        x.min <- min(x)
        x.max <- max(x)
        return((x-x.min)/(x.max-x.min))
    })
    beta.top <- 2*n.x.clim
    beta.mat <- matrix(exp(envelope.fit$par[1:beta.top]),byrow=TRUE,ncol=2)
    ax <- envelope.fit$par[(beta.top+1):(beta.top+n.x.clim)]
    az <- envelope.fit$par[beta.top+n.x.clim+2]
    beta0 <- az-exp(envelope.fit$par[beta.top+n.x.clim+1])
    if(!n.x.clim==1){
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
    if(missing(x.labels)){
        x.labels <- paste("Variable",1:n.x.clim)
    }
    x.plot.seq <- seq(x.plot.lims[1],x.plot.lims[2],length=len)
    x.plot <- array(NA,dim=c(len^n.x.clim,n.x.clim))
    for(i in 1:n.x.clim){
        x.plot[,i] <- rep(x.plot.seq,times=len^(i-1),each=len^(n.x.clim-i))
    }
    ax.mat <- matrix(rep(ax,each=nrow(x.plot)),ncol=n.x.clim)
    x.plot.a.centre <- x.plot-ax.mat
    x.mat.ind <- (x.plot.a.centre > 0)+1
    beta.mat.prod <- array(dim=c(nrow(x.plot),n.x.clim))
    for(i in 1:n.x.clim){
        beta.mat.prod[,i] <- beta.mat[i,x.mat.ind[,i]]
    }
    colnames(x.plot) <- x.clim.names
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
    palette("default")
    if(n.x.clim==1){
        x.plot.seq.backscaled.1 <- x.plot.seq*(max(x.clim)-min(x.clim))+min(x.clim)
        plot(x.plot.seq.backscaled.1,logityrange,type="l",xlab=x.labels[1],ylab="Probability of Occurrence",ylim=c(-0.1,1))
        points(x.clim,jitter(rep(-0.06,nrow(x.clim)),20),col=envelope.fit$y+1,cex=0.5)
    }else{
        if(missing(plot.vars)){
            plot.vars <- x.clim.names[1:2]
        }
        logityrange.tmp <- logityrange
        if(n.x.clim>2.5){
            x.plot.tmp <- x.plot
            x.clim.names.excluded <- x.clim.names[!(x.clim.names %in% plot.vars)]
            if(missing(fix.values)){
                fix.values <- round(rep(len/2,n.x.clim-2))
                names(fix.values) <- x.clim.names.excluded
            }else{
                for(i in 1:(n.x.clim-2)){
                    fix.values[x.clim.names.excluded[i]] <- round(((fix.values[x.clim.names.excluded[i]]-min(x.clim[,x.clim.names.excluded[i]]))/(max(x.clim[,x.clim.names.excluded[i]])-min(x.clim[,x.clim.names.excluded[i]])))*len)
                }
            }
            for(k in 1:(n.x.clim-2)){
                logityrange.tmp <- logityrange.tmp[x.plot[,x.clim.names.excluded[k]]==x.plot.seq[fix.values[x.clim.names.excluded[k]]]]
                x.plot.tmp <- x.plot.tmp[x.plot[,x.clim.names.excluded[k]]==x.plot.seq[fix.values[x.clim.names.excluded[k]]],]
            }
        }
        if(type=="persp"){
            op <- par()
            par(mar=c(1,0,1,0))
            temppersp <- persp(x.plot.seq,x.plot.seq,t(matrix(logityrange.tmp,ncol=len)),
                theta=20,phi=30,zlim=c(0,1),xlab=x.labels[plot.vars[1]],ylab=x.labels[plot.vars[2]],
                zlab="Probability of Occurrence")
            if(n.x.clim>2.5){
                if(close.points=="shade"){
                    points(trans3d(x.clim.std[,plot.vars[1]],x.clim.std[,plot.vars[2]],z=1,pmat=temppersp),
                        col=envelope.fit$y+1,
                        pch=1+15*(abs(rowSums(x.clim.std[,x.clim.names.excluded]-matrix(rep(fix.values[x.clim.names.excluded]/len,nrow(x.clim)),ncol=n.x.clim-2,byrow=TRUE)))<0.1))
                }else{
                    if(close.points=="plot"){
                        points(trans3d(x.clim.std[,plot.vars[1]],x.clim.std[,plot.vars[2]],z=1,pmat=temppersp),
                            cex=0+(abs(rowSums(x.clim.std[,x.clim.names.excluded]-matrix(rep(fix.values[x.clim.names.excluded]/len,nrow(x.clim)),ncol=n.x.clim-2,byrow=TRUE)))<0.1),
                            col=envelope.fit$y+1,pch=16)
                    }else{
                        points(trans3d(x.clim.std[,plot.vars[1]],x.clim.std[,plot.vars[2]],z=1,pmat=temppersp),
                            col=envelope.fit$y+1)
                    }
                }
            }else{
                points(trans3d(x.clim.std[,plot.vars[1]],x.clim.std[,plot.vars[2]],z=1,pmat=temppersp),
                    col=envelope.fit$y+1)
            }
            par(op)
        }else{ # type=="contour"
            if(missing(contour.levels)){
                contour.levels <- pretty(range(t(matrix(logityrange.tmp,ncol=len)), finite = TRUE), 10)
            }
            x.plot.seq.backscaled.1 <- x.plot.seq*(max(x.clim[,plot.vars[1]])-min(x.clim[,plot.vars[1]]))+min(x.clim[,plot.vars[1]])
            x.plot.seq.backscaled.2 <- x.plot.seq*(max(x.clim[,plot.vars[2]])-min(x.clim[,plot.vars[2]]))+min(x.clim[,plot.vars[2]])
            contour(x.plot.seq.backscaled.1,x.plot.seq.backscaled.2,t(matrix(logityrange.tmp,ncol=len)),
                xlab=x.labels[plot.vars[1]],ylab=x.labels[plot.vars[2]],levels=contour.levels)
            if(n.x.clim>2.5){
                if(close.points=="shade"){
                    points(x.clim[,plot.vars[1]],x.clim[,plot.vars[2]],col=envelope.fit$y+1,
                        pch=1+15*(abs(rowSums(x.clim.std[,x.clim.names.excluded]-matrix(rep(fix.values[x.clim.names.excluded]/len,nrow(x.clim)),ncol=n.x.clim-2,byrow=TRUE)))<0.1))
                }else{
                    if(close.points=="plot"){
                        points(x.clim[,plot.vars[1]],x.clim[,plot.vars[2]],col=envelope.fit$y+1,
                            pch=1+15*(abs(rowSums(x.clim.std[,x.clim.names.excluded]-matrix(rep(fix.values[x.clim.names.excluded]/len,nrow(x.clim)),ncol=n.x.clim-2,byrow=TRUE)))<0.1),
                            cex=0+(abs(rowSums(x.clim.std[,x.clim.names.excluded]-matrix(rep(fix.values[x.clim.names.excluded]/len,nrow(x.clim)),ncol=n.x.clim-2,byrow=TRUE)))<0.1))
                    }else{
                        points(x.clim[,plot.vars[1]],x.clim[,plot.vars[2]],col=envelope.fit$y+1)
                    }
                }
            }else{
                points(x.clim[,plot.vars[1]],x.clim[,plot.vars[2]],col=envelope.fit$y+1)
            }
        }
    }

}
