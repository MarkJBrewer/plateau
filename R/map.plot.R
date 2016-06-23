#' Plot maps showing predictions for a given set of climate data.
#'
#' \code{map.plot} plots the fitted maps of "climate suitability" for a given
#' set of climate variable values. Can also be use to plot maps of the presence/absence data.
#'
#' @param inputs If \code{plot.type} has value "Prediction", then \code{inputs} contains a fitted model
#' object as output by \code{fit.glm.env} or \code{fit.bugs.env}. If \code{plot.type} has value "Presence",
#' then \code{inputs} should be a numeric vector of length n and containing 0's and 1's (or FALSE and TRUE) marking absence
#' and presence respectively. The "Prediction" form is used for plotting shaded predictions for
#' fitted climate envelope models; the "Presence" form is used for making simple presence/absence plots.
#' For presence/absence plots, \code{x.clim.new}, and other similar arguments are ignored.
#' @param plot.type A string to indicate what kind of map to plot; takes values "Prediction" (the
#' default) or "Presence". (In fact, any string not identical to "Prediction" will lead to a Presence plot.)
#' @param x.clim.new A matrix, having the same columns as \code{inputs$x.clim}, and
#' containing the points in climate space at which to evaluate the
#' envelope and produce predictions in space.
#' @param x.nonclim.new The n by p2 matrix of non-climate covariates for prediction, only to
#' be specified if a corresponding \code{x.nonclim} element exists in \code{inputs} (default
#' \code{NULL}).
#' @param x.factor.new The n by p3 matrix of non-climate factors for prediction, only to
#' be specified if a corresponding \code{x.factor} element exists in \code{inputs} (default
#' \code{NULL}).
#' @param coordinates An n by 2 matrix or list or data frame, having elements
#'  (or columns) \code{long} and \code{lat} for longitude and latitude
#' respectively.
#' @param species.name The pretty name of the species for plotting.
#' @param scenario.name The name of the scenario (and year) for plotting.
#' @param save.PDF Logical to indicate whether to save as PDF (default
#' \code{FALSE}).
#' @param file.name If \code{save.PDF} is \code{TRUE}, the name of the pdf
#' file; if missing, this is created using the species and scenario names
#' (if available).
#'
#' @return Map of predictions of "climate suitability".
#' @seealso \code{link{envelope.plot}}
#' @export
map.plot <- function(inputs,plot.type="Prediction",x.clim.new,x.nonclim.new=NULL,
    x.factor.new=NULL,coordinates,species.name="",scenario.name="",save.PDF=FALSE,file.name){
    if(plot.type=="Prediction"){
      env.pars <- inputs$par
    }
    if(save.PDF & missing(file.name)){
        file.name <- paste("Map",species.name,scenario.name,".pdf",sep="")
    }
    if(plot.type=="Prediction"){
        x.clim <- inputs$x.clim
        n.x.clim <- ncol(x.clim.new)
        # Now standardise the climate variables by mapping onto [0,1]
        x.clim.std <- x.clim.new
        x.clim.min <- apply(x.clim,2,min)
        x.clim.max <- apply(x.clim,2,max)
        for(i in 1:n.x.clim){
            x.clim.std[,i] <- (x.clim.new[,i]-x.clim.min[i])/(x.clim.max[i]-x.clim.min[i])
        }
        x.envelope <- env.fn(env.pars,x.clim.std)$x.envelope
        nonlinpart <- 0
        if("x.nonclim" %in% names(inputs)){
            pars.nonclim <- coef(inputs$nonclim.glm)[names(x.nonclim)]
            x.nonclim <- inputs$x.nonclim
            x.nonclim.means <- colMeans(x.nonclim)
            x.nonclim.new.centred <- matrix(rep(x.nonclim.means,nrow(x.nonclim.new)),byrow=TRUE,nrow=nrow(x.nonclim.new))
            nonlinpart <- nonlinpart+x.nonclim.new.centred%*%pars.nonclim
        }
        if("x.factor" %in% names(inputs)){
            x.factor <- inputs$x.factor
            if(length(x.factor)==1){
              temp.nlevels <- nlevels(x.factor)
              x.factor <- (contr.sum(temp.nlevels)+contr.helmert(temp.nlevels))[x.factor,]
              x.factor.new <- (contr.sum(temp.nlevels)+contr.helmert(temp.nlevels))[x.factor.new,]
            }else{
              for(i in 1:ncol(x.factor)){
                temp.nlevels <- nlevels(x.factor[,i])
                x.factor[,i] <- (contr.sum(temp.nlevels)+contr.helmert(temp.nlevels))[x.factor[,i],]
                x.factor.new[,i] <- (contr.sum(temp.nlevels)+contr.helmert(temp.nlevels))[x.factor.new[,i],]
              }
            }
            pars.factor <- coef(inputs$nonclim.glm)[-names(x.nonclim)]
            nonlinpart <- nonlinpart+x.factor.new%*%pars.factor
        }  
        linpred <- x.envelope+nonlinpart
        plot.p <- 1/(1+exp(linpred))
        plot.p <- round(plot.p*1000)
        palette(grey(seq(0,1,len=1001)))
    }else{
        if(scenario.name==""){
            scenario.name <- "PA"
        }
        plot.p <- 15*inputs+1
    }
    if(save.PDF){pdf(file.name)}
    op <- par(no.readonly = TRUE)
    par(mar=c(1,1,1,1)+0.1)
    if(plot.type=="Prediction"){
        plot(mapproject(x=coordinates$long,y=coordinates$lat,projection='azequalarea'),
            pch=16,xaxt="n",yaxt="n",col=plot.p,xlab="",ylab="")
        text(-0.226,-0.23,paste(scenario.name,", ",species.name,sep=""),pos=4,cex=2)
        legend(-0.21,-0.25,legend=c("Probability 0.1","Probability 0.5","Probability 0.9"),pch=c(16,16,16),
            col=c(900,500,100))
    }else{
        plot(mapproject(x=coordinates$long,y=coordinates$lat,projection='azequalarea'),
            pch=plot.p,xaxt="n",yaxt="n",xlab="",ylab="")
        text(-0.226,-0.23,species.name,pos=4,cex=2)
        legend(-0.2,-0.25,legend=c("Recorded Presence","Recorded Absence"),pch=c(16,1))
    }
    par(op)
    palette("default")
    if(save.PDF){dev.off()}
}
