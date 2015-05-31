#' Plot maps showing predictions for a given set of climate data.
#'
#' \code{map.plot} plots the fitted maps of "climate suitability" for a given
#' set of climate variable values. Can also be use to plot maps of the presence/absence data.
#'
#' @param inputs If a numeric vector, then it contains the envelope parameters and has length 2p+p+2+p(p-1)/2.
#' If a logical vector, then it should have length n. The former is used for plotting shaded predictions for
#' the climate envelope models; the latter is used for making simple presence/absence plots. For presence/absence
#' plots, \code{x.clim.new} and \code{x.clim} are ignored.
#' @param x.clim.new A matrix, having the same columns as \code{x.clim}, and
#' containing the points in climate space at which to evaluate the
#' envelope and produce predictions in space.
#' @param x.clim The n by p matrix of climate covariates.
#' @param x.nonclim The n by p2 matrix of non-climate covariates (default
#' \code{NULL}).
#' @param pars.nonclim A vector of length p2 specifying the parameter estimates
#' to be used in prediction for the non-climate
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
map.plot <- function(inputs,x.clim.new,x.clim,x.nonclim=NULL,pars.nonclim=NULL,
    coordinates,species.name="",scenario.name="",save.PDF=FALSE,file.name){
    plot.type <- "Prediction"
    if(is.logical(inputs)){
        plot.type <- "Presences"
        scenario.name <- "PA"
    }
    if(save.PDF & missing(file.name)){
        file.name <- paste("Map",species.name,scenario.name,".pdf",sep="")
    }
    if(plot.type=="Prediction"){
        if(missing(x.clim)){
            x.clim <- x.clim.new
        }
        n.x.clim <- ncol(x.clim.new)
        # Now standardise the climate variables by mapping onto [0,1]
        x.clim.std <- x.clim.new
        x.clim.min <- apply(x.clim,2,min)
        x.clim.max <- apply(x.clim,2,max)
        for(i in 1:n.x.clim){
            x.clim.std[,i] <- (x.clim.new[,i]-x.clim.min[i])/(x.clim.max[i]-x.clim.min[i])
        }
        x.envelope <- env.fn(inputs,x.clim.std)$x.envelope
        # Check the columns of x.nonclim are in the same order as the coefficients
        if(!is.null(x.nonclim)){
            x.nonclim <- as.matrix(x.nonclim[names(pars.nonclim)])
            nonlinpart <- x.nonclim%*%pars.nonclim
        }else{
            nonlinpart <- 0
        }
        linpred <- x.envelope+nonlinpart
        plot.p <- 1/(1+exp(linpred))        plot.p <- round(plot.p*1000)
        palette(grey(seq(0,1,len=1001)))
    }else{
        plot.p <- 15*inputs+1
    }
    if(save.PDF){pdf(file.name)}
    op <- par()
    par(mar=c(1,1,1,1)+0.1)
    plot(mapproject(x=coordinates$long,y=coordinates$lat,proj='azequalarea'),
        pch=16,xaxt="n",yaxt="n",col=plot.p,
        xlab="",ylab="")
    if(plot.type=="Prediction"){
        text(-0.226,-0.23,paste(scenario.name,", ",species.name,sep=""),pos=4,cex=2)
        legend(-0.21,-0.25,legend=c("Probability 0.1","Probability 0.5","Probability 0.9"),pch=c(16,16,16),
            col=c(900,500,100))
    }else{
        text(-0.226,-0.23,species.name,pos=4,cex=2)
        legend(-0.2,-0.25,legend=c("Recorded Presence","Recorded Absence"),pch=c(16,1))
    }
    par(op)
    palette("default")
    if(save.PDF){dev.off()}
}
