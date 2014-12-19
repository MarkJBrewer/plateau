#' Plot maps showing predictions for a given set of climate data.
#'
#' \code{map.plot} plots the fitted maps of "climate suitability" for a given
#' set of climate variable values.
#'
#' @param pars The vector of envelope parameters, length 2p+p+2+p(p-1)/2.
#' @param x.clim A matrix, the points in climate space at which to evaluate the
#'  envelope and produce predictions in space.
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
map.plot <- function(pars,x.clim.new,x.clim,coordinates,species.name="",scenario.name="",
    save.PDF=FALSE,file.name){
    if(save.PDF & missing(file.name)){
        file.name <- paste("Map",species.name,scenario.name,".pdf",sep="")
    }
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
    x.envelope <- env.fn(pars,x.clim.std)$x.envelope
    plot.p <- 1/(1+exp(x.envelope))
    if(save.PDF){pdf(file.name)}
    op <- par()
    plot.p <- round(plot.p*1000)
    palette(grey(seq(0,1,len=1001)))
    par(mar=c(1,1,1,1)+0.1)
    plot(mapproject(x=coordinates$long,y=coordinates$lat,proj='azequalarea'),
        pch=16,xaxt="n",yaxt="n",col=plot.p,
        xlab="",ylab="")
    text(-0.226,-0.23,paste(scenario.name,", ",species.name,sep=""),pos=4,cex=2)
    legend(-0.21,-0.25,legend=c("Probability 0.1","Probability 0.5","Probability 0.9"),pch=c(16,16,16),
        col=c(900,500,100))
    par(op)
    palette("default")
    if(save.PDF){dev.off()}
}
