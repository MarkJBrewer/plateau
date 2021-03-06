% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit.bugs.nospatial.env.R
\name{fit.bugs.nospatial.env}
\alias{fit.bugs.nospatial.env}
\title{Fit a spatial, Bayesian, plateau envelope model via MCMC in WinBUGS}
\usage{
fit.bugs.nospatial.env(data, y, x.clim, x.nonclim = NULL, x.factor = NULL,
  prior.ax, prior.beta, prior.beta0.difference, constrain.beta,
  initial.pars.input, informative.priors = list(beta = FALSE, beta0 = FALSE,
  ax = FALSE), burnin = 5000, post.burnin = 1000, chains = 2, thin = 1,
  working.directory = NULL, silent = TRUE,
  bugs.directory = "C:/Program Files (x86)/WinBUGS14/",
  WinBUGS.debug = FALSE, WinBUGS.code = NULL,
  WinBUGS.code.location = NULL, no.starting.value = NULL,
  estimate.p = FALSE)
}
\arguments{
\item{data}{The data frame (with n rows) containing all the variables for analysis.}

\item{y}{A string denoting the binary response variable (taking values 0 or 1 for absence and
presence respectively); must correspond to a column name in the data frame
specied at \code{data}.}

\item{x.clim}{A vector (length p) of strings denoting which columns in the supplied data
frame correspond to the climate covariates; must correspond to column names in the data frame
specied at \code{data}.}

\item{x.nonclim}{A vector (length p2) of strings denoting which columns in the supplied data
frame correspond to the  non-climate covariates; must correspond to column names in the data frame
specied at \code{data}.}

\item{x.factor}{A vector (length p3) of strings denoting which columns in the supplied data
frame correspond to the non-climate factors; must correspond to column names in the data frame
specied at \code{data}.

The next four inputs are used to set identifiability constraints on the
modelling.}

\item{prior.ax}{A list of up to two p-vector objects, \code{mean} and
\code{var}, containing the means and variances of (informative) priors
for the apex x co-ordinates ax; if either \code{mean} or \code{var} are
omitted, the means and variances are set to 0.5 and 1 respectively.}

\item{prior.beta}{A list of two p by 2 matrix objects, \code{mean} and
\code{var}, containing the means and variances of (informative) priors
for the slopes beta; if omitted, the means and variances are
set to 100 and 10000 respectively.}

\item{prior.beta0.difference}{A list of two scalar objects, \code{mean} and
\code{var}, containing the mean and variance for the log-difference
between the apex az and plateau level beta0; if omitted,
the mean and variance are set to 1 and 10 respectively.}

\item{constrain.beta}{A p by 2 matrix of logicals in order to indicate
which beta parameters should be constrained to not vary too
much from its pair - set at most only one of these to \code{TRUE}
in each row. A value of \code{TRUE} in column 1 suggests that
there are not enough data to estimate the "up" (left-hand)
part of the envelope, and in column 2 correspondingly for
the "down" (right-hand) part of the envelope. If
\code{initial.pars.input} is not set, the code will work out which
betas should be constrained.}

\item{initial.pars.input}{Vector of length 2p+p+2+p(p-1)/2 containing
starting values for each parameter; if missing, the code
works out its own starting value. If \code{x.factor} is non-null, there will be a total of
the number of total levels across all the factors minus p3 extra parameters.}

\item{informative.priors}{List of logical scalars for which informative
priors should be used (from: beta, beta0, ax);
this option uses the function \code{generate.initial.values} to
generate the means, and is an alternative to specifying the
priors yourself using \code{prior.ax} above etc. Default has all
options set to \code{FALSE}.}

\item{burnin}{Scalar specifying the number of "burn-in" iterations to be
discarded (default 5000).}

\item{post.burnin}{Scalar specifying the number of subsequent iterations to
be retained (default 1000).}

\item{chains}{Scalar, number of parallel chains to run (default 2).}

\item{thin}{Scalar, the thinning to apply to the MCMC iterations (default 1).}

\item{working.directory}{String containing the location of the WinBUGS code
file; default NULL, in which case a temporary folder is used.}

\item{silent}{Logical flag denoting whether the function runs silently or
not. Default is \code{TRUE}.}

\item{bugs.directory}{String containing location of WinBUGS installation;
default "C:/Program Files (x86)/WinBUGS14/" is for Windows
64-bit machines.}

\item{WinBUGS.debug}{Logical flag as to whether to close WinBUGS after
running (default \code{FALSE} implies WinBUGS is closed).}

\item{WinBUGS.code}{You can supply your own code file, especially useful if
you want to use informative priors for external information. The default
value is \code{NULL}, and the code file must be in the folder as specified by
\code{WinBUGS.code.location}. The effect of a \code{NULL} value that the
function \code{write.bugs.model} is called, and the WinBUGS code file is
generated automatically.}

\item{WinBUGS.code.location}{If \code{WinBUGS.code} is not \code{NULL} then
this folder is examined for the code file specified.  The default value is
\code{NULL}, whereby the code looks for the file in the current directory.}

\item{no.starting.value}{A list of strings denoting the objects we should
not initialise, usually because we set them or calculate them in a bespoke
WinBUGS code file. These objects will be set to \code{NA} by this function;
to use, create objects such as
                          \code{list("beta[2,2]","ax[1]")}
for example.}

\item{estimate.p}{Logical flag specifying whether or not to retain samples
for the posterior probabilities of presence for each cell,
default \code{FALSE}; be aware using \code{TRUE} can result in slow
interaction between R and WinBUGS.}
}
\value{
A list object as returned by the \code{bugs()} function in
\pkg{R2WinBUGS}; see the help file for that function for further details.

The list is augmented by the response \code{y} and the climate variables
\code{x.clim}, primarily to aid the plotting functions, and by
the vector which.beta, of length p, which contains 1's and
2's to signify which set of beta parameters to use; see
the description of beta below.

The variables in WinBUGS are (i indexes the climate variable):
\describe{
\item{\code{beta[which.beta[i],i,1/2]}}{\code{which.beta} (see above)
                   specifies which of two sets of \code{beta} parameters to use.
                   There are two \code{betas} in order to facilitate constraining
                   either the \code{beta[i,1]} (up-slope) or \code{beta[i,2]} (down-slope)
                   is to be constrained to be close to its opposite. If
                   \code{beta[i,2]} is constrained, we'll use the first set (i.e.
                   \code{which.beta[i]==1}) as this models \code{beta[i,1]} directly and
                   the multiplicative difference between that and \code{beta[i,2]};
                   if \code{beta[i,1]} is constrained (\code{which.beta[i]==2}) then we
                   model \code{beta[i,2]} directly.}
\item{\code{ax}}{The x-coordinate of the apex for each climate variable.}
\item{\code{beta0}}{Scalar suggesting where the top-slicing should be applied.}
\item{\code{az}}{Scalar coordinate on the logit scale of the response axis apex.}
\item{\code{gamma}}{(Upper triangle) matrix of pairwise interaction parameters
                   between the climate variables; these are constrained
                   relative to the \code{betas} in order not to break the geometric
                   formula for a cone.}
}
}
\description{
\code{fit.bugs.nospatial.env} fits a spatial GLM with plateau envelope on the climate
covariates via MCMC in WinBUGS.
}

