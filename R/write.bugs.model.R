#' Write WinBUGS code to a temporary file
#'
#' @param n.x.clim The number of climate covariates (n); needed as the WinBUGS
#' code is slightly different for n=1.
#' @param n.x.nonclim The number of non-climate covariates; needed to know
#' whether or not to include a separate loop for such covariates.
#' @param n.x.factor The number of non-climate factors; needed to know
#' whether or not to include a separate loop for such factors
#' @param not.spatial Logical scalar specifying whether the model should include
#' a spatial intrinsic CAR component. Defaults to FALSE, i.e. there \emph{is} a
#' spatial component.
#' @param working.directory Working directory as supplied by the function \code{fit.bugs.env},
#' either from the user or auto-generated.
#' @export
write.bugs.model <- function(n.x.clim,n.x.nonclim=0,n.x.factor=0,not.spatial=FALSE,
                                     working.directory=NULL){
  TEMP_BUG <- "
  model{
    for(i in 1:N) {
      y[i] ~ dbern(p[i])
      for(e in 1:NEnv){
        x.clim.C[i,e] <- x.clim[i,e]-ax[e]
        x.clim.beta[i,e,1] <- -beta[e,1]*x.clim.C[i,e]
        x.clim.beta[i,e,2] <- beta[e,which.beta[e]]*x.clim.C[i,e]
        x.clim.Calc2[i,e] <- max(x.clim.beta[i,e,1],x.clim.beta[i,e,2])
        x.clim.Calc[i,e] <- x.clim.Calc2[i,e]*abs(x.clim.C[i,e])
      }
  
      $TEMP_CROSSSUM
  
      #TEMP_FACTOR1
  
      #TEMP_NOC1
  
      logit(p[i]) <- x.clim.Total[i] #+ u[u.clique.start[clique[i]]+clique.i[i]-1]
    }
    beta0.prec <- 1/beta0.var
    beta0.diff ~ dnorm(beta0.mu,beta0.prec)
    az ~ dnorm(0,0.1)
    beta0 <- az - exp(beta0.diff)
  
    #TEMP_SPATIAL
  
    tau <- car.tau
  
    #TEMP_FACTOR2
  
    #TEMP_NOC2
  
    $ALL_FOR
  
    $TEMP_GAMMA
  }
  "
  
  TEMP_FOR <-"
  for(e in $VAR_NUM : $VAR_NUM  ){
    ax.prec[e] <- 1/ax.var[e]
    ax[e] ~ dnorm(ax.mu[e],ax.prec[e])  I(-1,2)
    beta.prec[e,1] <- 1/beta.var[e,1]
    beta.prec[e,2] <- 1/beta.var[e,2]
  
    $BETA_1
  
    $BETA_2
  }"

  TEMP_CROSSSUM2 <-"
  for(e1 in 1:(NEnv-1)){
    for(e2 in 1:e1){
      x.clim.CalcCross[i,e1,e2] <- 0
    }
    for(e2 in (e1+1):NEnv){
      x.clim.CalcCross[i,e1,e2] <- gamma[e1,e2]*x.clim.C[i,e1]*x.clim.C[i,e2]
    }
  }
  for(e in 1:NEnv){
    x.clim.CalcCrossColSum[i,e] <- sum(x.clim.CalcCross[i,,e])
  }
  x.clim.CalcSum[i] <- sum(x.clim.Calc[i,1:NEnv])
  x.clim.CalcCrossSum[i] <- sum(x.clim.CalcCrossColSum[i,])
  x.clim.Total.tmp[i] <- az - sqrt(x.clim.CalcSum[i] + x.clim.CalcCrossSum[i])
  x.clim.Total[i] <- min(x.clim.Total.tmp[i], beta0)
  "
  
  TEMP_CROSSSUM1 <- "
  x.clim.CalcSum[i] <- sum(x.clim.Calc[i,1:NEnv])
  x.clim.Total.tmp[i] <- az - sqrt(x.clim.CalcSum[i])
  x.clim.Total[i] <- min(x.clim.Total.tmp[i], beta0)
  "
  
  TEMP_FACTOR1 <-"
  for(e in 1:Nfac){
    x.fac.Calc[i,e] <- facbeta[e,x.factor[i,e]]
  }
  x.fac.Total[i] <- sum(x.fac.Calc[i,1:Nfac])
  "
  
  TEMP_FACTOR2 <-"
  for(e in 1:Nfac){
    facbeta[e,1] <- -sum(facbeta[e,2:x.factor.lengths[e]])
    for(e1 in 2:x.factor.lengths[e]){
      facbeta[e,e1] ~ dnorm(0,0.001)
    }
  }"

  TEMP_NOC1 <- "
  for(e in 1:NnonEnv){
    x.nonclim.Calc[i,e] <- nonbeta[e]*(x.nonclim[i,e]-mean(x.nonclim[,e]))
  }
  x.nonclim.Total[i] <- sum(x.nonclim.Calc[i,1:NnonEnv])
  "
  
  TEMP_NOC2 <- "
  for(e in 1:NnonEnv){
    nonbeta[e] ~ dnorm(0,0.001)
  }"
  
  TEMP_GAMMA <- "
  for(e1 in 1:(NEnv-1)){
    gamma.temp[e1,1] <- 0
    for(e2 in (e1+1):NEnv){
      gamma.temp[e1,e2] ~ dunif(0,1)
      min.gamma.1[e1,e2] <- beta[e1,1]*beta[e2,1]
      min.gamma.2[e1,e2] <- beta[e1,which.beta[e1]]*beta[e2,which.beta[e2]]
      max.gamma.1[e1,e2] <- beta[e1,1]*beta[e2,which.beta[e2]]
      max.gamma.2[e1,e2] <- beta[e1,which.beta[e1]]*beta[e2,1]
      min.gamma[e1,e2] <- -(2/(NEnv-1))*sqrt(min(min.gamma.1[e1,e2],min.gamma.2[e1,e2]))
      max.gamma[e1,e2] <- (2/(NEnv-1))*sqrt(min(max.gamma.1[e1,e2],max.gamma.2[e1,e2]))
      gamma[e1,e2] <- min.gamma[e1,e2]+gamma.temp[e1,e2]*(max.gamma[e1,e2]-min.gamma[e1,e2])
    }
  }
  for(e2 in 1:NEnv){
    gamma.temp[NEnv,e2] <- 0
  }"

  TEMP_SPATIAL <- "
  # the following few lines are for spatial autocorrelation module
  for(i in 1:NonSingletonClique){
    clique.length[nonsingleton.clique.list[i]] <- u.clique.end[nonsingleton.clique.list[i]] - u.clique.start[nonsingleton.clique.list[i]] + 1
    adj.clique.length[nonsingleton.clique.list[i]] <- adj.clique.end[nonsingleton.clique.list[i]] - adj.clique.start[nonsingleton.clique.list[i]] + 1
    u[u.clique.start[nonsingleton.clique.list[i]]:u.clique.end[nonsingleton.clique.list[i]]] ~ car.normal(adj[adj.clique.start[nonsingleton.clique.list[i]]:adj.clique.end[nonsingleton.clique.list[i]]], weights[1:adj.clique.length[nonsingleton.clique.list[i]]], num[u.clique.start[nonsingleton.clique.list[i]]:u.clique.end[nonsingleton.clique.list[i]]], tau)
  }
  for(k in 1:adj.clique.end[NClique]) {
    weights[k] <- 1
  }
  "
  
  ALL_FOR <-""		
  for(j in 1:n.x.clim){  # loop for number of variables
    BETA_1 <- "beta[e,1] ~ dnorm(beta.mu[e,1],beta.prec[e,1])  I(0,)"
    BETA_2 <- "beta[e,2] ~ dnorm(beta.mu[e,2],beta.prec[e,2])  I(0,)"
    # revise FOR STATEMENT
    ED_FOR <- gsub("$VAR_NUM",j,TEMP_FOR,fixed=TRUE)
    ED_FOR <- gsub("$BETA_1",BETA_1,ED_FOR,fixed=TRUE)
    ED_FOR <- gsub("$BETA_2",BETA_2,ED_FOR,fixed=TRUE)
    ALL_FOR <- paste(ALL_FOR, ED_FOR,sep="\n")
  }
  #cat(ED_FOR)
  #cat(ALL_FOR)		
  ED_BUG <- gsub("$ALL_FOR",ALL_FOR,TEMP_BUG,fixed=TRUE)
  #cat(ED_BUG)
  
  # check 1 or >=2 climatic variables, if >=2, then include variable interactions
  #CROSSSUM
  #GAMMA
  if (n.x.clim==1){
    ED_BUG <- gsub("$TEMP_CROSSSUM",TEMP_CROSSSUM1,ED_BUG,fixed=TRUE)
    ED_BUG <- gsub("$TEMP_GAMMA","",ED_BUG,fixed=TRUE)
  } else {
    ED_BUG <- gsub("$TEMP_CROSSSUM",TEMP_CROSSSUM2,ED_BUG,fixed=TRUE)
    ED_BUG <- gsub("$TEMP_GAMMA",TEMP_GAMMA,ED_BUG,fixed=TRUE)
  }
  
  # check nonclim & factor
  if (n.x.nonclim > 0.5){
    ED_BUG <- gsub("#TEMP_NOC1",TEMP_NOC1,ED_BUG,fixed=TRUE)
    ED_BUG <- gsub("#TEMP_NOC2",TEMP_NOC2,ED_BUG,fixed=TRUE)
    ED_BUG <- gsub("logit(p[i]) <-",
                   "logit(p[i]) <- x.nonclim.Total[i] + ",
                   ED_BUG,fixed=TRUE)
  }
  if (n.x.factor > 0.5){
    ED_BUG <- gsub("#TEMP_FACTOR1",TEMP_FACTOR1,ED_BUG,fixed=TRUE)
    ED_BUG <- gsub("#TEMP_FACTOR2",TEMP_FACTOR2,ED_BUG,fixed=TRUE)
    ED_BUG <- gsub("logit(p[i]) <-",
                   "logit(p[i]) <- x.fac.Total[i] + ",
                   ED_BUG,fixed=TRUE)
  }
  
  # check if spatial autocorrelation is used
  if(!not.spatial){
    ED_BUG <- gsub("#TEMP_SPATIAL",TEMP_SPATIAL,ED_BUG,fixed=TRUE)
  }
  
  output_file_path <- paste(working.directory,"/WINBUGS_code.txt",sep="")
  writeLines(ED_BUG,output_file_path)
}

