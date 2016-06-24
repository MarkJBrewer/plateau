#' Write WinBUGS code to a temporary file
#'
#' @param n.x.clim The number of climate covariates (n); needed as the WinBUGS
#' code is slightly different for n=1.
#' @param n.x.nonclim The number of non-climate covariates; needed to know
#' whether or not to include a separate loop for such covariates.
#' @param n.x.factor The number of non-climate factors; needed to know
#' whether or not to include a separate loop for such factors
#' @export
write.bugs.model <- function(n.x.clim,n.x.nonclim=0,n.x.factor=0){

    # WinBUGS code as a function
    NClique <- car.tau <- nonbeta <- beta0 <- az <- ax.var <- beta0.diff <- beta0.var <- beta0.mu <- N <- NEnv <- NnonEnv <- Nfac <- NonSingletonClique <- 1
    u.clique.start <- u.clique.end <- nonsingleton.clique.list <- x.clim.Total <- x.clim.Total.tmp <- ax <- which.beta <- p <- y <- x.factor.lengths <- numeric(2)
    adj.clique.start <- adj.clique.end <- adj <- clique <- clique.u <- numeric(2)
    x.clim <- x.nonclim <- x.factor <- min.gamma.1 <- min.gamma.2 <- max.gamma.1 <- max.gamma.2 <- array(NA,dim=c(2,2))
    min.gamma <- max.gamma <- gamma.temp <- x.clim.CalcCrossColSum <- gamma <- x.clim.Calc <- x.clim.Calc2 <- x.clim.C <- beta <- facbeta <- x.fac.Calc <- x.fac.Total <- array(NA,dim=c(2,2))
    x.clim.CalcCross <- x.clim.beta <- array(NA,dim=c(2,2,2))
    u <- clique.i <- numeric(2)
    beta.var <- array(NA,dim=c(2,2))
    `logit<-` <- function(x){x}
    
    if(n.x.factor==0 && n.x.nonclim==0 && n.x.clim>1.5){
        envmodel <- function(){
            pi <- 3.141592
            piover2 <- pi/2

            for(i in 1:N) {
                y[i] ~ dbern(p[i])
                for(e in 1:NEnv){
                    x.clim.C[i,e] <- x.clim[i,e]-ax[e]
                    x.clim.beta[i,e,1] <- -beta[e,1]*x.clim.C[i,e]
                    x.clim.beta[i,e,2] <- beta[e,which.beta[e]]*x.clim.C[i,e]
                    x.clim.Calc2[i,e] <- max(x.clim.beta[i,e,1],x.clim.beta[i,e,2])
                    x.clim.Calc[i,e] <- x.clim.Calc2[i,e]*abs(x.clim.C[i,e])
                }
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
        	   logit(p[i]) <- x.clim.Total[i] + u[u.clique.start[clique[i]]+clique.i[i]-1]
            }
            beta0.prec <- 1/beta0.var
            beta0.diff ~ dnorm(beta0.mu,beta0.prec)
            az ~ dnorm(0,0.1)
            beta0 <- az - exp(beta0.diff)
            for(i in 1:NonSingletonClique){
              clique.length[nonsingleton.clique.list[i]] <- u.clique.end[nonsingleton.clique.list[i]] - u.clique.start[nonsingleton.clique.list[i]] + 1
              adj.clique.length[nonsingleton.clique.list[i]] <- adj.clique.end[nonsingleton.clique.list[i]] - adj.clique.start[nonsingleton.clique.list[i]] + 1
              u[u.clique.start[nonsingleton.clique.list[i]]:u.clique.end[nonsingleton.clique.list[i]]] ~ car.normal(adj[adj.clique.start[nonsingleton.clique.list[i]]:adj.clique.end[nonsingleton.clique.list[i]]], weights[1:adj.clique.length[nonsingleton.clique.list[i]]], num[u.clique.start[nonsingleton.clique.list[i]]:u.clique.end[nonsingleton.clique.list[i]]], tau)
            }
            for(k in 1:adj.clique.end[NClique]) {
        	   weights[k] <- 1
            }
            tau <- car.tau

            for(e in 1:NEnv){
                ax.prec[e] <- 1/ax.var[e]
                ax[e] ~ dnorm(ax.mu[e],ax.prec[e])%_%I(-1,2)
                beta.prec[e,1] <- 1/beta.var[e,1]
                beta.prec[e,2] <- 1/beta.var[e,2]
                beta[e,1] ~ dnorm(beta.mu[e,1],beta.prec[e,1])%_%I(0,)
                beta[e,2] ~ dnorm(beta.mu[e,2],beta.prec[e,2])%_%I(0,)
            }

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
            }

        }
    }
    if(n.x.factor==0 && n.x.nonclim==0 && n.x.clim==1){
        envmodel <- function(){
            pi <- 3.141592
            piover2 <- pi/2
            for(i in 1:N) {
                y[i] ~ dbern(p[i])
                for(e in 1:NEnv){
                    x.clim.C[i,e] <- x.clim[i,e]-ax
                    x.clim.beta[i,e,1] <- -beta[e,1]*x.clim.C[i,e]
                    x.clim.beta[i,e,2] <- beta[e,which.beta]*x.clim.C[i,e]
                    x.clim.Calc2[i,e] <- max(x.clim.beta[i,e,1],x.clim.beta[i,e,2])
                    x.clim.Calc[i,e] <- x.clim.Calc2[i,e]*abs(x.clim.C[i,e])
                }
                x.clim.CalcSum[i] <- sum(x.clim.Calc[i,1:NEnv])
                x.clim.Total.tmp[i] <- az - sqrt(x.clim.CalcSum[i])
                x.clim.Total[i] <- min(x.clim.Total.tmp[i], beta0)
                logit(p[i]) <- x.clim.Total[i] + u[u.clique.start[clique[i]]+clique.i[i]-1]
            }
            beta0.prec <- 1/beta0.var
            beta0.diff ~ dnorm(beta0.mu,beta0.prec)
            az ~ dnorm(0,0.1)
            beta0 <- az - exp(beta0.diff)
            for(i in 1:NonSingletonClique){
              clique.length[nonsingleton.clique.list[i]] <- u.clique.end[nonsingleton.clique.list[i]] - u.clique.start[nonsingleton.clique.list[i]] + 1
              adj.clique.length[nonsingleton.clique.list[i]] <- adj.clique.end[nonsingleton.clique.list[i]] - adj.clique.start[nonsingleton.clique.list[i]] + 1
              u[u.clique.start[nonsingleton.clique.list[i]]:u.clique.end[nonsingleton.clique.list[i]]] ~ car.normal(adj[adj.clique.start[nonsingleton.clique.list[i]]:adj.clique.end[nonsingleton.clique.list[i]]], weights[1:adj.clique.length[nonsingleton.clique.list[i]]], num[u.clique.start[nonsingleton.clique.list[i]]:u.clique.end[nonsingleton.clique.list[i]]], tau)            }
            for(k in 1:adj.clique.end[NClique]) {
                weights[k] <- 1
            }
            tau <- car.tau
            ax.prec <- 1/ax.var
            ax ~ dnorm(ax.mu,ax.prec)%_%I(-1,2)
            for(e in 1:NEnv){
                beta.prec[e,1] <- 1/beta.var[e,1]
                beta.prec[e,2] <- 1/beta.var[e,2]
                beta[e,1] ~ dnorm(beta.mu[e,1],beta.prec[e,1])%_%I(0,)
                beta[e,2] ~ dnorm(beta.mu[e,2],beta.prec[e,2])%_%I(0,)
            }
        }
    }
    if(n.x.factor==0 && n.x.nonclim>0.5 && n.x.clim>1.5){
        envmodel <- function(){
            pi <- 3.141592
            piover2 <- pi/2

            for(i in 1:N) {
                y[i] ~ dbern(p[i])
                for(e in 1:NEnv){
                    x.clim.C[i,e] <- x.clim[i,e]-ax[e]
                    x.clim.beta[i,e,1] <- -beta[e,1]*x.clim.C[i,e]
                    x.clim.beta[i,e,2] <- beta[e,which.beta[e]]*x.clim.C[i,e]
                    x.clim.Calc2[i,e] <- max(x.clim.beta[i,e,1],x.clim.beta[i,e,2])
                    x.clim.Calc[i,e] <- x.clim.Calc2[i,e]*abs(x.clim.C[i,e])
                }
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
                for(e in 1:NnonEnv){
                    x.nonclim.Calc[i,e] <- nonbeta[e]*(x.nonclim[i,e]-mean(x.nonclim[,e]))
                }
                x.nonclim.Total[i] <- sum(x.nonclim.Calc[i,1:NnonEnv])
                logit(p[i]) <- x.clim.Total[i] + x.nonclim.Total[i] + u[u.clique.start[clique[i]]+clique.i[i]-1]
            }
            beta0.prec <- 1/beta0.var
            beta0.diff ~ dnorm(beta0.mu,beta0.prec)
            az ~ dnorm(0,0.1)
            beta0 <- az - exp(beta0.diff)
            for(i in 1:NonSingletonClique){
              clique.length[nonsingleton.clique.list[i]] <- u.clique.end[nonsingleton.clique.list[i]] - u.clique.start[nonsingleton.clique.list[i]] + 1
              adj.clique.length[nonsingleton.clique.list[i]] <- adj.clique.end[nonsingleton.clique.list[i]] - adj.clique.start[nonsingleton.clique.list[i]] + 1
              u[u.clique.start[nonsingleton.clique.list[i]]:u.clique.end[nonsingleton.clique.list[i]]] ~ car.normal(adj[adj.clique.start[nonsingleton.clique.list[i]]:adj.clique.end[nonsingleton.clique.list[i]]], weights[1:adj.clique.length[nonsingleton.clique.list[i]]], num[u.clique.start[nonsingleton.clique.list[i]]:u.clique.end[nonsingleton.clique.list[i]]], tau)
            }
            for(k in 1:adj.clique.end[NClique]) {
        	   weights[k] <- 1
            }
            tau <- car.tau

            for(e in 1:NEnv){
                ax.prec[e] <- 1/ax.var[e]
                ax[e] ~ dnorm(ax.mu[e],ax.prec[e])%_%I(-1,2)
                beta.prec[e,1] <- 1/beta.var[e,1]
                beta.prec[e,2] <- 1/beta.var[e,2]
                beta[e,1] ~ dnorm(beta.mu[e,1],beta.prec[e,1])%_%I(0,)
                beta[e,2] ~ dnorm(beta.mu[e,2],beta.prec[e,2])%_%I(0,)
            }
            for(e in 1:NnonEnv){
                nonbeta[e] ~ dnorm(0,0.001)
            }

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
            }

        }
    }
    if(n.x.factor==0 && n.x.nonclim>0.5 && n.x.clim==1){
        envmodel <- function(){
            pi <- 3.141592
            piover2 <- pi/2
            for(i in 1:N) {
                y[i] ~ dbern(p[i])
                for(e in 1:NEnv){
                    x.clim.C[i,e] <- x.clim[i,e]-ax
                    x.clim.beta[i,e,1] <- -beta[e,1]*x.clim.C[i,e]
                    x.clim.beta[i,e,2] <- beta[e,which.beta]*x.clim.C[i,e]
                    x.clim.Calc2[i,e] <- max(x.clim.beta[i,e,1],x.clim.beta[i,e,2])
                    x.clim.Calc[i,e] <- x.clim.Calc2[i,e]*abs(x.clim.C[i,e])
                }
                x.clim.CalcSum[i] <- sum(x.clim.Calc[i,1:NEnv])
                x.clim.Total.tmp[i] <- az - sqrt(x.clim.CalcSum[i])
                x.clim.Total[i] <- min(x.clim.Total.tmp[i], beta0)
                for(e in 1:NnonEnv){
                    x.nonclim.Calc[i,e] <- nonbeta[e]*(x.nonclim[i,e]-mean(x.nonclim[,e]))
                }
                x.nonclim.Total[i] <- sum(x.nonclim.Calc[i,1:NnonEnv])
                logit(p[i]) <- x.clim.Total[i] + x.nonclim.Total[i] + u[u.clique.start[clique[i]]+clique.i[i]-1]
            }
            beta0.prec <- 1/beta0.var
            beta0.diff ~ dnorm(beta0.mu,beta0.prec)
            az ~ dnorm(0,0.1)
            beta0 <- az - exp(beta0.diff)
            for(i in 1:NonSingletonClique){
              clique.length[nonsingleton.clique.list[i]] <- u.clique.end[nonsingleton.clique.list[i]] - u.clique.start[nonsingleton.clique.list[i]] + 1
              adj.clique.length[nonsingleton.clique.list[i]] <- adj.clique.end[nonsingleton.clique.list[i]] - adj.clique.start[nonsingleton.clique.list[i]] + 1
              u[u.clique.start[nonsingleton.clique.list[i]]:u.clique.end[nonsingleton.clique.list[i]]] ~ car.normal(adj[adj.clique.start[nonsingleton.clique.list[i]]:adj.clique.end[nonsingleton.clique.list[i]]], weights[1:adj.clique.length[nonsingleton.clique.list[i]]], num[u.clique.start[nonsingleton.clique.list[i]]:u.clique.end[nonsingleton.clique.list[i]]], tau)
            }
            for(k in 1:adj.clique.end[NClique]) {
                weights[k] <- 1
            }
            tau <- car.tau
            ax.prec <- 1/ax.var
            ax ~ dnorm(ax.mu,ax.prec)%_%I(-1,2)
            for(e in 1:NEnv){
                beta.prec[e,1] <- 1/beta.var[e,1]
                beta.prec[e,2] <- 1/beta.var[e,2]
                beta[e,1] ~ dnorm(beta.mu[e,1],beta.prec[e,1])%_%I(0,)
                beta[e,2] ~ dnorm(beta.mu[e,2],beta.prec[e,2])%_%I(0,)
            }
            for(e in 1:NnonEnv){
                nonbeta[e] ~ dnorm(0,0.001)
            }
        }
    }
    
    # NEW
    
    if(n.x.factor>0.5 && n.x.nonclim==0 && n.x.clim>1.5){
      envmodel <- function(){
        pi <- 3.141592
        piover2 <- pi/2
        
        for(i in 1:N) {
          y[i] ~ dbern(p[i])
          for(e in 1:NEnv){
            x.clim.C[i,e] <- x.clim[i,e]-ax[e]
            x.clim.beta[i,e,1] <- -beta[e,1]*x.clim.C[i,e]
            x.clim.beta[i,e,2] <- beta[e,which.beta[e]]*x.clim.C[i,e]
            x.clim.Calc2[i,e] <- max(x.clim.beta[i,e,1],x.clim.beta[i,e,2])
            x.clim.Calc[i,e] <- x.clim.Calc2[i,e]*abs(x.clim.C[i,e])
          }
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
          for(e in 1:Nfac){
            x.fac.Calc[i,e] <- facbeta[e,x.factor[i,e]]
          }
          x.fac.Total[i] <- sum(x.fac.Calc[i,1:Nfac])
          logit(p[i]) <- x.clim.Total[i] + x.fac.Total[i] + u[u.clique.start[clique[i]]+clique.i[i]-1]
        }
        beta0.prec <- 1/beta0.var
        beta0.diff ~ dnorm(beta0.mu,beta0.prec)
        az ~ dnorm(0,0.1)
        beta0 <- az - exp(beta0.diff)
        for(i in 1:NonSingletonClique){
          clique.length[nonsingleton.clique.list[i]] <- u.clique.end[nonsingleton.clique.list[i]] - u.clique.start[nonsingleton.clique.list[i]] + 1
          adj.clique.length[nonsingleton.clique.list[i]] <- adj.clique.end[nonsingleton.clique.list[i]] - adj.clique.start[nonsingleton.clique.list[i]] + 1
          u[u.clique.start[nonsingleton.clique.list[i]]:u.clique.end[nonsingleton.clique.list[i]]] ~ car.normal(adj[adj.clique.start[nonsingleton.clique.list[i]]:adj.clique.end[nonsingleton.clique.list[i]]], weights[1:adj.clique.length[nonsingleton.clique.list[i]]], num[u.clique.start[nonsingleton.clique.list[i]]:u.clique.end[nonsingleton.clique.list[i]]], tau)
        }
        for(k in 1:adj.clique.end[NClique]) {
          weights[k] <- 1
        }
        tau <- car.tau
        for(e in 1:Nfac){
          facbeta[e,1] <- -sum(facbeta[e,2:x.factor.lengths[e]])
          for(e1 in 2:x.factor.lengths[e]){
            facbeta[e,e1] ~ dnorm(0,0.001)
          }
        }
        
        for(e in 1:NEnv){
          ax.prec[e] <- 1/ax.var[e]
          ax[e] ~ dnorm(ax.mu[e],ax.prec[e])%_%I(-1,2)
          beta.prec[e,1] <- 1/beta.var[e,1]
          beta.prec[e,2] <- 1/beta.var[e,2]
          beta[e,1] ~ dnorm(beta.mu[e,1],beta.prec[e,1])%_%I(0,)
          beta[e,2] ~ dnorm(beta.mu[e,2],beta.prec[e,2])%_%I(0,)
        }
        
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
        }
        
      }
    }
    if(n.x.factor>0.5 && n.x.nonclim==0 && n.x.clim==1){
      envmodel <- function(){
        pi <- 3.141592
        piover2 <- pi/2
        for(i in 1:N) {
          y[i] ~ dbern(p[i])
          for(e in 1:NEnv){
            x.clim.C[i,e] <- x.clim[i,e]-ax
            x.clim.beta[i,e,1] <- -beta[e,1]*x.clim.C[i,e]
            x.clim.beta[i,e,2] <- beta[e,which.beta]*x.clim.C[i,e]
            x.clim.Calc2[i,e] <- max(x.clim.beta[i,e,1],x.clim.beta[i,e,2])
            x.clim.Calc[i,e] <- x.clim.Calc2[i,e]*abs(x.clim.C[i,e])
          }
          x.clim.CalcSum[i] <- sum(x.clim.Calc[i,1:NEnv])
          x.clim.Total.tmp[i] <- az - sqrt(x.clim.CalcSum[i])
          x.clim.Total[i] <- min(x.clim.Total.tmp[i], beta0)
          for(e in 1:Nfac){
            x.fac.Calc[i,e] <- facbeta[e,x.factor[i,e]]
          }
          x.fac.Total[i] <- sum(x.fac.Calc[i,1:Nfac])
          logit(p[i]) <- x.clim.Total[i] + x.fac.Total[i] + u[u.clique.start[clique[i]]+clique.i[i]-1]
        }
        beta0.prec <- 1/beta0.var
        beta0.diff ~ dnorm(beta0.mu,beta0.prec)
        az ~ dnorm(0,0.1)
        beta0 <- az - exp(beta0.diff)
        for(i in 1:NonSingletonClique){
          clique.length[nonsingleton.clique.list[i]] <- u.clique.end[nonsingleton.clique.list[i]] - u.clique.start[nonsingleton.clique.list[i]] + 1
          adj.clique.length[nonsingleton.clique.list[i]] <- adj.clique.end[nonsingleton.clique.list[i]] - adj.clique.start[nonsingleton.clique.list[i]] + 1
          u[u.clique.start[nonsingleton.clique.list[i]]:u.clique.end[nonsingleton.clique.list[i]]] ~ car.normal(adj[adj.clique.start[nonsingleton.clique.list[i]]:adj.clique.end[nonsingleton.clique.list[i]]], weights[1:adj.clique.length[nonsingleton.clique.list[i]]], num[u.clique.start[nonsingleton.clique.list[i]]:u.clique.end[nonsingleton.clique.list[i]]], tau)            }
        for(k in 1:adj.clique.end[NClique]) {
          weights[k] <- 1
        }
        tau <- car.tau
        ax.prec <- 1/ax.var
        ax ~ dnorm(ax.mu,ax.prec)%_%I(-1,2)
        for(e in 1:NEnv){
          beta.prec[e,1] <- 1/beta.var[e,1]
          beta.prec[e,2] <- 1/beta.var[e,2]
          beta[e,1] ~ dnorm(beta.mu[e,1],beta.prec[e,1])%_%I(0,)
          beta[e,2] ~ dnorm(beta.mu[e,2],beta.prec[e,2])%_%I(0,)
        }
        for(e in 1:Nfac){
          facbeta[e,1] <- -sum(facbeta[e,2:x.factor.lengths[e]])
          for(e1 in 2:x.factor.lengths[e]){
            facbeta[e,e1] ~ dnorm(0,0.001)
          }
        }
      }
    }
    if(n.x.factor>0.5 && n.x.nonclim>0.5 && n.x.clim>1.5){
      envmodel <- function(){
        pi <- 3.141592
        piover2 <- pi/2
        
        for(i in 1:N) {
          y[i] ~ dbern(p[i])
          for(e in 1:NEnv){
            x.clim.C[i,e] <- x.clim[i,e]-ax[e]
            x.clim.beta[i,e,1] <- -beta[e,1]*x.clim.C[i,e]
            x.clim.beta[i,e,2] <- beta[e,which.beta[e]]*x.clim.C[i,e]
            x.clim.Calc2[i,e] <- max(x.clim.beta[i,e,1],x.clim.beta[i,e,2])
            x.clim.Calc[i,e] <- x.clim.Calc2[i,e]*abs(x.clim.C[i,e])
          }
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
          for(e in 1:NnonEnv){
            x.nonclim.Calc[i,e] <- nonbeta[e]*(x.nonclim[i,e]-mean(x.nonclim[,e]))
          }
          x.nonclim.Total[i] <- sum(x.nonclim.Calc[i,1:NnonEnv])
          for(e in 1:Nfac){
            x.fac.Calc[i,e] <- facbeta[e,x.factor[i,e]]
          }
          x.fac.Total[i] <- sum(x.fac.Calc[i,1:Nfac])
          logit(p[i]) <- x.clim.Total[i] + x.nonclim.Total[i] + x.fac.Total[i] + u[u.clique.start[clique[i]]+clique.i[i]-1]
        }
        beta0.prec <- 1/beta0.var
        beta0.diff ~ dnorm(beta0.mu,beta0.prec)
        az ~ dnorm(0,0.1)
        beta0 <- az - exp(beta0.diff)
        for(i in 1:NonSingletonClique){
          clique.length[nonsingleton.clique.list[i]] <- u.clique.end[nonsingleton.clique.list[i]] - u.clique.start[nonsingleton.clique.list[i]] + 1
          adj.clique.length[nonsingleton.clique.list[i]] <- adj.clique.end[nonsingleton.clique.list[i]] - adj.clique.start[nonsingleton.clique.list[i]] + 1
          u[u.clique.start[nonsingleton.clique.list[i]]:u.clique.end[nonsingleton.clique.list[i]]] ~ car.normal(adj[adj.clique.start[nonsingleton.clique.list[i]]:adj.clique.end[nonsingleton.clique.list[i]]], weights[1:adj.clique.length[nonsingleton.clique.list[i]]], num[u.clique.start[nonsingleton.clique.list[i]]:u.clique.end[nonsingleton.clique.list[i]]], tau)
        }
        for(k in 1:adj.clique.end[NClique]) {
          weights[k] <- 1
        }
        tau <- car.tau
        
        for(e in 1:NEnv){
          ax.prec[e] <- 1/ax.var[e]
          ax[e] ~ dnorm(ax.mu[e],ax.prec[e])%_%I(-1,2)
          beta.prec[e,1] <- 1/beta.var[e,1]
          beta.prec[e,2] <- 1/beta.var[e,2]
          beta[e,1] ~ dnorm(beta.mu[e,1],beta.prec[e,1])%_%I(0,)
          beta[e,2] ~ dnorm(beta.mu[e,2],beta.prec[e,2])%_%I(0,)
        }
        for(e in 1:NnonEnv){
          nonbeta[e] ~ dnorm(0,0.001)
        }
        for(e in 1:Nfac){
          facbeta[e,1] <- -sum(facbeta[e,2:x.factor.lengths[e]])
          for(e1 in 2:x.factor.lengths[e]){
            facbeta[e,e1] ~ dnorm(0,0.001)
          }
        }
        
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
        }
        
      }
    }
    if(n.x.factor>0.5 && n.x.nonclim>0.5 && n.x.clim==1){
      envmodel <- function(){
        pi <- 3.141592
        piover2 <- pi/2
        for(i in 1:N) {
          y[i] ~ dbern(p[i])
          for(e in 1:NEnv){
            x.clim.C[i,e] <- x.clim[i,e]-ax
            x.clim.beta[i,e,1] <- -beta[e,1]*x.clim.C[i,e]
            x.clim.beta[i,e,2] <- beta[e,which.beta]*x.clim.C[i,e]
            x.clim.Calc2[i,e] <- max(x.clim.beta[i,e,1],x.clim.beta[i,e,2])
            x.clim.Calc[i,e] <- x.clim.Calc2[i,e]*abs(x.clim.C[i,e])
          }
          x.clim.CalcSum[i] <- sum(x.clim.Calc[i,1:NEnv])
          x.clim.Total.tmp[i] <- az - sqrt(x.clim.CalcSum[i])
          x.clim.Total[i] <- min(x.clim.Total.tmp[i], beta0)
          for(e in 1:NnonEnv){
            x.nonclim.Calc[i,e] <- nonbeta[e]*(x.nonclim[i,e]-mean(x.nonclim[,e]))
          }
          x.nonclim.Total[i] <- sum(x.nonclim.Calc[i,1:NnonEnv])
          for(e in 1:Nfac){
            x.fac.Calc[i,e] <- facbeta[e,x.factor[i,e]]
          }
          x.fac.Total[i] <- sum(x.fac.Calc[i,1:Nfac])
          logit(p[i]) <- x.clim.Total[i] + x.nonclim.Total[i] + x.fac.Total[i] + u[u.clique.start[clique[i]]+clique.i[i]-1]
        }
        beta0.prec <- 1/beta0.var
        beta0.diff ~ dnorm(beta0.mu,beta0.prec)
        az ~ dnorm(0,0.1)
        beta0 <- az - exp(beta0.diff)
        for(i in 1:NonSingletonClique){
          clique.length[nonsingleton.clique.list[i]] <- u.clique.end[nonsingleton.clique.list[i]] - u.clique.start[nonsingleton.clique.list[i]] + 1
          adj.clique.length[nonsingleton.clique.list[i]] <- adj.clique.end[nonsingleton.clique.list[i]] - adj.clique.start[nonsingleton.clique.list[i]] + 1
          u[u.clique.start[nonsingleton.clique.list[i]]:u.clique.end[nonsingleton.clique.list[i]]] ~ car.normal(adj[adj.clique.start[nonsingleton.clique.list[i]]:adj.clique.end[nonsingleton.clique.list[i]]], weights[1:adj.clique.length[nonsingleton.clique.list[i]]], num[u.clique.start[nonsingleton.clique.list[i]]:u.clique.end[nonsingleton.clique.list[i]]], tau)
        }
        for(k in 1:adj.clique.end[NClique]) {
          weights[k] <- 1
        }
        tau <- car.tau
        ax.prec <- 1/ax.var
        ax ~ dnorm(ax.mu,ax.prec)%_%I(-1,2)
        for(e in 1:NEnv){
          beta.prec[e,1] <- 1/beta.var[e,1]
          beta.prec[e,2] <- 1/beta.var[e,2]
          beta[e,1] ~ dnorm(beta.mu[e,1],beta.prec[e,1])%_%I(0,)
          beta[e,2] ~ dnorm(beta.mu[e,2],beta.prec[e,2])%_%I(0,)
        }
        for(e in 1:NnonEnv){
          nonbeta[e] ~ dnorm(0,0.001)
        }
        for(e in 1:Nfac){
          facbeta[e,1] <- -sum(facbeta[e,2:x.factor.lengths[e]])
          for(e1 in 2:x.factor.lengths[e]){
            facbeta[e,e1] ~ dnorm(0,0.001)
          }
        }
      }
    }

    write.model(envmodel, "WinBUGSmodel.txt")
    return("WinBUGSmodel.txt")

}
