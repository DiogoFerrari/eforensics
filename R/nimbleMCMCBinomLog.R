rm(list=ls())
library(devtools)
#install_github("DiogoFerrari/eforensics")
library(eforensics)
library(nimble)
library(MCMCpack)
library(data.table)
#Function to simulate data from bl model
simulate_bl <- function(n, nCov, model)
{
  ## parameters
  ## ----------
  k1         = .5
  k2         = .8
  d         = nCov
  
  #pi        = LaplacesDemon::rdirichlet(1, c(1,1,1))
  #Replace this value of pi to change proportions
  pi        = c(.7,.2,.1)
  mu.iota.m = stats::runif(1,0,k1)
  mu.iota.s = stats::runif(1,0,k1)
  mu.chi.m  = stats::runif(1,k2,1)
  mu.chi.s  = stats::runif(1,k2,1)
  beta.tau  = stats::runif(n=nCov+1, -.25,.25)
  beta.nu   = stats::runif(n=nCov+1, -.25,.25)
  
  if (d>0) {
    x        = as.data.frame(MASS::mvrnorm(n, mu=rep(0,d), Sigma=diag(1,d)))
    names(x) = paste0('x',1:d, sep='')
    mu.tau   = as.matrix(cbind(1,x)) %*% beta.tau
    mu.nu    = as.matrix(cbind(1,x)) %*% beta.nu
  }else{
    mu.tau   = rep(stats::runif(1,.3,.7),n)
    mu.nu    = rep(stats::runif(1,.3,.7),n)
  }
  p.tau    = 1/(1+exp(-mu.tau))
  p.nu     = 1/(1+exp(-mu.nu))
  ## vector with true parameters
  true.theta = unlist(list(n=n, pi=pi,
                           beta.tau  = beta.tau,
                           beta.nu   = beta.nu,
                           mu.iota.m = mu.iota.m, 
                           mu.iota.s = mu.iota.s, 
                           mu.chi.m  = mu.chi.m, 
                           mu.chi.s  = mu.chi.s  
  ))
  
  ## data
  ## ----
  ## latent
  N  = base::sample(500:1000, n, replace=T)
  z  = base::sample(c(1,2,3), n, prob=pi, replace=T)
  tau = nu = iota.m = chi.m = iota.s = chi.s = NA
  for (i in 1:n)
  {
    tau[i]    = stats::rbinom(1, N[i], prob=p.tau[i])  /N[i]
    nu[i]     = stats::rbinom(1, N[i], prob=p.nu[i])   /N[i]
    iota.m[i] = stats::rbinom(1, N[i], prob=mu.iota.m) /N[i]
    chi.m[i]  = stats::rbinom(1, N[i], prob=mu.chi.m)  /N[i]
    iota.s[i] = stats::rbinom(1, N[i], prob=mu.iota.s) /N[i]
    chi.s[i]  = stats::rbinom(1, N[i], prob=mu.chi.s)  /N[i]
  }
  latent     = list(z=z,tau=tau,nu=nu,iota.m=iota.m,chi.m=chi.m,iota.s=iota.s,chi.s=chi.s)
  ## observed
  ## --------
  a = N * ((z==1) * (1 - tau) +
             (z==2) * (1 - tau) * (1 - iota.m) +
             (z==3) * (1 - tau) * (1 - chi.m) )
  w = N * ((z==1) * (tau * nu) +
             (z==2) * (tau*nu + iota.m*(1-tau) + iota.s*(tau)*(1-nu)  ) +
             (z==3) * (tau*nu + chi.m* (1-tau) + chi.s *(tau)*(1-nu)  ) )
  ## rounding w and a
  w = round(w,0)
  a = round(a,0)
  ## a = a + (N-w-a)
  if (d>0) {
    data = data.frame(cbind(w = w, a = a, N = N, x))
  }else{
    data = data.frame(cbind(w = w, a = a, N = N))
  }
  return(list(parameters=true.theta, latent=latent, data=data))
}
#Simulate data
sims <- simulate_bl(500,1,"bl")
dats <- sims$data
w <- dats$w
a <- dats$a
n <- length(w)
Xa <- cbind(1,dats$x1,dats$x2,dats$x3,dats$x4)
Xw <- cbind(1,dats$x1,dats$x2,dats$x3,dats$x4)
N <- dats$N
dxw <- ncol(Xw)
dxa <- ncol(Xa)
k <- .8
tau <- 10^(-3)
#Nimble Model Code - Updated to be code from "bl" model in ef_models.R from eforensics package location
#Minor changes to make vector and matrix references complete for C++ compiler
blmod <- nimbleCode(
  {
    ## Constants
    ## ---------
    #k = .8
    for (i in 1:3){
      sigma[i] <- 1
    }
    for (i in 1:dxw) {
      mu.beta[i]  <- 0
      for (j in 1:dxw) {
        tau.beta[i,j]  <- ifelse(i == j, tau, 0)
      }
    } 
    for (i in 1:dxa) {
      mu.gamma[i] <- 0
      for (j in 1:dxa) {
        tau.gamma[i,j] <- ifelse(i == j, tau, 0)
      }
    } 
    ##  Hyperpriors
    ## ------------
    beta[1:dxw]     ~ dmnorm(mu.beta[1:dxw],tau.beta[1:dxw,1:dxw])
    gamma[1:dxa]    ~ dmnorm(mu.gamma[1:dxa],tau.gamma[1:dxa,1:dxa])
    pi[1:3] ~ ddirch( sigma[1:3] )
    
    zeta.o   ~ dunif(0,k)
    zeta.a   ~ dunif(0,k)
    chi.o    ~ dunif(k,1)
    chi.a    ~ dunif(k,1)
    ## Priors
    ## ------
    for(j in 1:n){
      z[j]      ~ dcat( pi[1:3] )
      ## Calculate alpha and omega
      alpha[j] <- 1/(1 + exp( - (inprod(gamma[1:dxa],Xa[j,1:dxa]) ) ) )
      omega[j] <- 1/(1 + exp( - (inprod(beta[1:dxw], Xw[j,1:dxw]) ) ) )
      ## Set success probability for each s and e and draw some binomial random draw
      ss.o[j] ~ dbin(zeta.o,N[j])
      ss.a[j] ~ dbin(zeta.a,N[j])
      ee.o[j] ~ dbin(chi.o,N[j])
      ee.a[j] ~ dbin(chi.a,N[j])
      ## Given alpha, use that as the success probability and draw a and w from binomial distribution
      aa[j]  ~ dbin(alpha[j],N[j])  ## counts for turnout (A) 
      ww[j]  ~ dbin(omega[j],N[j])  ## counts for turnout (W)
      ## Divide everything by N_i to get the probabilities
      s.a[j] <- ss.a[j]/N[j] 		## avoiding 1's in order to compute w.check
      s.o[j] <- ss.o[j]/N[j]
      e.a[j] <- ee.a[j]/N[j]		## avoiding 1's in order to compute w.check
      e.o[j] <- ee.o[j]/N[j]
      a.prop[j]   <- aa[j]/N[j]
      w.prop[j]   <- ww[j]/N[j]
      
      ## Data model
      ## ----------
      ## Note: a.check in the model are integers, here they are proportion, i.e., (A.check/Nj = a.check ). Similarly, w.check = W.check/Nj
      ## H(Nj*a.check | z,s,e,theta )  a.check: success probability for A.check in each cluster
      a.check[j,1] <- a.prop[j]                  ## Zi=O=1  
      a.check[j,2] <- a.prop[j]*(1 - s.a[j])     ## Zi=I=2
      a.check[j,3] <- a.prop[j]*(1 - e.a[j])     ## Zi=E=3
      
      ## F(Nj*w.check | a.check,z,s,e,theta ) 
      w.check[j,1] <- w.prop[j]*(1 - a.check[j,1])									   		       ## Zi=O=1
      w.check[j,2] <- w.prop[j]*( (1-s.o[j])/(1-s.a[j]) )*(1-s.a[j]-a.check[j,2]) + a.check[j,2] * ( (s.a[j] - s.o[j])/(1-s.a[j]) ) + s.o[j] ## Zi=I=2
      w.check[j,3] <- w.prop[j]*( (1-e.o[j])/(1-e.a[j]) )*(1-e.a[j]-a.check[j,3]) + a.check[j,3] * ( (e.a[j] - e.o[j])/(1-e.a[j]) ) + e.o[j] ## Zi=E=3
      
      ## Likelihood
      a[j] ~ dbin(a.check[j,z[j]],N[j])
      w[j] ~ dbin(w.check[j,z[j]],N[j])
    }
  }
)

#Make feasible sets of inits
#Use function to draw from feasible ranges
#Update later to include prior info
make.bl.inits <- function(dxw,dxa,n,k,N){
  beta = rnorm(dxw,0,1)
  gamma = rnorm(dxa,0,1) 
  pi = c(.8,.1,.1)
  zeta.o = runif(1,0,k) 
  zeta.a = runif(1,0,k) 
  chi.o = runif(1,k,1) 
  chi.a = runif(1,k,1) 
  z = sample(c(1,2,3),n,replace=T)
  ss.o <- c()
  ss.a <- c()
  ee.o <- c()
  ee.a <- c()
  #aa <- c()
  #ww <- c()
  #alpha.t <- 1/(1 + exp(-gamma.i))
  #omega.t <- 1/(1 + exp(-beta.i))
  for(i in 1:n){
    ss.o[i] <- rbinom(1, N[i], zeta.o)
    ss.a[i] <- rbinom(1, N[i], zeta.a)
    ee.o[i] <- rbinom(1, N[i], chi.o)
    ee.a[i] <- rbinom(1, N[i], chi.a)
    #aa[i] <- rbinom(1, N[i], alpha.t)
    #ww[i] <- N[i] - aa.i[i]
  }
  return(list(beta = beta, gamma = gamma, pi = pi, zeta.o = zeta.o, zeta.a = zeta.a, chi.o = chi.o, chi.a = chi.a, z = z, ss.o = ss.o, ss.a = ss.a, ee.o = ee.o, ee.a = ee.a))
}
#Establish the initial nimble model
#Ignore any warnings that say model is not fully initialized.
#This just means that there are non-terminal nodes that don't have meaningful intiial values passed.
#These nodes are computed as functions of the set of initial values
blmodr <- nimbleModel(blmod, constant = list(n = n, k = k, dxa = dxa, dxw = dxw,tau = tau), data = list(w = w, a = a, N = N, Xa = Xa, Xw = Xw), inits = make.bl.inits(dxw = dxw, dxa = dxa, n = n, k = k, N = N))
#Configure monitors for the model
spec <- configureMCMC(blmodr, print=TRUE)
spec$resetMonitors()
#CHANGE THIS IF YOU WANT TO REDUCE RUN TIME
spec$addMonitors(c("a.check","w.check","pi","gamma","zeta.o","zeta.a","chi.o","chi.a","ss.o","ss.a","ee.o","ee.a","aa","ww","z"))
#Build the MCMC routine from the spec
blmodr.mcmc <- buildMCMC(spec)
#Compile everything into one "project"
#Note that this can take a while, but the cost is well worth the time savings when running the model
Cblmodr <- compileNimble(blmodr, blmodr.mcmc)
#Run the model
#This is the approach that is similar to the JAGS way
#There's a different approach to use later
mcmcr <- runMCMC(mcmc = Cblmodr$blmodr.mcmc, niter = 7500, nburnin = 7000, nchains = 2)
#Returns model matrix
#Rbind the list to get one model frame
mcmcr <- lapply(mcmcr,as.data.frame)
mcmcr <- rbindlist(mcmcr)
