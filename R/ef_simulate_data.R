
simulate_bl <- function(n, nCov, model)
{
    ## parameters
    ## ----------
    k = .7
    d = nCov
    pi               = LaplacesDemon::rdirichlet(1, c(1,1,1))
    mu.iota          = runif(1,0,k)
    mu.iota.alpha    = runif(1,0,k)
    mu.chi           = runif(1,k,1)
    mu.chi.alpha     = runif(1,k,1)
    ## sigma.iota       = runif(1,0,.1)
    ## sigma.iota.alpha = runif(1,0,.1)
    ## sigma.chi        = 0.075
    ## sigma.chi.alpha  = 0.075
    ## sigma.tau        = runif(1, 0,.1)
    ## sigma.nu         = runif(1, 0,.1)
    beta.tau           = runif(n=nCov+1, -3,3)
    beta.nu            = runif(n=nCov+1, -3,3)
    if (d>0) {
        x        = as.data.frame(MASS::mvrnorm(n, mu=rep(0,d), Sigma=diag(1,d)))
        names(x) = paste0('x',1:d, sep='')
        mu.tau   = as.matrix(cbind(1,x)) %*% beta.tau
        mu.nu    = as.matrix(cbind(1,x)) %*% beta.nu
    }else{
        mu.tau   = rep(1,n) * beta.tau
        mu.nu    = rep(1,n) * beta.nu
    }
    p.tau    = 1/(1+exp(-mu.tau))
    p.nu     = 1/(1+exp(-mu.nu))
    ## vector with true parameters
    true.theta = unlist(list(n=n, pi=pi,
                             beta.tau      = beta.tau,
                             beta.nu       = beta.nu,
                             ## p.tau         = p.tau,
                             ## p.nu          = p.nu,
                             mu.iota       = mu.iota      , 
                             mu.iota.alpha = mu.iota.alpha, 
                             mu.chi        = mu.chi       , 
                             mu.chi.alpha  = mu.chi.alpha  
                             ))
    
    ## data
    ## ----
    ## latent
    N  = sample(500:100, n, replace=T)
    z  = sample(c(1,2,3), n, prob=pi, replace=T)
    tau=nu=iota=chi=iota.alpha=chi.alpha=NA
    for (i in 1:n)
    {
        tau[i]        = rbinom(1, N[i], prob=p.tau[i])/N[i]
        nu[i]         = rbinom(1, N[i], prob=p.nu[i])/N[i]
        iota[i]       = rbinom(1, N[i], mu.iota)/N[i]
        chi[i]        = rbinom(1, N[i], mu.chi)/N[i]
        iota.alpha[i] = rbinom(1, N[i], mu.iota.alpha)/N[i]
        chi.alpha[i]  = rbinom(1, N[i], mu.chi.alpha)/N[i]
    }
    latent     = list(z=z,tau=tau,nu=nu,iota=iota,chi=chi,iota.alpha=iota.alpha,chi.alpha=chi.alpha)
    ## observed
    ## --------
    a = N * ((z==1) * (1 - tau) +
             (z==2) * (1 - tau) * (1 - iota) +
             (z==3) * (1 - tau) * (1 - chi) )
    w = N * ((z==1) * (tau * nu) +
             (z==2) * (tau*nu + iota*(1-tau) + iota.alpha*(tau)*(1-nu)  ) +
             (z==3) * (tau*nu + chi*(1-tau)  + chi.alpha *(tau)*(1-nu)  ) )
    ## rounding w and a
    w = round(w,0)
    a = round(a,0)
    a = a + (N-w-a)
    if (d>0) {
        data = data.frame(cbind(w = w, a = a, N = N, x))
    }else{
        data = data.frame(cbind(w = w, a = a, N = N))
    }
    return(list(parameters=true.theta, latent=latent, data=data))
}

simulate_rn_no_alpha <- function(n, nCov, model)
{
    ## parameters
    ## ----------
    k = .7
    d = nCov
    alpha            = NA
    pi               = LaplacesDemon::rdirichlet(1, c(1,1,1))
    mu.iota          = runif(1,0,k)
    mu.iota.alpha    = runif(1,0,k)
    mu.chi           = runif(1,k,1)
    mu.chi.alpha     = runif(1,k,1)
    sigma.iota       = runif(1,0,.1)
    sigma.iota.alpha = runif(1,0,.1)
    sigma.chi        = 0.075
    sigma.chi.alpha  = 0.075
    sigma.tau        = runif(1, 0,.1)
    sigma.nu         = runif(1, 0,.1)
    mu.tau           = rnorm(n,.5,.1)
    mu.nu            = rnorm(n,.5,.1)
    if (d>0) {
        x = as.data.frame(MASS::mvrnorm(n, mu=rep(0,d), Sigma=diag(1,d)))
        names(x) = paste0('x',1:d, sep='')
        beta.tau         = lm(mu.tau~.,data=x)$coeff
        beta.nu          = lm(mu.nu~.,data=x)$coeff
    }else{
        beta.tau         = mu.tau[1]
        beta.nu          = mu.nu[1]
    }
    ## vector with true parameters
    true.theta=unlist(list(n=n, pi=pi, sigma.tau=sigma.tau, sigma.nu=sigma.nu, beta.tau=beta.tau, beta.nu=beta.nu,
                           mu.iota       = mu.iota      , sigma.iota       = sigma.iota, 
                           mu.iota.alpha = mu.iota.alpha, sigma.iota.alpha = sigma.iota.alpha,
                           mu.chi        = mu.chi       , sigma.chi        = sigma.chi, 
                           mu.chi.alpha  = mu.chi.alpha , sigma.chi.alpha  = sigma.chi.alpha,
                           alpha=alpha))

    ## data
    ## ----
    ## latent
    z          = sample(c(1,2,3), n, prob=pi, replace=T)
    tau        = msm::rtnorm(n, mu.tau,  sigma.tau,   0,1)
    nu         = msm::rtnorm(n, mu.nu,   sigma.nu,    0,1)
    iota       = msm::rtnorm(n, mu.iota, sigma.iota,  0,1)
    chi        = msm::rtnorm(n, mu.chi,  sigma.chi ,  0,1)
    iota.alpha = msm::rtnorm(n, mu.iota.alpha, sigma.iota.alpha,  0,1)
    chi.alpha  = msm::rtnorm(n, mu.chi.alpha,  sigma.chi.alpha ,  0,1)
    latent     = list(z=z,tau=tau,nu=nu,iota=iota,chi=chi,iota.alpha=iota.alpha,chi.alpha=chi.alpha)
    ## observed
    ## --------
    a = (z==1) * (1 - tau) +
        (z==2) * (1 - tau) * (1 - iota) +
        (z==3) * (1 - tau) * (1 - chi)
    w = (z==1) * (tau * nu) +
        (z==2) * (tau*nu + iota*(1-tau) + iota.alpha*(tau)*(1-nu)  ) +
        (z==3) * (tau*nu + chi*(1-tau)  + chi.alpha *(tau)*(1-nu)  )
    if (d>0) {
        data = data.frame(cbind(w = w, a = a, x))
    }else{
        data = data.frame(cbind(w = w, a = a))
    }

    ## computing the likelihood
    ## ------------------------
    ## ln.lambda.g.tau.inv <- c()
    ## ln.lambda.g.nu.inv  <- c()
    ## ln.lambda.neg       <- c()
    ## loglik              <- c()

    ## detJ.a.inv          <- c()
    ## detJ.w.inv          <- c()
    ## g.tau.inv           <- c()
    ## g.nu.inv            <- c()
    ## k.tau               <- c()
    ## k.nu                <- c()
    ## M                   <- 20
    ## for(i in 1:n){
    ##     detJ.a.inv[i] <-  (z[i] == 1) *  1   +
    ##         (z[i] == 2) * ( 1/(1 - iota[i])  ) +
    ##         (z[i] == 3) * ( 1/(1 - chi[i])  )
    ##     detJ.w.inv[i] <-  (z[i] == 1) * ( 1/ (1 - a[i]) ) +
    ##         (z[i] == 2) * ( (1 - iota[i])/( (1 - iota.alpha[i])*(1 - iota[i] - a[i]) ) ) +
    ##         (z[i] == 3) * ( (1 - chi[i])/( (1 - chi.alpha[i])*(1 - chi[i] - a[i]) ) ) 
    ##     g.tau.inv[i] <-   (z[i] == 1) * (1 - a[i] ) +           
    ##         (z[i] == 2) * (1 - a[i]/(1 - iota[i]) ) +
    ##         (z[i] == 3) * (1 - a[i]/(1 - chi[i]) ) 
    ##     g.nu.inv[i] <-
    ##         (z[i] == 1) * ( w[i]*(1/(1 - a[i])) ) +
    ##         (z[i] == 2) * ( w[i]*(1/(1 - iota[i] - a[i]))*((1 - iota[i])/(1-iota.alpha[i])) - iota.alpha[i]/(1 - iota.alpha[i]) - (a[i]*iota[i]) / ((1 - iota[i] - a[i])*(1-iota.alpha[i])) ) +
    ##         (z[i] == 3) * ( w[i]*(1/(1 - chi[i]  - a[i]))*((1 - chi[i]) /(1-chi.alpha[i]) ) - chi.alpha[i] /(1 - chi.alpha[i])  - (a[i]*chi[i])  / ((1 - chi[i]  - a[i])*(1-chi.alpha[i])) )
    ##     k.tau[i] <-  pnorm(1, mu.tau[i], sigma.tau) - pnorm(0, mu.tau[i], sigma.tau) 
    ##     k.nu[i]  <-  pnorm(1, mu.nu[i],  sigma.nu)  - pnorm(0, mu.nu[i],  sigma.nu)  

    ##     ln.lambda.g.tau.inv[i] <- log( dnorm(g.tau.inv[i], mu.tau[i], sigma.tau) )  - log(k.tau[i])
    ##     ln.lambda.g.nu.inv[i]  <- log( dnorm(g.nu.inv[i], mu.nu[i]  , sigma.nu)  )  - log(k.nu[i])
    ##     ln.lambda.neg[i]       <- - ln.lambda.g.tau.inv[i] - ln.lambda.g.nu.inv[i] + M
    ##     loglik[i]              <- ppois(0,ln.lambda.neg[i])

    ## }

    ## return(list(parameters=true.theta, latent=latent, data=data, lik    = data.frame(loglik=loglik, ln.lambda.neg=ln.lambda.neg,k.tau=k.tau,k.nu=k.nu,z=latent$z)))
    return(list(parameters=true.theta, latent=latent, data=data))
}

## {{{ docs }}}
#' Simulate Election data
#'
#' This function simulates data sets from the Election Forensics Finite Mixture models. 
#'
#'
#' @param n integer, the number of sample points representing observations of election units (e.g., precinct, ballot boxes)
#' @param nCov number of covariates affecting both turnout and votes for the winner in each election unit
#' @param model see \code{\link{eforensics}}
#'
#' @return The function returns a list with a data frame with the simulated data and a sublist with the parameters used to generate the data
#'
#' @export
## }}}
ef_simulateData <- function(n=2000,  nCov=0, model)
{
    if(model=='bl')         {return(simulate_bl(n,nCov,model))}
    if(model=='rn_no_alpha'){return(simulate_rn_no_alpha(n,nCov,model))}
}
