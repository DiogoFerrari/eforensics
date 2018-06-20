
simulate_bl <- function(n, nCov, model)
{
    ## parameters
    ## ----------
    k1         = .5
    k2         = .8
    d         = nCov

    pi        = LaplacesDemon::rdirichlet(1, c(1,1,1))
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
    sim_data = list(parameters=true.theta, latent=latent, data=data)
    class(sim_data) = 'eforensics_sim_data'
    return(sim_data)
}

simulate_rn_no_alpha <- function(n, nCov, model)
{
    ## parameters
    ## ----------
    k1           = .5
    k2           = .8
    d            = nCov

    pi           = LaplacesDemon::rdirichlet(1, c(1,1,1))

    mu.iota.m    = stats::runif(1,0,k1)
    mu.iota.s    = stats::runif(1,0,k1)
    mu.chi.m     = stats::runif(1,k2,1)
    mu.chi.s     = stats::runif(1,k2,1)
    mu.tau       = stats::rnorm(n,.5,.1)
    mu.nu        = stats::rnorm(n,.5,.1)

    sigma.iota.m = stats::runif(1,0,.1)
    sigma.iota.s = stats::runif(1,0,.1)
    sigma.chi.m  = 0.075
    sigma.chi.s  = 0.075
    sigma.tau    = stats::runif(1, 0,.1)
    sigma.nu     = stats::runif(1, 0,.1)

    if (d>0) {
        x = as.data.frame(MASS::mvrnorm(n, mu=rep(0,d), Sigma=diag(1,d)))
        names(x) = paste0('x',1:d, sep='')
        beta.tau         = stats::lm(mu.tau~.,data=x)$coeff
        beta.nu          = stats::lm(mu.nu~.,data=x)$coeff
    }else{
        beta.tau         = mu.tau[1]
        beta.nu          = mu.nu[1]
    }

    ## vector with true parameters
    true.theta=unlist(list(n=n, pi=pi,
                           sigma.tau    = sigma.tau, 
                           sigma.nu     = sigma.nu, 
                           beta.tau     = beta.tau, 
                           beta.nu      = beta.nu,
                           mu.iota.m    = mu.iota.m      ,
                           sigma.iota.m = sigma.iota.m, 
                           mu.iota.s    = mu.iota.s, 
                           sigma.iota.s = sigma.iota.s,
                           mu.chi.m     = mu.chi.m       , 
                           sigma.chi.m  = sigma.chi.m, 
                           mu.chi.s     = mu.chi.s , 
                           sigma.chi.s  = sigma.chi.s))

    ## data
    ## ----
    ## latent
    z         = base::sample(c(1,2,3), n, prob=pi, replace=T)
    tau       = msm::rtnorm(n, mu.tau,  sigma.tau,   0,1)
    nu        = msm::rtnorm(n, mu.nu,   sigma.nu,    0,1)
    iota.m    = msm::rtnorm(n, mu.iota.m, sigma.iota.m,  0,1)
    chi.m     = msm::rtnorm(n, mu.chi.m,  sigma.chi.m ,  0,1)
    iota.s    = msm::rtnorm(n, mu.iota.s, sigma.iota.s,  0,1)
    chi.s     = msm::rtnorm(n, mu.chi.s,  sigma.chi.s ,  0,1)
    latent    = list(z=z,tau=tau,nu=nu,iota.m=iota.m,chi.m=chi.m,iota.s=iota.s,chi.s=chi.s)
    ## observed
    ## --------
    a = (z==1) * (1 - tau) +
        (z==2) * (1 - tau) * (1 - iota.m) +
        (z==3) * (1 - tau) * (1 - chi.m)
    w = (z==1) * (tau * nu) +
        (z==2) * (tau*nu + iota.m*(1-tau) + iota.s*(tau)*(1-nu)  ) +
        (z==3) * (tau*nu + chi.m *(1-tau) + chi.s *(tau)*(1-nu)  )
    if (d>0) {
        data = data.frame(cbind(w = w, a = a, x))
    }else{
        data = data.frame(cbind(w = w, a = a))
    }

    sim_data = list(parameters=true.theta, latent=latent, data=data)
    class(sim_data) = 'eforensics_sim_data'
    return(sim_data)
}

simulate_rn <- function(n, nCov, model)
{
    ## parameters
    ## ----------
    k1         = .5
    k2         = .8
    d          = nCov

    alpha      = stats::runif(1, 0,2)
    pi         = LaplacesDemon::rdirichlet(1, c(1,1,1))

    mu.iota.m    = stats::runif(1,0,k1)
    mu.chi.m     = stats::runif(1,k2,1)
    mu.tau       = stats::rnorm(n,.5,.1)
    mu.nu        = stats::rnorm(n,.5,.1)

    sigma.iota.m = stats::runif(1,0,.1)
    sigma.chi.m  = 0.075
    sigma.tau    = stats::runif(1, 0,.1)
    sigma.nu     = stats::runif(1, 0,.1)

    if (d>0) {
        x = as.data.frame(MASS::mvrnorm(n, mu=rep(0,d), Sigma=diag(1,d)))
        names(x) = paste0('x',1:d, sep='')
        beta.tau         = stats::lm(mu.tau~.,data=x)$coeff
        beta.nu          = stats::lm(mu.nu~.,data=x)$coeff
    }else{
        beta.tau         = mu.tau[1]
        beta.nu          = mu.nu[1]
    }

    ## vector with true parameters
    true.theta=unlist(list(n=n, pi=pi, sigma.tau=sigma.tau, sigma.nu=sigma.nu, beta.tau=beta.tau, beta.nu=beta.nu,
                           mu.iota.m       = mu.iota.m      , sigma.iota.m       = sigma.iota.m, 
                           ## mu.iota.m.alpha = mu.iota.m.alpha, sigma.iota.m.alpha = sigma.iota.m.alpha,
                           mu.chi.m        = mu.chi.m       , sigma.chi.m        = sigma.chi.m, 
                           ## mu.chi.m.alpha  = mu.chi.m.alpha , sigma.chi.m.alpha  = sigma.chi.m.alpha,
                           alpha=alpha))

    ## data
    ## ----
    ## latent
    z          = base::sample(c(1,2,3), n, prob=pi, replace=T)
    tau        = msm::rtnorm(n, mu.tau,  sigma.tau,   0,1)
    nu         = msm::rtnorm(n, mu.nu,   sigma.nu,    0,1)
    iota.m     = msm::rtnorm(n, mu.iota.m, sigma.iota.m,  0,1)
    chi.m      = msm::rtnorm(n, mu.chi.m,  sigma.chi.m ,  0,1)
    latent     = list(z=z,tau=tau,nu=nu,iota.m=iota.m,chi.m=chi.m)
    ## observed
    ## --------
    a = (z==1) * (1 - tau) +
        (z==2) * (1 - tau) * (1 - iota.m) +
        (z==3) * (1 - tau) * (1 - chi.m)
    w = (z==1) * (tau * nu) +
        (z==2) * (tau*nu + iota.m*(1-tau) + (iota.m^alpha)*(tau)*(1-nu)  ) +
        (z==3) * (tau*nu + chi.m *(1-tau) + (chi.m ^alpha)*(tau)*(1-nu)  )
    if (d>0) {
        data = data.frame(cbind(w = w, a = a, x))
    }else{
        data = data.frame(cbind(w = w, a = a))
    }

    sim_data = list(parameters=true.theta, latent=latent, data=data)
    class(sim_data) = 'eforensics_sim_data'
    return(sim_data)
}

## {{{ docs }}}
#' Simulate Election data
#'
#' This function simulates data sets from the Election Forensics Finite Mixture models. 
#'
#'
#' @param n integer, the number of sample points representing observations of election units (e.g., precinct, ballot boxes)
#' @param nCov number of covariates affecting both turnout and votes for the winner in each election unit
#' @inheritParams eforensics
#'
#' @return The function returns a list with a data frame with the simulated data and a sublist with the parameters used to generate the data
#'
#' @export
## }}}
ef_simulateData <- function(n=2000,  nCov=0, model)
{
    if(model=='bl')         {return(simulate_bl(n,nCov,model))}
    if(model=='rn_no_alpha'){return(simulate_rn_no_alpha(n,nCov,model))}
    if(model=='rn')         {return(simulate_rn(n,nCov,model))}
}
