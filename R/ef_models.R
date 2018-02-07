#' @export
rn_no_alpha <- function()
{
   model = "
data{
	for(i in 1:n){
	      zeros[i] <- 0
	}
}
model{
	## Constants
	## ---------
	M	  = 25             ## auxiliary variable for the zero trick (makes the - loglikelihood > 0)
	k         = .7
	a.alpha	  = 1
	b.alpha	  = 1
	a.sigma	  = 1
	b.sigma	  = 1
	sigma.chi = 0.075
	sigma.chi.alpha = 0.075
	for (i in 1:3){
	    psi[i]<-1             ## parameters of the dist of the mixing probabilities  
	}   
	for (i in 1:dxa) {        ## d is the number of covars, and if it is zero, we use only the intercept. That is why we add 1 here
	    mu.beta.tau[i]  <- 0                       
            for (j in 1:dxa) {
            	sigma.beta.tau[i,j] <- ifelse(i == j, 10^(2), 0)
	    }
	}
	for (i in 1:dxw) {        ## d is the number of covars, and if it is zero, we use only the intercept. That is why we add 1 here
	    mu.beta.nu[i]   <- 0
            for (j in 1:dxw) {
	    	sigma.beta.nu[i,j]  <- ifelse(i == j, 10^(2), 0)
	    }
	}

	## Hyperpriors
	## -----------
	pi	         ~ ddirch( psi )                        ## mixing probabilities
   	beta.tau	 ~ dmnorm.vcov(mu.beta.tau, sigma.beta.tau)  ## linear coefficients of expectation of turnout
   	beta.nu	         ~ dmnorm.vcov(mu.beta.nu, sigma.beta.nu)    ## linear coefficients of expectation of votes for the winner
	sigma.iota	 ~ dgamma(a.sigma, b.sigma)
	sigma.iota.alpha ~ dgamma(a.sigma, b.sigma)
	sigma.nu	 ~ dgamma(a.sigma, b.sigma)
	sigma.tau	 ~ dgamma(a.sigma, b.sigma)
	alpha	         ~ dgamma(a.alpha, b.alpha)

	mu.iota		~ dunif(0,k)
	mu.iota.alpha	~ dunif(0,k)
	mu.chi		~ dunif(k,1)
	mu.chi.alpha	~ dunif(k,1)
	for(i in 1:n){
	    ## linear transformation of the parameters of tau and nu
	    mu.tau[i]   = inprod(beta.tau, Xa[i,])
	    mu.nu[i]    = inprod(beta.nu , Xw[i,])

	    ## Priors
	    ## ------        
	    z[i]          ~ dcat( pi )

	    iota[i]       ~  dnorm(mu.iota      , 1/pow(sigma.iota,2))       T(0,.9999)            ## truncated normal
	    iota.alpha[i] ~  dnorm(mu.iota.alpha, 1/pow(sigma.iota.alpha,2)) T(0,.9999)            ## truncated normal

	    chi[i]        ~ dnorm(mu.chi      , 1/pow(sigma.chi,2))       T(0,.9999)
	    chi.alpha[i]  ~ dnorm(mu.chi.alpha, 1/pow(sigma.chi.alpha,2)) T(0,.9999)

	    ## Data model
	    ## ----------
            g.tau.inv[i] <-   (z[i] == 1) * (1 - a[i] ) +           
                              (z[i] == 2) * (1 - a[i]/(1 - iota[i]) ) +
                              (z[i] == 3) * (1 - a[i]/(1 - chi[i]) ) 
            g.nu.inv[i] <-    (z[i] == 1) * ( w[i]*(1/(1 - a[i])) ) +
                              (z[i] == 2) * ( w[i]*(1/(1 - iota[i] - a[i]))*((1 - iota[i])/(1-iota.alpha[i])) - iota.alpha[i]/(1 - iota.alpha[i]) - (a[i]*iota[i]) / ((1 - iota[i] - a[i])*(1-iota.alpha[i])) ) +
                              (z[i] == 3) * ( w[i]*(1/(1 - chi[i]  - a[i]))*((1 - chi[i]) /(1-chi.alpha[i]) ) - chi.alpha[i] /(1 - chi.alpha[i])  - (a[i]*chi[i])  / ((1 - chi[i]  - a[i])*(1-chi.alpha[i])) )
           k.tau[i] <-  pnorm(1, mu.tau[i], 1/(sigma.tau^2)) - pnorm(0, mu.tau[i], 1/(sigma.tau^2)) 
           k.nu[i]  <-  pnorm(1, mu.nu[i],  1/(sigma.nu^2))  - pnorm(0, mu.nu[i],  1/(sigma.nu^2))  

           ln.lambda.g.tau.inv[i] <- (0 <= g.tau.inv[i] && g.tau.inv[i] <= 1) * ( logdensity.norm(g.tau.inv[i], mu.tau[i], 1/(sigma.tau^2)) - log(k.tau[i]))
      	   ln.lambda.g.nu.inv[i]  <- (0 <= g.nu.inv[i]  && g.nu.inv[i]  <= 1) * ( logdensity.norm(g.nu.inv[i], mu.nu[i]  , 1/(sigma.nu^2))  - log(k.nu[i]) )
	   lambda[i]              <- - (ln.lambda.g.tau.inv[i] + ln.lambda.g.nu.inv[i]) + M

	    ## Likelihood
	    zeros[i] ~ dpois(lambda[i])

    }
} "
   invisible(model)
}

#' @export
bl <- function()
{
    "
model{
   ## Constants
   ## ---------
   k = .8
   for (i in 1:3){
        sigma[i]<-1
    }
    tau <- 10^(-3)
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
    beta     ~ dmnorm(mu.beta,tau.beta)
    gamma    ~ dmnorm(mu.gamma,tau.gamma)
    pi[1:3] ~ ddirch( sigma )

    zeta.o   ~ dunif(0,k)
    zeta.a   ~ dunif(0,k)
    chi.o    ~ dunif(k,1)
    chi.a    ~ dunif(k,1)
    ## Priors
    ## ------
    for(j in 1:n){
	z[j]      ~ dcat( pi[1:3] )
        ## Calculate alpha and omega
	alpha[j] <- 1/(1 + exp( - (inprod(gamma,Xa[j,]) ) ) )
       	omega[j] <- 1/(1 + exp( - (inprod(beta, Xw[j,]) ) ) )
	## Set success probability for each s and e and draw some binomial random draw
        ss.o[j] ~ dbin(zeta.o,N[j])
        ss.a[j] ~ dbin(zeta.a,N[j])
        ee.o[j] ~ dbin(chi.o,N[j])
        ee.a[j] ~ dbin(chi.a,N[j])
        ## Given alpha, use that as the success probability and draw a and w from binomial distribution
        aa[j]  ~ dbin(alpha[j],N[j])  ## counts for turnout (A) 
        ww[j]  ~ dbin(omega[j],N[j])  ## counts for turnout (W)
        ## Divide everything by N_i to get the probabilities
    	s.a[j] <- ifelse( ss.a[j]/N[j] == 1, .999, ss.a[j]/N[j]) 		## avoiding 1's in order to compute w.check
	s.o[j] <- ss.o[j]/N[j]
	e.a[j] <- ifelse( ee.a[j]/N[j] == 1, .999, ee.a[j]/N[j])		## avoiding 1's in order to compute w.check
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
"
}

#' @export
rn <- function()
{
    "## data for the zero trick to compute the likelihood
data{
	for(i in 1:n){
	      zeros[i] <- 0
	}
}

model{
	## Constants
	## ---------
	M	  = 25             ## auxiliary variable for the zero trick (makes the - loglikelihood > 0)
	k         = .7
	a.alpha	  = 1
	b.alpha	  = 1
	a.sigma	  = 1
	b.sigma	  = 1
	sigma.chi = 0.075
	sigma.chi.alpha = 0.075
	for (i in 1:3){
	    psi[i]<-1             ## parameters of the dist of the mixing probabilities  
	}   
	for (i in 1:dxa) {        ## d is the number of covars, and if it is zero, we use only the intercept. That is why we add 1 here
	    mu.beta.tau[i]  <- 0                       
            for (j in 1:dxa) {
            	sigma.beta.tau[i,j] <- ifelse(i == j, 10^(2), 0)
	    }
	}
	for (i in 1:dxw) {        ## d is the number of covars, and if it is zero, we use only the intercept. That is why we add 1 here
	    mu.beta.nu[i]   <- 0
            for (j in 1:dxw) {
	    	sigma.beta.nu[i,j]  <- ifelse(i == j, 10^(2), 0)
	    }
	}

	## Hyperpriors
	## -----------
	pi	         ~ ddirch( psi )                        ## mixing probabilities
   	beta.tau	 ~ dmnorm.vcov(mu.beta.tau, sigma.beta.tau)  ## linear coefficients of expectation of turnout
   	beta.nu	         ~ dmnorm.vcov(mu.beta.nu, sigma.beta.nu)    ## linear coefficients of expectation of votes for the winner
	sigma.iota	 ~ dgamma(a.sigma, b.sigma)
	sigma.iota.alpha ~ dgamma(a.sigma, b.sigma)
	sigma.nu	 ~ dgamma(a.sigma, b.sigma)
	sigma.tau	 ~ dgamma(a.sigma, b.sigma)
	alpha	         ~ dgamma(a.alpha, b.alpha)

	mu.iota		~ dunif(0,k)
	mu.iota.alpha	~ dunif(0,k)
	mu.chi		~ dunif(k,1)
	mu.chi.alpha	~ dunif(k,1)
	for(i in 1:n){
	    ## linear transformation of the parameters of tau and nu
	    mu.tau[i]   = inprod(beta.tau, Xa[i,])
	    mu.nu[i]    = inprod(beta.nu , Xw[i,])

	    ## Priors
	    ## ------        
	    z[i]          ~ dcat( pi )

	    iota[i]       ~  dnorm(mu.iota      , 1/pow(sigma.iota,2))       T(0,.9999)            ## truncated normal
	    iota.alpha[i] ~  dnorm(mu.iota.alpha, 1/pow(sigma.iota.alpha,2)) T(0,.9999)            ## truncated normal

	    chi[i]        ~ dnorm(mu.chi      , 1/pow(sigma.chi,2))       T(0,.9999)
	    chi.alpha[i]  ~ dnorm(mu.chi.alpha, 1/pow(sigma.chi.alpha,2)) T(0,.9999)

	    ## Data model
	    ## ----------
            g.tau.inv[i] <-   (z[i] == 1) * (1 - a[i] ) +           
                              (z[i] == 2) * (1 - a[i]/(1 - iota[i]) ) +
                              (z[i] == 3) * (1 - a[i]/(1 - chi[i]) ) 
            g.nu.inv[i] <-    (z[i] == 1) * ( w[i]*(1/(1 - a[i])) ) +
                              (z[i] == 2) * ( w[i]*(1/(1 - iota[i] - a[i]))*((1 - iota[i])/(1-iota.alpha[i])) - iota.alpha[i]/(1 - iota.alpha[i]) - (a[i]*iota[i]) / ((1 - iota[i] - a[i])*(1-iota.alpha[i])) ) +
                              (z[i] == 3) * ( w[i]*(1/(1 - chi[i]  - a[i]))*((1 - chi[i]) /(1-chi.alpha[i]) ) - chi.alpha[i] /(1 - chi.alpha[i])  - (a[i]*chi[i])  / ((1 - chi[i]  - a[i])*(1-chi.alpha[i])) )

	   g.tau.inv[i] ~ dnorm(mu.tau[i], 1/pow(sigma.tau,2)) T(0,1)
	   g.nu.inv[i]  ~ dnorm(mu.nu[i] , 1/pow(sigma.nu ,2)) T(0,1)
           # ln.lambda.g.tau.inv[i] <- logdensity.norm(g.tau.inv[i], mu.tau[i], 1/(sigma.tau^2)) ## - log(k.tau[i])
      	   # ln.lambda.g.nu.inv[i]  <- logdensity.norm(g.nu.inv[i], mu.nu[i]  , 1/(sigma.nu^2))  ## - log(k.nu[i])
	   # lambda[i]              <- - (ln.lambda.g.tau.inv[i] + ln.lambda.g.nu.inv[i]) + M

	    ## Likelihood
	    # zeros[i] ~ dpois(lambda[i])
    }
}

"
}
