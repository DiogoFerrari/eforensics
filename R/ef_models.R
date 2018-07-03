
bl <- function()
{
"model{
	## Constants
	## ---------
	k         = .7
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

        mu.iota.m   ~ dunif(0,k)
        mu.iota.s   ~ dunif(0,k)
        mu.chi.m    ~ dunif(k,1)
        mu.chi.s    ~ dunif(k,1)

        for(j in 1:n){
            ## linear transformation of the parameters of tau and nu
            mu.tau[j] <- 1/(1 + exp( - (inprod(beta.tau,Xa[j,]) ) ) )
            mu.nu[j]  <- 1/(1 + exp( - (inprod(beta.nu, Xw[j,]) ) ) )

            ## Priors
            ## ------
  	    Z[j]      ~ dcat( pi )

	    ## Set success probability for each iota (incremental fraue) and chi (extreme fraud) for both cases os manufactured and stolen votes
            N.iota.s[j] ~ dbin(mu.iota.s,N[j])
            N.iota.m[j] ~ dbin(mu.iota.m,N[j])
            N.chi.s[j]  ~ dbin(mu.chi.s, N[j])
            N.chi.m[j]  ~ dbin(mu.chi.m, N[j])
            ## Computing the propostions
       	    iota.m[j] <- ifelse( N.iota.m[j]/N[j] == 1, .999, N.iota.m[j]/N[j]) 		## avoiding 1's in order to compute w.check
	    iota.s[j] <- N.iota.s[j]/N[j]
	    chi.m[j]  <- ifelse( N.chi.m[j]/N[j] == 1, .999, N.chi.m[j]/N[j])		## avoiding 1's in order to compute w.check
	    chi.s[j]  <- N.chi.s[j]/N[j]

            ## ## Given alpha, use that as the success probability and draw a and w from binomial distribution
            ## N.tau[j]  ~ dbin(mu.tau[j], N[j])  ## counts for turnout 
            ## N.nu[j]   ~ dbin(mu.nu[j] ,N[j])   ## counts for votes for the winner 
            ## ## Computing the propostions
	    ## tau[j]  <- N.tau[j]/N[j]
            ## nu[j]   <- N.nu[j] /N[j]
		
            ## Data model
            ## ----------
            p.a[j] = (Z[j] == 1) * (1 - mu.tau[j]) +
                     (Z[j] == 2) * (1 - mu.tau[j]) * (1 - iota.m[j]) +
                     (Z[j] == 3) * (1 - mu.tau[j]) * (1 - chi.m[j])

            p.w[j]  = (Z[j] == 1) * ( mu.nu[j] * (1 - a[j]/N[j]) ) +
                      (Z[j] == 2) * ( mu.nu[j] * ( (1-iota.s[j]) / (1-iota.m[j]) )*( 1-iota.m[j]-a[j]/N[j]) + a[j]/N[j] * ( (iota.m[j] - iota.s[j]) / (1-iota.m[j]) ) + iota.s[j] ) +
                      (Z[j] == 3) * ( mu.nu[j] * ( (1-chi.s[j]) /  (1-chi.m[j])  )*( 1-chi.m[j] -a[j]/N[j]) + a[j]/N[j] * ( (chi.m[j] - chi.s[j]  ) / (1-chi.m[j])  ) + chi.s[j]  )

           ## Likelihood
           a[j] ~ dbin(p.a[j],N[j])
           w[j] ~ dbin(p.w[j],N[j])
       }
}
"
}

bl.vd <- function()
{
  "model{
  ## Constants
  ## ---------
  k         = .7
  for(i in 1:3){
    psi.p[i] <- ifelse(i == 1, 1, .5)
  }
  for (i in 1:3){
    psi.i[i] ~ dbin(psi.p[i],1)
    psi.v[i] ~ dgamma(1,1)
    psi[i] <- psi.i[i]*psi.v[i]
  }

  for(i in 1:3){
    pi[i]	         <-  psi[i]/sum(psi[1:3]) ## mixing probabilities
  }
  pi1 <- pi[1]
  pi2 <- pi[1] + pi[2]

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
  beta.tau	 ~ dmnorm.vcov(mu.beta.tau, sigma.beta.tau)  ## linear coefficients of expectation of turnout
  beta.nu	         ~ dmnorm.vcov(mu.beta.nu, sigma.beta.nu)    ## linear coefficients of expectation of votes for the winner
  
  mu.iota.m   ~ dunif(0,k)
  mu.iota.s   ~ dunif(0,k)
  mu.chi.m    ~ dunif(k,1)
  mu.chi.s    ~ dunif(k,1)
  
  for(j in 1:n){
  ## linear transformation of the parameters of tau and nu
  mu.tau[j] <- 1/(1 + exp( - (inprod(beta.tau,Xa[j,]) ) ) )
  mu.nu[j]  <- 1/(1 + exp( - (inprod(beta.nu, Xw[j,]) ) ) )
  
  ## Priors
  ## ------
  Zn[j] ~ dunif(0,1)
  Zn1[j] <- ifelse(Zn[j] <= pi1, 0, 1)
  Zn2[j] <- ifelse(Zn[j] > pi2, 1 , 0)
  Z[j]  <-  Zn1[j] + Zn2[j] + 1
  
  ## Set success probability for each iota (incremental fraue) and chi (extreme fraud) for both cases os manufactured and stolen votes
  N.iota.s[j] ~ dbin(mu.iota.s,N[j])
  N.iota.m[j] ~ dbin(mu.iota.m,N[j])
  N.chi.s[j]  ~ dbin(mu.chi.s, N[j])
  N.chi.m[j]  ~ dbin(mu.chi.m, N[j])
  ## Computing the propostions
  iota.m[j] <- ifelse( N.iota.m[j]/N[j] == 1, .999, N.iota.m[j]/N[j]) 		## avoiding 1's in order to compute w.check
  iota.s[j] <- N.iota.s[j]/N[j]
  chi.m[j]  <- ifelse( N.chi.m[j]/N[j] == 1, .999, N.chi.m[j]/N[j])		## avoiding 1's in order to compute w.check
  chi.s[j]  <- N.chi.s[j]/N[j]
  
  ## ## Given alpha, use that as the success probability and draw a and w from binomial distribution
  ## N.tau[j]  ~ dbin(mu.tau[j], N[j])  ## counts for turnout 
  ## N.nu[j]   ~ dbin(mu.nu[j] ,N[j])   ## counts for votes for the winner 
  ## ## Computing the propostions
  ## tau[j]  <- N.tau[j]/N[j]
  ## nu[j]   <- N.nu[j] /N[j]
  
  ## Data model
  ## ----------
  p.a[j] = (Z[j] == 1) * (1 - mu.tau[j]) +
  (Z[j] == 2) * (1 - mu.tau[j]) * (1 - iota.m[j]) +
  (Z[j] == 3) * (1 - mu.tau[j]) * (1 - chi.m[j])
  
  p.w[j]  = (Z[j] == 1) * ( mu.nu[j] * (1 - a[j]/N[j]) ) +
  (Z[j] == 2) * ( mu.nu[j] * ( (1-iota.s[j]) / (1-iota.m[j]) )*( 1-iota.m[j]-a[j]/N[j]) + a[j]/N[j] * ( (iota.m[j] - iota.s[j]) / (1-iota.m[j]) ) + iota.s[j] ) +
  (Z[j] == 3) * ( mu.nu[j] * ( (1-chi.s[j]) /  (1-chi.m[j])  )*( 1-chi.m[j] -a[j]/N[j]) + a[j]/N[j] * ( (chi.m[j] - chi.s[j]  ) / (1-chi.m[j])  ) + chi.s[j]  )
  
  ## Likelihood
  a[j] ~ dbin(p.a[j],N[j])
  w[j] ~ dbin(p.w[j],N[j])
  }
}
"
}


bl_working <- function()
{
"model{
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

## =====================================================
## restricted normals
## =====================================================


rn <- function()
{
"data{
	for(i in 1:n){
	    ## zeros[i] <- 0
            ones[i] <- 1
	}
}
model{
	## Constants
	## ---------
	M	  = 1000             ## auxiliary variable for the zero trick (makes the - loglikelihood > 0)
	k         = .7
	a.alpha	  = 1
	b.alpha	  = 1
	a.sigma	  = 1
	b.sigma	  = 1
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
	alpha	         ~ dgamma(a.alpha, b.alpha)

	sigma.iota.m	 ~ dgamma(a.sigma, b.sigma)
	sigma.chi        = 0.075
	sigma.nu	 ~ dgamma(a.sigma, b.sigma)
	sigma.tau	 ~ dgamma(a.sigma, b.sigma)

	## sigma.iota.m	 ~ dunif(0,.1)
	## sigma.chi        = 0.075
	## sigma.nu	 ~ dunif(0,.1)
	## sigma.tau	 ~ dunif(0,.1)

	mu.iota.m		~ dunif(0,k)
	mu.chi.m		~ dunif(k,1)

	for(i in 1:n){
	    ## linear transformation of the parameters of tau and nu
	    mu.tau[i]   = inprod(beta.tau, Xa[i,])
	    mu.nu[i]    = inprod(beta.nu , Xw[i,])

	    ## Priors
	    ## ------        
	    Z[i]          ~  dcat( pi )

	    iota.m[i]     ~  dnorm(mu.iota.m , 1/pow(sigma.iota.m,2)) T(0,.9999)     
	    chi.m[i]      ~  dnorm(mu.chi.m  , 1/pow(sigma.chi,2) ) T(0,.9999)
	    iota.alpha[i] <- pow(iota.m[i], alpha)
	    chi.alpha[i]  <- pow(chi.m[i], alpha)

	    ## Data model
	    ## ----------
            detJ.a.inv[i] <-  (Z[i] == 1) *  1   +
                              (Z[i] == 2) * ( 1/(1 - iota.m[i])  ) +
                              (Z[i] == 3) * ( 1/(1 - chi.m[i])  )
            detJ.w.inv[i] <-  (Z[i] == 1) * ( 1/ (1 - min(a[i],.9999)) ) +
                              (Z[i] == 2) * ( (1 - iota.m[i])/( (1 - iota.alpha[i])*(1 - iota.m[i] - a[i]) ) ) +
                              (Z[i] == 3) * ( (1 - chi.m[i]) /( (1 - chi.alpha[i] )*(1 - chi.m[i]  - a[i]) ) ) 

            g.inv.tau[i] <-   (Z[i] == 1) * (1 -  a[i] ) +           
                              (Z[i] == 2) * (1 - a[i]/(1 - iota.m[i]) )  +
                              (Z[i] == 3) * (1 - a[i]/(1 - chi.m[i] ) ) 
            g.inv.nu[i]  <-   (Z[i] == 1) * ( w[i]*(1/(1 - min(a[i], .9999))) ) +
                              (Z[i] == 2) * ( w[i]*(1/(1 - iota.m[i] - min(a[i], .9999)))*((1 - iota.m[i])/(1-iota.alpha[i])) - iota.alpha[i]/(1 - iota.alpha[i]) - (a[i]*iota.m[i]) / ((1 - iota.m[i] - min(a[i], .9999))*(1-iota.alpha[i])) ) +
                              (Z[i] == 3) * ( w[i]*(1/(1 - chi.m[i]  - min(a[i], .9999)))*((1 - chi.m[i]) /(1-chi.alpha[i]) ) - chi.alpha[i] /(1 - chi.alpha[i])  - (a[i]*chi.m[i])  / ((1 - chi.m[i]  - min(a[i], .9999))*(1-chi.alpha[i])) )

            g.inv.tau.scaled[i] = (g.inv.tau[i] - mu.tau[i]) / sigma.tau
            b.tau[i]            = (  1    - mu.tau[i]) / sigma.tau
            a.tau[i]            = (  0    - mu.tau[i]) / sigma.tau
            g.inv.nu.scaled[i]  = (g.inv.nu[i]  - mu.nu[i])  / sigma.nu
            b.nu[i]             = (  1    - mu.nu[i])  / sigma.nu
            a.nu[i]             = (  0    - mu.nu[i])  / sigma.nu

            k.tau[i] <-  pnorm(b.tau[i], 0, 1) - pnorm(a.tau[i], 0, 1) 
            k.nu[i]  <-  pnorm(b.nu[i] , 0, 1) - pnorm(a.nu[i] , 0, 1)  

            f.a[i]   <- (0 <= g.inv.tau[i] && g.inv.tau[i] <=1) * (  (dnorm(g.inv.tau.scaled[i], 0, 1)*detJ.a.inv[i]) / ( sigma.tau * k.tau[i] )  )
      	    f.w[i]   <- (0 <= g.inv.nu[i]  && g.inv.nu[i]  <=1) * (  (dnorm(g.inv.nu.scaled[i] , 0, 1)*detJ.w.inv[i]) / ( sigma.nu  * k.nu[i]  )  )
            f.w.a[i] <- ( f.a[i] * f.w[i] ) / M

            ones[i] ~ dbern(f.w.a[i])



    }
} 
"}


normal <- function()
{
"data{
	for(i in 1:n){
            ones[i] <- 1
	}
}
model{
	## Constants
	## ---------
	M	  = 10             ## auxiliary variable for the zero trick (makes the - loglikelihood > 0)
	k         = .7
	a.alpha	  = 1
	b.alpha	  = 1
	a.sigma	  = 1
	b.sigma	  = 1
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
	alpha	         ~ dgamma(a.alpha, b.alpha)

	sigma.iota.m	 ~ dgamma(a.sigma, b.sigma)
	sigma.chi        = 0.075
	sigma.nu	 ~ dgamma(a.sigma, b.sigma)
	sigma.tau	 ~ dgamma(a.sigma, b.sigma)

	## sigma.iota.m	 ~ dunif(0,.1)
	## sigma.chi        = 0.075
	## sigma.nu	 ~ dunif(0,.1)
	## sigma.tau	 ~ dunif(0,.1)

	mu.iota.m		~ dunif(0,k)
	mu.chi.m		~ dunif(k,1)

	for(i in 1:n){
	    ## linear transformation of the parameters of tau and nu
	    mu.tau[i]   = inprod(beta.tau, Xa[i,])
	    mu.nu[i]    = inprod(beta.nu , Xw[i,])

	    ## Priors
	    ## ------        
	    Z[i]          ~  dcat( pi )

	    iota.m[i]     ~  dnorm(mu.iota.m , 1/pow(sigma.iota.m,2)) T(0,.9999)     
	    chi.m[i]      ~  dnorm(mu.chi.m  , 1/pow(sigma.chi,2) ) T(0,.9999)
	    iota.alpha[i] <- pow(iota.m[i], alpha)
	    chi.alpha[i]  <- pow(chi.m[i], alpha)

	    ## Data model
	    ## ----------
            detJ.a.inv[i] <-  (Z[i] == 1) *  1   +
                              (Z[i] == 2) * ( 1/(1 - iota.m[i])  ) +
                              (Z[i] == 3) * ( 1/(1 - chi.m[i])  )
            detJ.w.inv[i] <-  (Z[i] == 1) * ( 1/ (1 - min(a[i],.9999)) ) +
                              (Z[i] == 2) * ( (1 - iota.m[i])/( (1 - iota.alpha[i])*(1 - iota.m[i] - a[i]) ) ) +
                              (Z[i] == 3) * ( (1 - chi.m[i]) /( (1 - chi.alpha[i] )*(1 - chi.m[i]  - a[i]) ) ) 

            g.inv.tau[i] <-   (Z[i] == 1) * (1 -  a[i] ) +           
                              (Z[i] == 2) * (1 - a[i]/(1 - iota.m[i]) )  +
                              (Z[i] == 3) * (1 - a[i]/(1 - chi.m[i] ) ) 
            g.inv.nu[i]  <-   (Z[i] == 1) * ( w[i]*(1/(1 - min(a[i], .9999))) ) +
                              (Z[i] == 2) * ( w[i]*(1/(1 - iota.m[i] - min(a[i], .9999)))*((1 - iota.m[i])/(1-iota.alpha[i])) - iota.alpha[i]/(1 - iota.alpha[i]) - (a[i]*iota.m[i]) / ((1 - iota.m[i] - min(a[i], .9999))*(1-iota.alpha[i])) ) +
                              (Z[i] == 3) * ( w[i]*(1/(1 - chi.m[i]  - min(a[i], .9999)))*((1 - chi.m[i]) /(1-chi.alpha[i]) ) - chi.alpha[i] /(1 - chi.alpha[i])  - (a[i]*chi.m[i])  / ((1 - chi.m[i]  - min(a[i], .9999))*(1-chi.alpha[i])) )

            g.inv.tau.scaled[i] = (g.inv.tau[i] - mu.tau[i]) / sigma.tau
            b.tau[i]            = (  1    - mu.tau[i]) / sigma.tau
            a.tau[i]            = (  0    - mu.tau[i]) / sigma.tau
            g.inv.nu.scaled[i]  = (g.inv.nu[i]  - mu.nu[i])  / sigma.nu
            b.nu[i]             = (  1    - mu.nu[i])  / sigma.nu
            a.nu[i]             = (  0    - mu.nu[i])  / sigma.nu

            k.tau[i] <-  pnorm(b.tau[i], 0, 1) - pnorm(a.tau[i], 0, 1) 
            k.nu[i]  <-  pnorm(b.nu[i] , 0, 1) - pnorm(a.nu[i] , 0, 1)  

            f.a[i]   <- (0 <= g.inv.tau[i] && g.inv.tau[i] <=1) * (dnorm(g.inv.tau.scaled[i], 0, 1))
      	    f.w[i]   <- (0 <= g.inv.nu[i]  && g.inv.nu[i]  <=1) * (dnorm(g.inv.nu.scaled[i] , 0, 1))
            f.w.a[i] <- ( f.a[i] * f.w[i] ) / M

            ones[i] ~ dbern(f.w.a[i])



    }
} 
"}


rn_no_scaled <- function()
{
"data{
	for(i in 1:n){
	    ## zeros[i] <- 0
            ones[i] <- 1
	}
}
model{
	## Constants
	## ---------
	M	  = 10             ## auxiliary variable for the zero trick (makes the - loglikelihood > 0)
	k         = .7
	a.alpha	  = 1
	b.alpha	  = 1
	a.sigma	  = 1
	b.sigma	  = 1
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
	alpha	         ~ dgamma(a.alpha, b.alpha)

	sigma.iota.m	 ~ dgamma(a.sigma, b.sigma)
	sigma.chi        = 0.075
	sigma.nu	 ~ dgamma(a.sigma, b.sigma)
	sigma.tau	 ~ dgamma(a.sigma, b.sigma)

	## sigma.iota.m	 ~ dunif(0,.5)
	## sigma.chi        = 0.075
	## sigma.nu	 ~ dunif(0,.5)
	## sigma.tau	 ~ dunif(0,.5)

	mu.iota.m		~ dunif(0,k)
	mu.chi.m		~ dunif(k,1)

	for(i in 1:n){
	    ## linear transformation of the parameters of tau and nu
	    mu.tau[i]   = inprod(beta.tau, Xa[i,])
	    mu.nu[i]    = inprod(beta.nu , Xw[i,])

	    ## Priors
	    ## ------        
	    Z[i]          ~  dcat( pi )

	    iota.m[i]     ~  dnorm(mu.iota.m , 1/pow(sigma.iota.m,2)) T(0,.9999)     
	    chi.m[i]      ~  dnorm(mu.chi.m  , 1/pow(sigma.chi,2) ) T(0,.9999)
	    iota.alpha[i] <- pow(iota.m[i], alpha)
	    chi.alpha[i]  <- pow(chi.m[i], alpha)

	    ## Data model
	    ## ----------
            detJ.a.inv[i] <-  (Z[i] == 1) *  1   +
                              (Z[i] == 2) * ( 1/(1 - iota.m[i])  ) +
                              (Z[i] == 3) * ( 1/(1 - chi.m[i])  )
            detJ.w.inv[i] <-  (Z[i] == 1) * ( 1/ (1 - min(a[i],.9999)) ) +
                              (Z[i] == 2) * ( (1 - iota.m[i])/( (1 - iota.alpha[i])*(1 - iota.m[i] - a[i]) ) ) +
                              (Z[i] == 3) * ( (1 - chi.m[i]) /( (1 - chi.alpha[i] )*(1 - chi.m[i]  - a[i]) ) ) 

            g.inv.tau[i] <-   (Z[i] == 1) * (1 -  a[i] ) +           
                              (Z[i] == 2) * (1 - a[i]/(1 - iota.m[i]) )  +
                              (Z[i] == 3) * (1 - a[i]/(1 - chi.m[i] ) ) 
            g.inv.nu[i]  <-   (Z[i] == 1) * ( w[i]*(1/(1 - min(a[i], .9999))) ) +
                              (Z[i] == 2) * ( w[i]*(1/(1 - iota.m[i] - min(a[i], .9999)))*((1 - iota.m[i])/(1-iota.alpha[i])) - iota.alpha[i]/(1 - iota.alpha[i]) - (a[i]*iota.m[i]) / ((1 - iota.m[i] - min(a[i], .9999))*(1-iota.alpha[i])) ) +
                              (Z[i] == 3) * ( w[i]*(1/(1 - chi.m[i]  - min(a[i], .9999)))*((1 - chi.m[i]) /(1-chi.alpha[i]) ) - chi.alpha[i] /(1 - chi.alpha[i])  - (a[i]*chi.m[i])  / ((1 - chi.m[i]  - min(a[i], .9999))*(1-chi.alpha[i])) )

            g.inv.tau.scaled[i] = (g.inv.tau[i] - mu.tau[i]) / sigma.tau
            b.tau[i]            = (  1    - mu.tau[i]) / sigma.tau
            a.tau[i]            = (  0    - mu.tau[i]) / sigma.tau
            g.inv.nu.scaled[i]  = (g.inv.nu[i]  - mu.nu[i])  / sigma.nu
            b.nu[i]             = (  1    - mu.nu[i])  / sigma.nu
            a.nu[i]             = (  0    - mu.nu[i])  / sigma.nu

            k.tau[i] <-  pnorm(b.tau[i], 0, 1) - pnorm(a.tau[i], 0, 1) 
            k.nu[i]  <-  pnorm(b.nu[i] , 0, 1) - pnorm(a.nu[i] , 0, 1)  

            f.a[i]   <- (0 <= g.inv.tau[i] && g.inv.tau[i] <=1) * (  (dnorm(g.inv.tau.scaled[i], 0, 1)) / ( sigma.tau * k.tau[i] )  )
      	    f.w[i]   <- (0 <= g.inv.nu[i]  && g.inv.nu[i]  <=1) * (  (dnorm(g.inv.nu.scaled[i] , 0, 1)) / ( sigma.nu  * k.nu[i]  )  )
            f.w.a[i] <- ( f.a[i] * f.w[i] ) / M

            ones[i] ~ dbern(f.w.a[i])



    }
} 
"}


rn_no_alpha <- function()
{
   model = "
data{
	for(i in 1:n){
              ones[i]  <- 1
	}
}
model{
	## Constants
	## ---------
	M	  = 10             ## auxiliary variable for the zero trick (makes the - loglikelihood > 0)
	k         = .7
	a.alpha	  = 1
	b.alpha	  = 1
	a.sigma	  = 1
	b.sigma	  = 1
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
	sigma.nu	 ~ dgamma(a.sigma, b.sigma)
	sigma.tau	 ~ dgamma(a.sigma, b.sigma)


	mu.iota.m	~ dunif(0,k)
	mu.iota.s	~ dunif(0,k)
	mu.chi.m	~ dunif(k,1)
	mu.chi.s	~ dunif(k,1)
	sigma.iota.m	~ dgamma(a.sigma, b.sigma)
	sigma.iota.s    ~ dgamma(a.sigma, b.sigma)
	sigma.chi.m     = 0.075
	sigma.chi.s     = 0.075


	for(i in 1:n){
	    ## linear transformation of the parameters of tau and nu
	    mu.tau[i]   = inprod(beta.tau, Xa[i,])
	    mu.nu[i]    = inprod(beta.nu , Xw[i,])

	    ## Priors
	    ## ------        
	    Z[i]       ~  dcat( pi )

	    iota.m[i]  ~  dnorm(mu.iota.m , 1/pow(sigma.iota.m,2)) T(0,.9999)            ## truncated normal
	    iota.s[i]  ~  dnorm(mu.iota.s , 1/pow(sigma.iota.s,2)) T(0,.9999)            ## truncated normal
	    chi.m[i]   ~  dnorm(mu.chi.m  , 1/pow(sigma.chi.m,2) ) T(0,.9999)
	    chi.s[i]   ~  dnorm(mu.chi.s  , 1/pow(sigma.chi.s,2) ) T(0,.9999)

	    ## Data model
	    ## ----------
            tau[i] <-   (Z[i] == 1) * (1 -  min(a[i],.9999) ) +           
                        (Z[i] == 2) * (1 - (min(a[i], .9999)/(1 - iota.m[i])) ) +
                        (Z[i] == 3) * (1 - (min(a[i], .9999)/(1 - chi.m[i] )) ) 
            nu[i]  <-   (Z[i] == 1) * ( w[i]*(1/(1 - min(a[i], .9999))) ) +
                        (Z[i] == 2) * ( w[i]*(1/(1 - iota.m[i] - min(a[i], .9999)))*((1 - iota.m[i])/(1-iota.s[i])) - iota.s[i]/(1 - iota.s[i]) - (a[i]*iota.m[i]) / ((1 - iota.m[i] - min(a[i], .9999))*(1-iota.s[i])) ) +
                        (Z[i] == 3) * ( w[i]*(1/(1 - chi.m[i]  - min(a[i], .9999)))*((1 - chi.m[i]) /(1-chi.s[i]) ) - chi.s[i] /(1 - chi.s[i])  - (a[i]*chi.m[i])  / ((1 - chi.m[i]  - min(a[i], .9999))*(1-chi.s[i])) )

            ## k.tau[i] <-  pnorm(1, mu.tau[i], 1/(sigma.tau^2)) - pnorm(0, mu.tau[i], 1/(sigma.tau^2)) 
            ## k.nu[i]  <-  pnorm(1, mu.nu[i],  1/(sigma.nu^2))  - pnorm(0, mu.nu[i],  1/(sigma.nu^2))  

            ## p.tau[i] <- (0 <= tau[i] && tau[i] <=1) * dnorm(tau[i], mu.tau[i], 1/(sigma.tau^2)) / k.tau[i]
      	    ## p.nu[i]  <- (0 <= nu[i]  && nu[i]  <=1) * dnorm(nu[i] , mu.nu[i] , 1/(sigma.nu^2) ) / k.nu[i]
            ## p[i]     <- ( p.tau[i] * p.nu[i] ) / M

            tau.scaled[i] = (tau[i] - mu.tau[i]) / sigma.tau
            b.tau[i]      = (  1    - mu.tau[i]) / sigma.tau
            a.tau[i]      = (  0    - mu.tau[i]) / sigma.tau
            nu.scaled[i]  = (nu[i]  - mu.nu[i]) / sigma.nu
            b.nu[i]       = (  1    - mu.nu[i]) / sigma.nu
            a.nu[i]       = (  0    - mu.nu[i]) / sigma.nu

            k.tau[i] <-  pnorm(b.tau[i], 0, 1) - pnorm(a.tau[i], 0, 1) 
            k.nu[i]  <-  pnorm(b.nu[i] , 0, 1) - pnorm(a.nu[i] , 0, 1)  

            p.tau[i] <- (0 <= tau[i] && tau[i] <=1) * (  dnorm(tau.scaled[i], 0, 1) / ( sigma.tau * k.tau[i] )  )
      	    p.nu[i]  <- (0 <= nu[i]  && nu[i]  <=1) * (  dnorm(nu.scaled[i] , 0, 1) / ( sigma.nu  * k.nu[i]  )  )
            p[i]     <- ( p.tau[i] * p.nu[i] ) / M

            ones[i] ~ dbern(p[i])


    }

    
} "
   invisible(model)
}

## --------------------
## separted likelihoods
## --------------------



rn_sep <- function()
{
"data{
	for(i in 1:n){
            ones.w[i] <- 1
            ones.a[i] <- 1
	}
}
model{
	## Constants
	## ---------
	M	  = 10             ## auxiliary variable for the zero trick (makes the - loglikelihood > 0)
	k         = .7
	a.alpha	  = 1
	b.alpha	  = 1
	a.sigma	  = 1
	b.sigma	  = 1
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
	alpha	         ~ dgamma(a.alpha, b.alpha)

	sigma.iota.m	 ~ dgamma(a.sigma, b.sigma)
	sigma.chi        = 0.075
	sigma.nu	 ~ dgamma(a.sigma, b.sigma)
	sigma.tau	 ~ dgamma(a.sigma, b.sigma)

	## sigma.iota.m	 ~ dunif(0,.1)
	## sigma.chi        = 0.075
	## sigma.nu	 ~ dunif(0,.1)
	## sigma.tau	 ~ dunif(0,.1)

	mu.iota.m		~ dunif(0,k)
	mu.chi.m		~ dunif(k,1)

	for(i in 1:n){
	    ## linear transformation of the parameters of tau and nu
	    mu.tau[i]   = inprod(beta.tau, Xa[i,])
	    mu.nu[i]    = inprod(beta.nu , Xw[i,])

	    ## Priors
	    ## ------        
	    Z[i]          ~  dcat( pi )

	    iota.m[i]     ~  dnorm(mu.iota.m , 1/pow(sigma.iota.m,2)) T(0,.9999)     
	    chi.m[i]      ~  dnorm(mu.chi.m  , 1/pow(sigma.chi,2) ) T(0,.9999)
	    iota.alpha[i] <- pow(iota.m[i], alpha)
	    chi.alpha[i]  <- pow(chi.m[i], alpha)

	    ## Data model
	    ## ----------
            detJ.a.inv[i] <-  (Z[i] == 1) *  1   +
                              (Z[i] == 2) * ( 1/(1 - iota.m[i])  ) +
                              (Z[i] == 3) * ( 1/(1 - chi.m[i])  )
            detJ.w.inv[i] <-  (Z[i] == 1) * ( 1/ (1 - min(a[i],.9999)) ) +
                              (Z[i] == 2) * ( (1 - iota.m[i])/( (1 - iota.alpha[i])*(1 - iota.m[i] - a[i]) ) ) +
                              (Z[i] == 3) * ( (1 - chi.m[i]) /( (1 - chi.alpha[i] )*(1 - chi.m[i]  - a[i]) ) ) 

            g.inv.tau[i] <-   (Z[i] == 1) * (1 -  a[i] ) +           
                              (Z[i] == 2) * (1 - a[i]/(1 - iota.m[i]) )  +
                              (Z[i] == 3) * (1 - a[i]/(1 - chi.m[i] ) ) 
            g.inv.nu[i]  <-   (Z[i] == 1) * ( w[i]*(1/(1 - min(a[i], .9999))) ) +
                              (Z[i] == 2) * ( w[i]*(1/(1 - iota.m[i] - min(a[i], .9999)))*((1 - iota.m[i])/(1-iota.alpha[i])) - iota.alpha[i]/(1 - iota.alpha[i]) - (a[i]*iota.m[i]) / ((1 - iota.m[i] - min(a[i], .9999))*(1-iota.alpha[i])) ) +
                              (Z[i] == 3) * ( w[i]*(1/(1 - chi.m[i]  - min(a[i], .9999)))*((1 - chi.m[i]) /(1-chi.alpha[i]) ) - chi.alpha[i] /(1 - chi.alpha[i])  - (a[i]*chi.m[i])  / ((1 - chi.m[i]  - min(a[i], .9999))*(1-chi.alpha[i])) )

            g.inv.tau.scaled[i] = (g.inv.tau[i] - mu.tau[i]) / sigma.tau
            b.tau[i]            = (  1    - mu.tau[i]) / sigma.tau
            a.tau[i]            = (  0    - mu.tau[i]) / sigma.tau
            g.inv.nu.scaled[i]  = (g.inv.nu[i]  - mu.nu[i])  / sigma.nu
            b.nu[i]             = (  1    - mu.nu[i])  / sigma.nu
            a.nu[i]             = (  0    - mu.nu[i])  / sigma.nu

            k.tau[i] <-  pnorm(b.tau[i], 0, 1) - pnorm(a.tau[i], 0, 1) 
            k.nu[i]  <-  pnorm(b.nu[i] , 0, 1) - pnorm(a.nu[i] , 0, 1)  

            f.a[i]   <- ((0 <= g.inv.tau[i] && g.inv.tau[i] <=1) * (  (dnorm(g.inv.tau.scaled[i], 0, 1)*detJ.a.inv[i]) / ( sigma.tau * k.tau[i] )  )) / M
      	    f.w[i]   <- ((0 <= g.inv.nu[i]  && g.inv.nu[i]  <=1) * (  (dnorm(g.inv.nu.scaled[i] , 0, 1)*detJ.w.inv[i]) / ( sigma.nu  * k.nu[i]  )  )) / M

            ones.w[i] ~ dbern(f.a[i])
            ones.a[i] ~ dbern(f.w[i])



    }
} 
"}


normal_sep <- function()
{
"data{
	for(i in 1:n){
            ones.w[i] <- 1
            ones.a[i] <- 1
	}
}
model{
	## Constants
	## ---------
	M	  = 10             ## auxiliary variable for the zero trick (makes the - loglikelihood > 0)
	k         = .7
	a.alpha	  = 1
	b.alpha	  = 1
	a.sigma	  = 1
	b.sigma	  = 1
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
	alpha	         ~ dgamma(a.alpha, b.alpha)

	sigma.iota.m	 ~ dgamma(a.sigma, b.sigma)
	sigma.chi        = 0.075
	sigma.nu	 ~ dgamma(a.sigma, b.sigma)
	sigma.tau	 ~ dgamma(a.sigma, b.sigma)

	## sigma.iota.m	 ~ dunif(0,.1)
	## sigma.chi        = 0.075
	## sigma.nu	 ~ dunif(0,.1)
	## sigma.tau	 ~ dunif(0,.1)

	mu.iota.m		~ dunif(0,k)
	mu.chi.m		~ dunif(k,1)

	for(i in 1:n){
	    ## linear transformation of the parameters of tau and nu
	    mu.tau[i]   = inprod(beta.tau, Xa[i,])
	    mu.nu[i]    = inprod(beta.nu , Xw[i,])

	    ## Priors
	    ## ------        
	    Z[i]          ~  dcat( pi )

	    iota.m[i]     ~  dnorm(mu.iota.m , 1/pow(sigma.iota.m,2)) T(0,.9999)     
	    chi.m[i]      ~  dnorm(mu.chi.m  , 1/pow(sigma.chi,2) ) T(0,.9999)
	    iota.alpha[i] <- pow(iota.m[i], alpha)
	    chi.alpha[i]  <- pow(chi.m[i], alpha)

	    ## Data model
	    ## ----------
            detJ.a.inv[i] <-  (Z[i] == 1) *  1   +
                              (Z[i] == 2) * ( 1/(1 - iota.m[i])  ) +
                              (Z[i] == 3) * ( 1/(1 - chi.m[i])  )
            detJ.w.inv[i] <-  (Z[i] == 1) * ( 1/ (1 - min(a[i],.9999)) ) +
                              (Z[i] == 2) * ( (1 - iota.m[i])/( (1 - iota.alpha[i])*(1 - iota.m[i] - a[i]) ) ) +
                              (Z[i] == 3) * ( (1 - chi.m[i]) /( (1 - chi.alpha[i] )*(1 - chi.m[i]  - a[i]) ) ) 

            g.inv.tau[i] <-   (Z[i] == 1) * (1 -  a[i] ) +           
                              (Z[i] == 2) * (1 - a[i]/(1 - iota.m[i]) )  +
                              (Z[i] == 3) * (1 - a[i]/(1 - chi.m[i] ) ) 
            g.inv.nu[i]  <-   (Z[i] == 1) * ( w[i]*(1/(1 - min(a[i], .9999))) ) +
                              (Z[i] == 2) * ( w[i]*(1/(1 - iota.m[i] - min(a[i], .9999)))*((1 - iota.m[i])/(1-iota.alpha[i])) - iota.alpha[i]/(1 - iota.alpha[i]) - (a[i]*iota.m[i]) / ((1 - iota.m[i] - min(a[i], .9999))*(1-iota.alpha[i])) ) +
                              (Z[i] == 3) * ( w[i]*(1/(1 - chi.m[i]  - min(a[i], .9999)))*((1 - chi.m[i]) /(1-chi.alpha[i]) ) - chi.alpha[i] /(1 - chi.alpha[i])  - (a[i]*chi.m[i])  / ((1 - chi.m[i]  - min(a[i], .9999))*(1-chi.alpha[i])) )

            g.inv.tau.scaled[i] = (g.inv.tau[i] - mu.tau[i]) / sigma.tau
            b.tau[i]            = (  1    - mu.tau[i]) / sigma.tau
            a.tau[i]            = (  0    - mu.tau[i]) / sigma.tau
            g.inv.nu.scaled[i]  = (g.inv.nu[i]  - mu.nu[i])  / sigma.nu
            b.nu[i]             = (  1    - mu.nu[i])  / sigma.nu
            a.nu[i]             = (  0    - mu.nu[i])  / sigma.nu

            k.tau[i] <-  pnorm(b.tau[i], 0, 1) - pnorm(a.tau[i], 0, 1) 
            k.nu[i]  <-  pnorm(b.nu[i] , 0, 1) - pnorm(a.nu[i] , 0, 1)  

            f.a[i]   <- ( (0 <= g.inv.tau[i] && g.inv.tau[i] <=1) * (dnorm(g.inv.tau.scaled[i], 0, 1) ) )/M
      	    f.w[i]   <- ( (0 <= g.inv.nu[i]  && g.inv.nu[i]  <=1) * (dnorm(g.inv.nu.scaled[i] , 0, 1) ) )/M

            ones.w[i] ~ dbern(f.a[i])
            ones.a[i] ~ dbern(f.w[i])

    }
} 
"}


rn_no_scaled_sep <- function()
{
"data{
	for(i in 1:n){
            ones.w[i] <- 1
            ones.a[i] <- 1
	}
}
model{
	## Constants
	## ---------
	M	  = 10             ## auxiliary variable for the zero trick (makes the - loglikelihood > 0)
	k         = .7
	a.alpha	  = 1
	b.alpha	  = 1
	a.sigma	  = 1
	b.sigma	  = 1
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
	alpha	         ~ dgamma(a.alpha, b.alpha)

	sigma.iota.m	 ~ dgamma(a.sigma, b.sigma)
	sigma.chi        = 0.075
	sigma.nu	 ~ dgamma(a.sigma, b.sigma)
	sigma.tau	 ~ dgamma(a.sigma, b.sigma)

	## sigma.iota.m	 ~ dunif(0,.5)
	## sigma.chi        = 0.075
	## sigma.nu	 ~ dunif(0,.5)
	## sigma.tau	 ~ dunif(0,.5)

	mu.iota.m		~ dunif(0,k)
	mu.chi.m		~ dunif(k,1)

	for(i in 1:n){
	    ## linear transformation of the parameters of tau and nu
	    mu.tau[i]   = inprod(beta.tau, Xa[i,])
	    mu.nu[i]    = inprod(beta.nu , Xw[i,])

	    ## Priors
	    ## ------        
	    Z[i]          ~  dcat( pi )

	    iota.m[i]     ~  dnorm(mu.iota.m , 1/pow(sigma.iota.m,2)) T(0,.9999)     
	    chi.m[i]      ~  dnorm(mu.chi.m  , 1/pow(sigma.chi,2) ) T(0,.9999)
	    iota.alpha[i] <- pow(iota.m[i], alpha)
	    chi.alpha[i]  <- pow(chi.m[i], alpha)

	    ## Data model
	    ## ----------
            detJ.a.inv[i] <-  (Z[i] == 1) *  1   +
                              (Z[i] == 2) * ( 1/(1 - iota.m[i])  ) +
                              (Z[i] == 3) * ( 1/(1 - chi.m[i])  )
            detJ.w.inv[i] <-  (Z[i] == 1) * ( 1/ (1 - min(a[i],.9999)) ) +
                              (Z[i] == 2) * ( (1 - iota.m[i])/( (1 - iota.alpha[i])*(1 - iota.m[i] - a[i]) ) ) +
                              (Z[i] == 3) * ( (1 - chi.m[i]) /( (1 - chi.alpha[i] )*(1 - chi.m[i]  - a[i]) ) ) 

            g.inv.tau[i] <-   (Z[i] == 1) * (1 -  a[i] ) +           
                              (Z[i] == 2) * (1 - a[i]/(1 - iota.m[i]) )  +
                              (Z[i] == 3) * (1 - a[i]/(1 - chi.m[i] ) ) 
            g.inv.nu[i]  <-   (Z[i] == 1) * ( w[i]*(1/(1 - min(a[i], .9999))) ) +
                              (Z[i] == 2) * ( w[i]*(1/(1 - iota.m[i] - min(a[i], .9999)))*((1 - iota.m[i])/(1-iota.alpha[i])) - iota.alpha[i]/(1 - iota.alpha[i]) - (a[i]*iota.m[i]) / ((1 - iota.m[i] - min(a[i], .9999))*(1-iota.alpha[i])) ) +
                              (Z[i] == 3) * ( w[i]*(1/(1 - chi.m[i]  - min(a[i], .9999)))*((1 - chi.m[i]) /(1-chi.alpha[i]) ) - chi.alpha[i] /(1 - chi.alpha[i])  - (a[i]*chi.m[i])  / ((1 - chi.m[i]  - min(a[i], .9999))*(1-chi.alpha[i])) )

            g.inv.tau.scaled[i] = (g.inv.tau[i] - mu.tau[i]) / sigma.tau
            b.tau[i]            = (  1    - mu.tau[i]) / sigma.tau
            a.tau[i]            = (  0    - mu.tau[i]) / sigma.tau
            g.inv.nu.scaled[i]  = (g.inv.nu[i]  - mu.nu[i])  / sigma.nu
            b.nu[i]             = (  1    - mu.nu[i])  / sigma.nu
            a.nu[i]             = (  0    - mu.nu[i])  / sigma.nu

            k.tau[i] <-  pnorm(b.tau[i], 0, 1) - pnorm(a.tau[i], 0, 1) 
            k.nu[i]  <-  pnorm(b.nu[i] , 0, 1) - pnorm(a.nu[i] , 0, 1)  

            f.a[i]   <- ( (0 <= g.inv.tau[i] && g.inv.tau[i] <=1) * (  (dnorm(g.inv.tau.scaled[i], 0, 1)) / ( sigma.tau * k.tau[i] )  ) ) /M
      	    f.w[i]   <- ( (0 <= g.inv.nu[i]  && g.inv.nu[i]  <=1) * (  (dnorm(g.inv.nu.scaled[i] , 0, 1)) / ( sigma.nu  * k.nu[i]  )  ) ) /M

            ones.w[i] ~ dbern(f.a[i])
            ones.a[i] ~ dbern(f.w[i])



    }
} 
"}


rn_no_alpha_sep <- function()
{
   model = "
data{
	for(i in 1:n){
            ones.w[i] <- 1
            ones.a[i] <- 1
	}
}
model{
	## Constants
	## ---------
	M	  = 10             ## auxiliary variable for the zero trick (makes the - loglikelihood > 0)
	k         = .7
	a.alpha	  = 1
	b.alpha	  = 1
	a.sigma	  = 1
	b.sigma	  = 1
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
	sigma.nu	 ~ dgamma(a.sigma, b.sigma)
	sigma.tau	 ~ dgamma(a.sigma, b.sigma)


	mu.iota.m	~ dunif(0,k)
	mu.iota.s	~ dunif(0,k)
	mu.chi.m	~ dunif(k,1)
	mu.chi.s	~ dunif(k,1)
	sigma.iota.m	~ dgamma(a.sigma, b.sigma)
	sigma.iota.s    ~ dgamma(a.sigma, b.sigma)
	sigma.chi.m     = 0.075
	sigma.chi.s     = 0.075


	for(i in 1:n){
	    ## linear transformation of the parameters of tau and nu
	    mu.tau[i]   = inprod(beta.tau, Xa[i,])
	    mu.nu[i]    = inprod(beta.nu , Xw[i,])

	    ## Priors
	    ## ------        
	    Z[i]       ~  dcat( pi )

	    iota.m[i]  ~  dnorm(mu.iota.m , 1/pow(sigma.iota.m,2)) T(0,.9999)            ## truncated normal
	    iota.s[i]  ~  dnorm(mu.iota.s , 1/pow(sigma.iota.s,2)) T(0,.9999)            ## truncated normal
	    chi.m[i]   ~  dnorm(mu.chi.m  , 1/pow(sigma.chi.m,2) ) T(0,.9999)
	    chi.s[i]   ~  dnorm(mu.chi.s  , 1/pow(sigma.chi.s,2) ) T(0,.9999)

	    ## Data model
	    ## ----------
            tau[i] <-   (Z[i] == 1) * (1 -  min(a[i],.9999) ) +           
                        (Z[i] == 2) * (1 - (min(a[i], .9999)/(1 - iota.m[i])) ) +
                        (Z[i] == 3) * (1 - (min(a[i], .9999)/(1 - chi.m[i] )) ) 
            nu[i]  <-   (Z[i] == 1) * ( w[i]*(1/(1 - min(a[i], .9999))) ) +
                        (Z[i] == 2) * ( w[i]*(1/(1 - iota.m[i] - min(a[i], .9999)))*((1 - iota.m[i])/(1-iota.s[i])) - iota.s[i]/(1 - iota.s[i]) - (a[i]*iota.m[i]) / ((1 - iota.m[i] - min(a[i], .9999))*(1-iota.s[i])) ) +
                        (Z[i] == 3) * ( w[i]*(1/(1 - chi.m[i]  - min(a[i], .9999)))*((1 - chi.m[i]) /(1-chi.s[i]) ) - chi.s[i] /(1 - chi.s[i])  - (a[i]*chi.m[i])  / ((1 - chi.m[i]  - min(a[i], .9999))*(1-chi.s[i])) )

            ## k.tau[i] <-  pnorm(1, mu.tau[i], 1/(sigma.tau^2)) - pnorm(0, mu.tau[i], 1/(sigma.tau^2)) 
            ## k.nu[i]  <-  pnorm(1, mu.nu[i],  1/(sigma.nu^2))  - pnorm(0, mu.nu[i],  1/(sigma.nu^2))  

            ## p.tau[i] <- (0 <= tau[i] && tau[i] <=1) * dnorm(tau[i], mu.tau[i], 1/(sigma.tau^2)) / k.tau[i]
      	    ## p.nu[i]  <- (0 <= nu[i]  && nu[i]  <=1) * dnorm(nu[i] , mu.nu[i] , 1/(sigma.nu^2) ) / k.nu[i]
            ## p[i]     <- ( p.tau[i] * p.nu[i] ) / M

            tau.scaled[i] = (tau[i] - mu.tau[i]) / sigma.tau
            b.tau[i]      = (  1    - mu.tau[i]) / sigma.tau
            a.tau[i]      = (  0    - mu.tau[i]) / sigma.tau
            nu.scaled[i]  = (nu[i]  - mu.nu[i]) / sigma.nu
            b.nu[i]       = (  1    - mu.nu[i]) / sigma.nu
            a.nu[i]       = (  0    - mu.nu[i]) / sigma.nu

            k.tau[i] <-  pnorm(b.tau[i], 0, 1) - pnorm(a.tau[i], 0, 1) 
            k.nu[i]  <-  pnorm(b.nu[i] , 0, 1) - pnorm(a.nu[i] , 0, 1)  

            p.tau[i] <- ( (0 <= tau[i] && tau[i] <=1) * (  dnorm(tau.scaled[i], 0, 1) / ( sigma.tau * k.tau[i] )  ) ) / M
      	    p.nu[i]  <- ( (0 <= nu[i]  && nu[i]  <=1) * (  dnorm(nu.scaled[i] , 0, 1) / ( sigma.nu  * k.nu[i]  )  ) ) / M

            ones.w[i] ~ dbern(p.tau[i])
            ones.a[i] ~ dbern(p.nu[i])


    }

    
} "
   invisible(model)
}

##=====================================================
## varying dimensions models
##=====================================================

## =====================================================
## restricted normals
## =====================================================


rn.vd <- function()
{
  "data{
  for(i in 1:n){
  ## zeros[i] <- 0
  ones[i] <- 1
  }
}
model{
## Constants
## ---------
M	  = 1000             ## auxiliary variable for the zero trick (makes the - loglikelihood > 0)
k         = .7
a.alpha	  = 1
b.alpha	  = 1
a.sigma	  = 1
b.sigma	  = 1
for(i in 1:3){
    psi.p[i] <- ifelse(i == 1, 1, .5)
}
for (i in 1:3){
  psi.i[i] ~ dbin(psi.p[i],1)
  psi.v[i] ~ dgamma(1,1)
  psi[i] <- psi.i[i]*psi.v[i]
}

for(i in 1:3){
  pi[i]	         <-  psi[i]/sum(psi[1:3]) ## mixing probabilities
}
pi1 <- pi[1]
pi2 <- pi[1] + pi[2]   
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
## mixing probabilities
beta.tau	 ~ dmnorm.vcov(mu.beta.tau, sigma.beta.tau)  ## linear coefficients of expectation of turnout
beta.nu	         ~ dmnorm.vcov(mu.beta.nu, sigma.beta.nu)    ## linear coefficients of expectation of votes for the winner
alpha	         ~ dgamma(a.alpha, b.alpha)

sigma.iota.m	 ~ dgamma(a.sigma, b.sigma)
sigma.chi        = 0.075
sigma.nu	 ~ dgamma(a.sigma, b.sigma)
sigma.tau	 ~ dgamma(a.sigma, b.sigma)

## sigma.iota.m	 ~ dunif(0,.1)
## sigma.chi        = 0.075
## sigma.nu	 ~ dunif(0,.1)
## sigma.tau	 ~ dunif(0,.1)

mu.iota.m		~ dunif(0,k)
mu.chi.m		~ dunif(k,1)

for(i in 1:n){
## linear transformation of the parameters of tau and nu
mu.tau[i]   = inprod(beta.tau, Xa[i,])
mu.nu[i]    = inprod(beta.nu , Xw[i,])

## Priors
## ------        
Zn[i] ~ dunif(0,1)
  Zn1[i] <- ifelse(Zn[i] <= pi1, 0, 1)
  Zn2[i] <- ifelse(Zn[i] > pi2, 1 , 0)
  Z[i]  <-  Zn1[i] + Zn2[i] + 1

iota.m[i]     ~  dnorm(mu.iota.m , 1/pow(sigma.iota.m,2)) T(0,.9999)     
chi.m[i]      ~  dnorm(mu.chi.m  , 1/pow(sigma.chi,2) ) T(0,.9999)
iota.alpha[i] <- pow(iota.m[i], alpha)
chi.alpha[i]  <- pow(chi.m[i], alpha)

## Data model
## ----------
detJ.a.inv[i] <-  (Z[i] == 1) *  1   +
(Z[i] == 2) * ( 1/(1 - iota.m[i])  ) +
(Z[i] == 3) * ( 1/(1 - chi.m[i])  )
detJ.w.inv[i] <-  (Z[i] == 1) * ( 1/ (1 - min(a[i],.9999)) ) +
(Z[i] == 2) * ( (1 - iota.m[i])/( (1 - iota.alpha[i])*(1 - iota.m[i] - a[i]) ) ) +
(Z[i] == 3) * ( (1 - chi.m[i]) /( (1 - chi.alpha[i] )*(1 - chi.m[i]  - a[i]) ) ) 

g.inv.tau[i] <-   (Z[i] == 1) * (1 -  a[i] ) +           
(Z[i] == 2) * (1 - a[i]/(1 - iota.m[i]) )  +
(Z[i] == 3) * (1 - a[i]/(1 - chi.m[i] ) ) 
g.inv.nu[i]  <-   (Z[i] == 1) * ( w[i]*(1/(1 - min(a[i], .9999))) ) +
(Z[i] == 2) * ( w[i]*(1/(1 - iota.m[i] - min(a[i], .9999)))*((1 - iota.m[i])/(1-iota.alpha[i])) - iota.alpha[i]/(1 - iota.alpha[i]) - (a[i]*iota.m[i]) / ((1 - iota.m[i] - min(a[i], .9999))*(1-iota.alpha[i])) ) +
(Z[i] == 3) * ( w[i]*(1/(1 - chi.m[i]  - min(a[i], .9999)))*((1 - chi.m[i]) /(1-chi.alpha[i]) ) - chi.alpha[i] /(1 - chi.alpha[i])  - (a[i]*chi.m[i])  / ((1 - chi.m[i]  - min(a[i], .9999))*(1-chi.alpha[i])) )

g.inv.tau.scaled[i] = (g.inv.tau[i] - mu.tau[i]) / sigma.tau
b.tau[i]            = (  1    - mu.tau[i]) / sigma.tau
a.tau[i]            = (  0    - mu.tau[i]) / sigma.tau
g.inv.nu.scaled[i]  = (g.inv.nu[i]  - mu.nu[i])  / sigma.nu
b.nu[i]             = (  1    - mu.nu[i])  / sigma.nu
a.nu[i]             = (  0    - mu.nu[i])  / sigma.nu

k.tau[i] <-  pnorm(b.tau[i], 0, 1) - pnorm(a.tau[i], 0, 1) 
k.nu[i]  <-  pnorm(b.nu[i] , 0, 1) - pnorm(a.nu[i] , 0, 1)  

f.a[i]   <- (0 <= g.inv.tau[i] && g.inv.tau[i] <=1) * (  (dnorm(g.inv.tau.scaled[i], 0, 1)*detJ.a.inv[i]) / ( sigma.tau * k.tau[i] )  )
f.w[i]   <- (0 <= g.inv.nu[i]  && g.inv.nu[i]  <=1) * (  (dnorm(g.inv.nu.scaled[i] , 0, 1)*detJ.w.inv[i]) / ( sigma.nu  * k.nu[i]  )  )
f.w.a[i] <- ( f.a[i] * f.w[i] ) / M

ones[i] ~ dbern(f.w.a[i])



}
} 
"}


normal.vd <- function()
{
  "data{
  for(i in 1:n){
  ones[i] <- 1
  }
}
model{
## Constants
## ---------
M	  = 10             ## auxiliary variable for the zero trick (makes the - loglikelihood > 0)
k         = .7
a.alpha	  = 1
b.alpha	  = 1
a.sigma	  = 1
b.sigma	  = 1
for(i in 1:3){
  psi.p[i] <- ifelse(i == 1, 1, .5)
}
for (i in 1:3){
  psi.i[i] ~ dbin(psi.p[i],1)
  psi.v[i] ~ dgamma(1,1)
  psi[i] <- psi.i[i]*psi.v[i]
}

for(i in 1:3){
  pi[i]	         <-  psi[i]/sum(psi[1:3]) ## mixing probabilities
}
pi1 <- pi[1]
pi2 <- pi[1] + pi[2]   
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
beta.tau	 ~ dmnorm.vcov(mu.beta.tau, sigma.beta.tau)  ## linear coefficients of expectation of turnout
beta.nu	         ~ dmnorm.vcov(mu.beta.nu, sigma.beta.nu)    ## linear coefficients of expectation of votes for the winner
alpha	         ~ dgamma(a.alpha, b.alpha)

sigma.iota.m	 ~ dgamma(a.sigma, b.sigma)
sigma.chi        = 0.075
sigma.nu	 ~ dgamma(a.sigma, b.sigma)
sigma.tau	 ~ dgamma(a.sigma, b.sigma)

## sigma.iota.m	 ~ dunif(0,.1)
## sigma.chi        = 0.075
## sigma.nu	 ~ dunif(0,.1)
## sigma.tau	 ~ dunif(0,.1)

mu.iota.m		~ dunif(0,k)
mu.chi.m		~ dunif(k,1)

for(i in 1:n){
## linear transformation of the parameters of tau and nu
mu.tau[i]   = inprod(beta.tau, Xa[i,])
mu.nu[i]    = inprod(beta.nu , Xw[i,])

## Priors
## ------        
Zn[i] ~ dunif(0,1)
  Zn1[i] <- ifelse(Zn[i] <= pi1, 0, 1)
  Zn2[i] <- ifelse(Zn[i] > pi2, 1 , 0)
  Z[i]  <-  Zn1[i] + Zn2[i] + 1

iota.m[i]     ~  dnorm(mu.iota.m , 1/pow(sigma.iota.m,2)) T(0,.9999)     
chi.m[i]      ~  dnorm(mu.chi.m  , 1/pow(sigma.chi,2) ) T(0,.9999)
iota.alpha[i] <- pow(iota.m[i], alpha)
chi.alpha[i]  <- pow(chi.m[i], alpha)

## Data model
## ----------
detJ.a.inv[i] <-  (Z[i] == 1) *  1   +
(Z[i] == 2) * ( 1/(1 - iota.m[i])  ) +
(Z[i] == 3) * ( 1/(1 - chi.m[i])  )
detJ.w.inv[i] <-  (Z[i] == 1) * ( 1/ (1 - min(a[i],.9999)) ) +
(Z[i] == 2) * ( (1 - iota.m[i])/( (1 - iota.alpha[i])*(1 - iota.m[i] - a[i]) ) ) +
(Z[i] == 3) * ( (1 - chi.m[i]) /( (1 - chi.alpha[i] )*(1 - chi.m[i]  - a[i]) ) ) 

g.inv.tau[i] <-   (Z[i] == 1) * (1 -  a[i] ) +           
(Z[i] == 2) * (1 - a[i]/(1 - iota.m[i]) )  +
(Z[i] == 3) * (1 - a[i]/(1 - chi.m[i] ) ) 
g.inv.nu[i]  <-   (Z[i] == 1) * ( w[i]*(1/(1 - min(a[i], .9999))) ) +
(Z[i] == 2) * ( w[i]*(1/(1 - iota.m[i] - min(a[i], .9999)))*((1 - iota.m[i])/(1-iota.alpha[i])) - iota.alpha[i]/(1 - iota.alpha[i]) - (a[i]*iota.m[i]) / ((1 - iota.m[i] - min(a[i], .9999))*(1-iota.alpha[i])) ) +
(Z[i] == 3) * ( w[i]*(1/(1 - chi.m[i]  - min(a[i], .9999)))*((1 - chi.m[i]) /(1-chi.alpha[i]) ) - chi.alpha[i] /(1 - chi.alpha[i])  - (a[i]*chi.m[i])  / ((1 - chi.m[i]  - min(a[i], .9999))*(1-chi.alpha[i])) )

g.inv.tau.scaled[i] = (g.inv.tau[i] - mu.tau[i]) / sigma.tau
b.tau[i]            = (  1    - mu.tau[i]) / sigma.tau
a.tau[i]            = (  0    - mu.tau[i]) / sigma.tau
g.inv.nu.scaled[i]  = (g.inv.nu[i]  - mu.nu[i])  / sigma.nu
b.nu[i]             = (  1    - mu.nu[i])  / sigma.nu
a.nu[i]             = (  0    - mu.nu[i])  / sigma.nu

k.tau[i] <-  pnorm(b.tau[i], 0, 1) - pnorm(a.tau[i], 0, 1) 
k.nu[i]  <-  pnorm(b.nu[i] , 0, 1) - pnorm(a.nu[i] , 0, 1)  

f.a[i]   <- (0 <= g.inv.tau[i] && g.inv.tau[i] <=1) * (dnorm(g.inv.tau.scaled[i], 0, 1))
f.w[i]   <- (0 <= g.inv.nu[i]  && g.inv.nu[i]  <=1) * (dnorm(g.inv.nu.scaled[i] , 0, 1))
f.w.a[i] <- ( f.a[i] * f.w[i] ) / M

ones[i] ~ dbern(f.w.a[i])



}
} 
"}


rn_no_scaled.vd <- function()
{
  "data{
  for(i in 1:n){
  ## zeros[i] <- 0
  ones[i] <- 1
  }
}
model{
## Constants
## ---------
M	  = 10             ## auxiliary variable for the zero trick (makes the - loglikelihood > 0)
k         = .7
a.alpha	  = 1
b.alpha	  = 1
a.sigma	  = 1
b.sigma	  = 1
for(i in 1:3){
  psi.p[i] <- ifelse(i == 1, 1, .5)
}
for (i in 1:3){
  psi.i[i] ~ dbin(psi.p[i],1)
  psi.v[i] ~ dgamma(1,1)
  psi[i] <- psi.i[i]*psi.v[i]
}

for(i in 1:3){
  pi[i]	         <-  psi[i]/sum(psi[1:3]) ## mixing probabilities
}
pi1 <- pi[1]
pi2 <- pi[1] + pi[2]   
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
beta.tau	 ~ dmnorm.vcov(mu.beta.tau, sigma.beta.tau)  ## linear coefficients of expectation of turnout
beta.nu	         ~ dmnorm.vcov(mu.beta.nu, sigma.beta.nu)    ## linear coefficients of expectation of votes for the winner
alpha	         ~ dgamma(a.alpha, b.alpha)

sigma.iota.m	 ~ dgamma(a.sigma, b.sigma)
sigma.chi        = 0.075
sigma.nu	 ~ dgamma(a.sigma, b.sigma)
sigma.tau	 ~ dgamma(a.sigma, b.sigma)

## sigma.iota.m	 ~ dunif(0,.5)
## sigma.chi        = 0.075
## sigma.nu	 ~ dunif(0,.5)
## sigma.tau	 ~ dunif(0,.5)

mu.iota.m		~ dunif(0,k)
mu.chi.m		~ dunif(k,1)

for(i in 1:n){
## linear transformation of the parameters of tau and nu
mu.tau[i]   = inprod(beta.tau, Xa[i,])
mu.nu[i]    = inprod(beta.nu , Xw[i,])

## Priors
## ------        
Zn[i] ~ dunif(0,1)
  Zn1[i] <- ifelse(Zn[i] <= pi1, 0, 1)
  Zn2[i] <- ifelse(Zn[i] > pi2, 1 , 0)
  Z[i]  <-  Zn1[i] + Zn2[i] + 1

iota.m[i]     ~  dnorm(mu.iota.m , 1/pow(sigma.iota.m,2)) T(0,.9999)     
chi.m[i]      ~  dnorm(mu.chi.m  , 1/pow(sigma.chi,2) ) T(0,.9999)
iota.alpha[i] <- pow(iota.m[i], alpha)
chi.alpha[i]  <- pow(chi.m[i], alpha)

## Data model
## ----------
detJ.a.inv[i] <-  (Z[i] == 1) *  1   +
(Z[i] == 2) * ( 1/(1 - iota.m[i])  ) +
(Z[i] == 3) * ( 1/(1 - chi.m[i])  )
detJ.w.inv[i] <-  (Z[i] == 1) * ( 1/ (1 - min(a[i],.9999)) ) +
(Z[i] == 2) * ( (1 - iota.m[i])/( (1 - iota.alpha[i])*(1 - iota.m[i] - a[i]) ) ) +
(Z[i] == 3) * ( (1 - chi.m[i]) /( (1 - chi.alpha[i] )*(1 - chi.m[i]  - a[i]) ) ) 

g.inv.tau[i] <-   (Z[i] == 1) * (1 -  a[i] ) +           
(Z[i] == 2) * (1 - a[i]/(1 - iota.m[i]) )  +
(Z[i] == 3) * (1 - a[i]/(1 - chi.m[i] ) ) 
g.inv.nu[i]  <-   (Z[i] == 1) * ( w[i]*(1/(1 - min(a[i], .9999))) ) +
(Z[i] == 2) * ( w[i]*(1/(1 - iota.m[i] - min(a[i], .9999)))*((1 - iota.m[i])/(1-iota.alpha[i])) - iota.alpha[i]/(1 - iota.alpha[i]) - (a[i]*iota.m[i]) / ((1 - iota.m[i] - min(a[i], .9999))*(1-iota.alpha[i])) ) +
(Z[i] == 3) * ( w[i]*(1/(1 - chi.m[i]  - min(a[i], .9999)))*((1 - chi.m[i]) /(1-chi.alpha[i]) ) - chi.alpha[i] /(1 - chi.alpha[i])  - (a[i]*chi.m[i])  / ((1 - chi.m[i]  - min(a[i], .9999))*(1-chi.alpha[i])) )

g.inv.tau.scaled[i] = (g.inv.tau[i] - mu.tau[i]) / sigma.tau
b.tau[i]            = (  1    - mu.tau[i]) / sigma.tau
a.tau[i]            = (  0    - mu.tau[i]) / sigma.tau
g.inv.nu.scaled[i]  = (g.inv.nu[i]  - mu.nu[i])  / sigma.nu
b.nu[i]             = (  1    - mu.nu[i])  / sigma.nu
a.nu[i]             = (  0    - mu.nu[i])  / sigma.nu

k.tau[i] <-  pnorm(b.tau[i], 0, 1) - pnorm(a.tau[i], 0, 1) 
k.nu[i]  <-  pnorm(b.nu[i] , 0, 1) - pnorm(a.nu[i] , 0, 1)  

f.a[i]   <- (0 <= g.inv.tau[i] && g.inv.tau[i] <=1) * (  (dnorm(g.inv.tau.scaled[i], 0, 1)) / ( sigma.tau * k.tau[i] )  )
f.w[i]   <- (0 <= g.inv.nu[i]  && g.inv.nu[i]  <=1) * (  (dnorm(g.inv.nu.scaled[i] , 0, 1)) / ( sigma.nu  * k.nu[i]  )  )
f.w.a[i] <- ( f.a[i] * f.w[i] ) / M

ones[i] ~ dbern(f.w.a[i])



}
} 
"}


rn_no_alpha.vd <- function()
{
  model = "
  data{
  for(i in 1:n){
  ones[i]  <- 1
  }
  }
  model{
  ## Constants
  ## ---------
  M	  = 10             ## auxiliary variable for the zero trick (makes the - loglikelihood > 0)
  k         = .7
  a.alpha	  = 1
  b.alpha	  = 1
  a.sigma	  = 1
  b.sigma	  = 1
  for(i in 1:3){
    psi.p[i] <- ifelse(i == 1, 1, .5)
  }
  for (i in 1:3){
    psi.i[i] ~ dbin(psi.p[i],1)
    psi.v[i] ~ dgamma(1,1)
    psi[i] <- psi.i[i]*psi.v[i]
  }

  for(i in 1:3){
    pi[i]	         <-  psi[i]/sum(psi[1:3]) ## mixing probabilities
  }
  pi1 <- pi[1]
  pi2 <- pi[1] + pi[2]   
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
  
  beta.tau	 ~ dmnorm.vcov(mu.beta.tau, sigma.beta.tau)  ## linear coefficients of expectation of turnout
  beta.nu	         ~ dmnorm.vcov(mu.beta.nu, sigma.beta.nu)    ## linear coefficients of expectation of votes for the winner
  sigma.nu	 ~ dgamma(a.sigma, b.sigma)
  sigma.tau	 ~ dgamma(a.sigma, b.sigma)
  
  
  mu.iota.m	~ dunif(0,k)
  mu.iota.s	~ dunif(0,k)
  mu.chi.m	~ dunif(k,1)
  mu.chi.s	~ dunif(k,1)
  sigma.iota.m	~ dgamma(a.sigma, b.sigma)
  sigma.iota.s    ~ dgamma(a.sigma, b.sigma)
  sigma.chi.m     = 0.075
  sigma.chi.s     = 0.075
  
  
  for(i in 1:n){
  ## linear transformation of the parameters of tau and nu
  mu.tau[i]   = inprod(beta.tau, Xa[i,])
  mu.nu[i]    = inprod(beta.nu , Xw[i,])
  
  ## Priors
  ## ------        
  Zn[i] ~ dunif(0,1)
  Zn1[i] <- ifelse(Zn[i] <= pi1, 0, 1)
  Zn2[i] <- ifelse(Zn[i] > pi2, 1 , 0)
  Z[i]  <-  Zn1[i] + Zn2[i] + 1
  
  iota.m[i]  ~  dnorm(mu.iota.m , 1/pow(sigma.iota.m,2)) T(0,.9999)            ## truncated normal
  iota.s[i]  ~  dnorm(mu.iota.s , 1/pow(sigma.iota.s,2)) T(0,.9999)            ## truncated normal
  chi.m[i]   ~  dnorm(mu.chi.m  , 1/pow(sigma.chi.m,2) ) T(0,.9999)
  chi.s[i]   ~  dnorm(mu.chi.s  , 1/pow(sigma.chi.s,2) ) T(0,.9999)
  
  ## Data model
  ## ----------
  tau[i] <-   (Z[i] == 1) * (1 -  min(a[i],.9999) ) +           
  (Z[i] == 2) * (1 - (min(a[i], .9999)/(1 - iota.m[i])) ) +
  (Z[i] == 3) * (1 - (min(a[i], .9999)/(1 - chi.m[i] )) ) 
  nu[i]  <-   (Z[i] == 1) * ( w[i]*(1/(1 - min(a[i], .9999))) ) +
  (Z[i] == 2) * ( w[i]*(1/(1 - iota.m[i] - min(a[i], .9999)))*((1 - iota.m[i])/(1-iota.s[i])) - iota.s[i]/(1 - iota.s[i]) - (a[i]*iota.m[i]) / ((1 - iota.m[i] - min(a[i], .9999))*(1-iota.s[i])) ) +
  (Z[i] == 3) * ( w[i]*(1/(1 - chi.m[i]  - min(a[i], .9999)))*((1 - chi.m[i]) /(1-chi.s[i]) ) - chi.s[i] /(1 - chi.s[i])  - (a[i]*chi.m[i])  / ((1 - chi.m[i]  - min(a[i], .9999))*(1-chi.s[i])) )
  
  ## k.tau[i] <-  pnorm(1, mu.tau[i], 1/(sigma.tau^2)) - pnorm(0, mu.tau[i], 1/(sigma.tau^2)) 
  ## k.nu[i]  <-  pnorm(1, mu.nu[i],  1/(sigma.nu^2))  - pnorm(0, mu.nu[i],  1/(sigma.nu^2))  
  
  ## p.tau[i] <- (0 <= tau[i] && tau[i] <=1) * dnorm(tau[i], mu.tau[i], 1/(sigma.tau^2)) / k.tau[i]
  ## p.nu[i]  <- (0 <= nu[i]  && nu[i]  <=1) * dnorm(nu[i] , mu.nu[i] , 1/(sigma.nu^2) ) / k.nu[i]
  ## p[i]     <- ( p.tau[i] * p.nu[i] ) / M
  
  tau.scaled[i] = (tau[i] - mu.tau[i]) / sigma.tau
  b.tau[i]      = (  1    - mu.tau[i]) / sigma.tau
  a.tau[i]      = (  0    - mu.tau[i]) / sigma.tau
  nu.scaled[i]  = (nu[i]  - mu.nu[i]) / sigma.nu
  b.nu[i]       = (  1    - mu.nu[i]) / sigma.nu
  a.nu[i]       = (  0    - mu.nu[i]) / sigma.nu
  
  k.tau[i] <-  pnorm(b.tau[i], 0, 1) - pnorm(a.tau[i], 0, 1) 
  k.nu[i]  <-  pnorm(b.nu[i] , 0, 1) - pnorm(a.nu[i] , 0, 1)  
  
  p.tau[i] <- (0 <= tau[i] && tau[i] <=1) * (  dnorm(tau.scaled[i], 0, 1) / ( sigma.tau * k.tau[i] )  )
  p.nu[i]  <- (0 <= nu[i]  && nu[i]  <=1) * (  dnorm(nu.scaled[i] , 0, 1) / ( sigma.nu  * k.nu[i]  )  )
  p[i]     <- ( p.tau[i] * p.nu[i] ) / M
  
  ones[i] ~ dbern(p[i])
  
  
  }
  
  
  } "
   invisible(model)
}

## --------------------
## separted likelihoods
## --------------------



rn_sep.vd <- function()
{
  "data{
  for(i in 1:n){
  ones.w[i] <- 1
  ones.a[i] <- 1
  }
}
model{
## Constants
## ---------
M	  = 10             ## auxiliary variable for the zero trick (makes the - loglikelihood > 0)
k         = .7
a.alpha	  = 1
b.alpha	  = 1
a.sigma	  = 1
b.sigma	  = 1
for(i in 1:3){
    psi.p[i] <- ifelse(i == 1, 1, .5)
  }
  for (i in 1:3){
    psi.i[i] ~ dbin(psi.p[i],1)
    psi.v[i] ~ dgamma(1,1)
    psi[i] <- psi.i[i]*psi.v[i]
  }

  for(i in 1:3){
    pi[i]	         <-  psi[i]/sum(psi[1:3]) ## mixing probabilities
  }
  pi1 <- pi[1]
  pi2 <- pi[1] + pi[2]
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
beta.tau	 ~ dmnorm.vcov(mu.beta.tau, sigma.beta.tau)  ## linear coefficients of expectation of turnout
beta.nu	         ~ dmnorm.vcov(mu.beta.nu, sigma.beta.nu)    ## linear coefficients of expectation of votes for the winner
alpha	         ~ dgamma(a.alpha, b.alpha)

sigma.iota.m	 ~ dgamma(a.sigma, b.sigma)
sigma.chi        = 0.075
sigma.nu	 ~ dgamma(a.sigma, b.sigma)
sigma.tau	 ~ dgamma(a.sigma, b.sigma)

## sigma.iota.m	 ~ dunif(0,.1)
## sigma.chi        = 0.075
## sigma.nu	 ~ dunif(0,.1)
## sigma.tau	 ~ dunif(0,.1)

mu.iota.m		~ dunif(0,k)
mu.chi.m		~ dunif(k,1)

for(i in 1:n){
## linear transformation of the parameters of tau and nu
mu.tau[i]   = inprod(beta.tau, Xa[i,])
mu.nu[i]    = inprod(beta.nu , Xw[i,])

## Priors
## ------        
Zn[i] ~ dunif(0,1)
  Zn1[i] <- ifelse(Zn[i] <= pi1, 0, 1)
  Zn2[i] <- ifelse(Zn[i] > pi2, 1 , 0)
  Z[i]  <-  Zn1[i] + Zn2[i] + 1

iota.m[i]     ~  dnorm(mu.iota.m , 1/pow(sigma.iota.m,2)) T(0,.9999)     
chi.m[i]      ~  dnorm(mu.chi.m  , 1/pow(sigma.chi,2) ) T(0,.9999)
iota.alpha[i] <- pow(iota.m[i], alpha)
chi.alpha[i]  <- pow(chi.m[i], alpha)

## Data model
## ----------
detJ.a.inv[i] <-  (Z[i] == 1) *  1   +
(Z[i] == 2) * ( 1/(1 - iota.m[i])  ) +
(Z[i] == 3) * ( 1/(1 - chi.m[i])  )
detJ.w.inv[i] <-  (Z[i] == 1) * ( 1/ (1 - min(a[i],.9999)) ) +
(Z[i] == 2) * ( (1 - iota.m[i])/( (1 - iota.alpha[i])*(1 - iota.m[i] - a[i]) ) ) +
(Z[i] == 3) * ( (1 - chi.m[i]) /( (1 - chi.alpha[i] )*(1 - chi.m[i]  - a[i]) ) ) 

g.inv.tau[i] <-   (Z[i] == 1) * (1 -  a[i] ) +           
(Z[i] == 2) * (1 - a[i]/(1 - iota.m[i]) )  +
(Z[i] == 3) * (1 - a[i]/(1 - chi.m[i] ) ) 
g.inv.nu[i]  <-   (Z[i] == 1) * ( w[i]*(1/(1 - min(a[i], .9999))) ) +
(Z[i] == 2) * ( w[i]*(1/(1 - iota.m[i] - min(a[i], .9999)))*((1 - iota.m[i])/(1-iota.alpha[i])) - iota.alpha[i]/(1 - iota.alpha[i]) - (a[i]*iota.m[i]) / ((1 - iota.m[i] - min(a[i], .9999))*(1-iota.alpha[i])) ) +
(Z[i] == 3) * ( w[i]*(1/(1 - chi.m[i]  - min(a[i], .9999)))*((1 - chi.m[i]) /(1-chi.alpha[i]) ) - chi.alpha[i] /(1 - chi.alpha[i])  - (a[i]*chi.m[i])  / ((1 - chi.m[i]  - min(a[i], .9999))*(1-chi.alpha[i])) )

g.inv.tau.scaled[i] = (g.inv.tau[i] - mu.tau[i]) / sigma.tau
b.tau[i]            = (  1    - mu.tau[i]) / sigma.tau
a.tau[i]            = (  0    - mu.tau[i]) / sigma.tau
g.inv.nu.scaled[i]  = (g.inv.nu[i]  - mu.nu[i])  / sigma.nu
b.nu[i]             = (  1    - mu.nu[i])  / sigma.nu
a.nu[i]             = (  0    - mu.nu[i])  / sigma.nu

k.tau[i] <-  pnorm(b.tau[i], 0, 1) - pnorm(a.tau[i], 0, 1) 
k.nu[i]  <-  pnorm(b.nu[i] , 0, 1) - pnorm(a.nu[i] , 0, 1)  

f.a[i]   <- ((0 <= g.inv.tau[i] && g.inv.tau[i] <=1) * (  (dnorm(g.inv.tau.scaled[i], 0, 1)*detJ.a.inv[i]) / ( sigma.tau * k.tau[i] )  )) / M
f.w[i]   <- ((0 <= g.inv.nu[i]  && g.inv.nu[i]  <=1) * (  (dnorm(g.inv.nu.scaled[i] , 0, 1)*detJ.w.inv[i]) / ( sigma.nu  * k.nu[i]  )  )) / M

ones.w[i] ~ dbern(f.a[i])
ones.a[i] ~ dbern(f.w[i])



}
} 
"}


normal_sep.vd <- function()
{
  "data{
  for(i in 1:n){
  ones.w[i] <- 1
  ones.a[i] <- 1
  }
}
model{
## Constants
## ---------
M	  = 10             ## auxiliary variable for the zero trick (makes the - loglikelihood > 0)
k         = .7
a.alpha	  = 1
b.alpha	  = 1
a.sigma	  = 1
b.sigma	  = 1
for(i in 1:3){
    psi.p[i] <- ifelse(i == 1, 1, .5)
  }
  for (i in 1:3){
    psi.i[i] ~ dbin(psi.p[i],1)
    psi.v[i] ~ dgamma(1,1)
    psi[i] <- psi.i[i]*psi.v[i]
  }

  for(i in 1:3){
    pi[i]	         <-  psi[i]/sum(psi[1:3]) ## mixing probabilities
  }
  pi1 <- pi[1]
  pi2 <- pi[1] + pi[2] 
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
beta.tau	 ~ dmnorm.vcov(mu.beta.tau, sigma.beta.tau)  ## linear coefficients of expectation of turnout
beta.nu	         ~ dmnorm.vcov(mu.beta.nu, sigma.beta.nu)    ## linear coefficients of expectation of votes for the winner
alpha	         ~ dgamma(a.alpha, b.alpha)

sigma.iota.m	 ~ dgamma(a.sigma, b.sigma)
sigma.chi        = 0.075
sigma.nu	 ~ dgamma(a.sigma, b.sigma)
sigma.tau	 ~ dgamma(a.sigma, b.sigma)

## sigma.iota.m	 ~ dunif(0,.1)
## sigma.chi        = 0.075
## sigma.nu	 ~ dunif(0,.1)
## sigma.tau	 ~ dunif(0,.1)

mu.iota.m		~ dunif(0,k)
mu.chi.m		~ dunif(k,1)

for(i in 1:n){
## linear transformation of the parameters of tau and nu
mu.tau[i]   = inprod(beta.tau, Xa[i,])
mu.nu[i]    = inprod(beta.nu , Xw[i,])

## Priors
## ------        
Zn[i] ~ dunif(0,1)
  Zn1[i] <- ifelse(Zn[i] <= pi1, 0, 1)
  Zn2[i] <- ifelse(Zn[i] > pi2, 1 , 0)
  Z[i]  <-  Zn1[i] + Zn2[i] + 1

iota.m[i]     ~  dnorm(mu.iota.m , 1/pow(sigma.iota.m,2)) T(0,.9999)     
chi.m[i]      ~  dnorm(mu.chi.m  , 1/pow(sigma.chi,2) ) T(0,.9999)
iota.alpha[i] <- pow(iota.m[i], alpha)
chi.alpha[i]  <- pow(chi.m[i], alpha)

## Data model
## ----------
detJ.a.inv[i] <-  (Z[i] == 1) *  1   +
(Z[i] == 2) * ( 1/(1 - iota.m[i])  ) +
(Z[i] == 3) * ( 1/(1 - chi.m[i])  )
detJ.w.inv[i] <-  (Z[i] == 1) * ( 1/ (1 - min(a[i],.9999)) ) +
(Z[i] == 2) * ( (1 - iota.m[i])/( (1 - iota.alpha[i])*(1 - iota.m[i] - a[i]) ) ) +
(Z[i] == 3) * ( (1 - chi.m[i]) /( (1 - chi.alpha[i] )*(1 - chi.m[i]  - a[i]) ) ) 

g.inv.tau[i] <-   (Z[i] == 1) * (1 -  a[i] ) +           
(Z[i] == 2) * (1 - a[i]/(1 - iota.m[i]) )  +
(Z[i] == 3) * (1 - a[i]/(1 - chi.m[i] ) ) 
g.inv.nu[i]  <-   (Z[i] == 1) * ( w[i]*(1/(1 - min(a[i], .9999))) ) +
(Z[i] == 2) * ( w[i]*(1/(1 - iota.m[i] - min(a[i], .9999)))*((1 - iota.m[i])/(1-iota.alpha[i])) - iota.alpha[i]/(1 - iota.alpha[i]) - (a[i]*iota.m[i]) / ((1 - iota.m[i] - min(a[i], .9999))*(1-iota.alpha[i])) ) +
(Z[i] == 3) * ( w[i]*(1/(1 - chi.m[i]  - min(a[i], .9999)))*((1 - chi.m[i]) /(1-chi.alpha[i]) ) - chi.alpha[i] /(1 - chi.alpha[i])  - (a[i]*chi.m[i])  / ((1 - chi.m[i]  - min(a[i], .9999))*(1-chi.alpha[i])) )

g.inv.tau.scaled[i] = (g.inv.tau[i] - mu.tau[i]) / sigma.tau
b.tau[i]            = (  1    - mu.tau[i]) / sigma.tau
a.tau[i]            = (  0    - mu.tau[i]) / sigma.tau
g.inv.nu.scaled[i]  = (g.inv.nu[i]  - mu.nu[i])  / sigma.nu
b.nu[i]             = (  1    - mu.nu[i])  / sigma.nu
a.nu[i]             = (  0    - mu.nu[i])  / sigma.nu

k.tau[i] <-  pnorm(b.tau[i], 0, 1) - pnorm(a.tau[i], 0, 1) 
k.nu[i]  <-  pnorm(b.nu[i] , 0, 1) - pnorm(a.nu[i] , 0, 1)  

f.a[i]   <- ( (0 <= g.inv.tau[i] && g.inv.tau[i] <=1) * (dnorm(g.inv.tau.scaled[i], 0, 1) ) )/M
f.w[i]   <- ( (0 <= g.inv.nu[i]  && g.inv.nu[i]  <=1) * (dnorm(g.inv.nu.scaled[i] , 0, 1) ) )/M

ones.w[i] ~ dbern(f.a[i])
ones.a[i] ~ dbern(f.w[i])

}
} 
"}


rn_no_scaled_sep.vd <- function()
{
  "data{
  for(i in 1:n){
  ones.w[i] <- 1
  ones.a[i] <- 1
  }
}
model{
## Constants
## ---------
M	  = 10             ## auxiliary variable for the zero trick (makes the - loglikelihood > 0)
k         = .7
a.alpha	  = 1
b.alpha	  = 1
a.sigma	  = 1
b.sigma	  = 1
for(i in 1:3){
    psi.p[i] <- ifelse(i == 1, 1, .5)
  }
  for (i in 1:3){
    psi.i[i] ~ dbin(psi.p[i],1)
    psi.v[i] ~ dgamma(1,1)
    psi[i] <- psi.i[i]*psi.v[i]
  }

  for(i in 1:3){
    pi[i]	         <-  psi[i]/sum(psi[1:3]) ## mixing probabilities
  }
  pi1 <- pi[1]
  pi2 <- pi[1] + pi[2]  
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
beta.tau	 ~ dmnorm.vcov(mu.beta.tau, sigma.beta.tau)  ## linear coefficients of expectation of turnout
beta.nu	         ~ dmnorm.vcov(mu.beta.nu, sigma.beta.nu)    ## linear coefficients of expectation of votes for the winner
alpha	         ~ dgamma(a.alpha, b.alpha)

sigma.iota.m	 ~ dgamma(a.sigma, b.sigma)
sigma.chi        = 0.075
sigma.nu	 ~ dgamma(a.sigma, b.sigma)
sigma.tau	 ~ dgamma(a.sigma, b.sigma)

## sigma.iota.m	 ~ dunif(0,.5)
## sigma.chi        = 0.075
## sigma.nu	 ~ dunif(0,.5)
## sigma.tau	 ~ dunif(0,.5)

mu.iota.m		~ dunif(0,k)
mu.chi.m		~ dunif(k,1)

for(i in 1:n){
## linear transformation of the parameters of tau and nu
mu.tau[i]   = inprod(beta.tau, Xa[i,])
mu.nu[i]    = inprod(beta.nu , Xw[i,])

## Priors
## ------        
Zn[i] ~ dunif(0,1)
  Zn1[i] <- ifelse(Zn[i] <= pi1, 0, 1)
  Zn2[i] <- ifelse(Zn[i] > pi2, 1 , 0)
  Z[i]  <-  Zn1[i] + Zn2[i] + 1

iota.m[i]     ~  dnorm(mu.iota.m , 1/pow(sigma.iota.m,2)) T(0,.9999)     
chi.m[i]      ~  dnorm(mu.chi.m  , 1/pow(sigma.chi,2) ) T(0,.9999)
iota.alpha[i] <- pow(iota.m[i], alpha)
chi.alpha[i]  <- pow(chi.m[i], alpha)

## Data model
## ----------
detJ.a.inv[i] <-  (Z[i] == 1) *  1   +
(Z[i] == 2) * ( 1/(1 - iota.m[i])  ) +
(Z[i] == 3) * ( 1/(1 - chi.m[i])  )
detJ.w.inv[i] <-  (Z[i] == 1) * ( 1/ (1 - min(a[i],.9999)) ) +
(Z[i] == 2) * ( (1 - iota.m[i])/( (1 - iota.alpha[i])*(1 - iota.m[i] - a[i]) ) ) +
(Z[i] == 3) * ( (1 - chi.m[i]) /( (1 - chi.alpha[i] )*(1 - chi.m[i]  - a[i]) ) ) 

g.inv.tau[i] <-   (Z[i] == 1) * (1 -  a[i] ) +           
(Z[i] == 2) * (1 - a[i]/(1 - iota.m[i]) )  +
(Z[i] == 3) * (1 - a[i]/(1 - chi.m[i] ) ) 
g.inv.nu[i]  <-   (Z[i] == 1) * ( w[i]*(1/(1 - min(a[i], .9999))) ) +
(Z[i] == 2) * ( w[i]*(1/(1 - iota.m[i] - min(a[i], .9999)))*((1 - iota.m[i])/(1-iota.alpha[i])) - iota.alpha[i]/(1 - iota.alpha[i]) - (a[i]*iota.m[i]) / ((1 - iota.m[i] - min(a[i], .9999))*(1-iota.alpha[i])) ) +
(Z[i] == 3) * ( w[i]*(1/(1 - chi.m[i]  - min(a[i], .9999)))*((1 - chi.m[i]) /(1-chi.alpha[i]) ) - chi.alpha[i] /(1 - chi.alpha[i])  - (a[i]*chi.m[i])  / ((1 - chi.m[i]  - min(a[i], .9999))*(1-chi.alpha[i])) )

g.inv.tau.scaled[i] = (g.inv.tau[i] - mu.tau[i]) / sigma.tau
b.tau[i]            = (  1    - mu.tau[i]) / sigma.tau
a.tau[i]            = (  0    - mu.tau[i]) / sigma.tau
g.inv.nu.scaled[i]  = (g.inv.nu[i]  - mu.nu[i])  / sigma.nu
b.nu[i]             = (  1    - mu.nu[i])  / sigma.nu
a.nu[i]             = (  0    - mu.nu[i])  / sigma.nu

k.tau[i] <-  pnorm(b.tau[i], 0, 1) - pnorm(a.tau[i], 0, 1) 
k.nu[i]  <-  pnorm(b.nu[i] , 0, 1) - pnorm(a.nu[i] , 0, 1)  

f.a[i]   <- ( (0 <= g.inv.tau[i] && g.inv.tau[i] <=1) * (  (dnorm(g.inv.tau.scaled[i], 0, 1)) / ( sigma.tau * k.tau[i] )  ) ) /M
f.w[i]   <- ( (0 <= g.inv.nu[i]  && g.inv.nu[i]  <=1) * (  (dnorm(g.inv.nu.scaled[i] , 0, 1)) / ( sigma.nu  * k.nu[i]  )  ) ) /M

ones.w[i] ~ dbern(f.a[i])
ones.a[i] ~ dbern(f.w[i])



}
} 
"}


rn_no_alpha_sep.vd <- function()
{
  model = "
  data{
  for(i in 1:n){
  ones.w[i] <- 1
  ones.a[i] <- 1
  }
  }
  model{
  ## Constants
  ## ---------
  M	  = 10             ## auxiliary variable for the zero trick (makes the - loglikelihood > 0)
  k         = .7
  a.alpha	  = 1
  b.alpha	  = 1
  a.sigma	  = 1
  b.sigma	  = 1
  for(i in 1:3){
    psi.p[i] <- ifelse(i == 1, 1, .5)
  }
  for (i in 1:3){
    psi.i[i] ~ dbin(psi.p[i],1)
    psi.v[i] ~ dgamma(1,1)
    psi[i] <- psi.i[i]*psi.v[i]
  }

  for(i in 1:3){
    pi[i]	         <-  psi[i]/sum(psi[1:3]) ## mixing probabilities
  }
  pi1 <- pi[1]
  pi2 <- pi[1] + pi[2]  
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
  
  beta.tau	 ~ dmnorm.vcov(mu.beta.tau, sigma.beta.tau)  ## linear coefficients of expectation of turnout
  beta.nu	         ~ dmnorm.vcov(mu.beta.nu, sigma.beta.nu)    ## linear coefficients of expectation of votes for the winner
  sigma.nu	 ~ dgamma(a.sigma, b.sigma)
  sigma.tau	 ~ dgamma(a.sigma, b.sigma)
  
  
  mu.iota.m	~ dunif(0,k)
  mu.iota.s	~ dunif(0,k)
  mu.chi.m	~ dunif(k,1)
  mu.chi.s	~ dunif(k,1)
  sigma.iota.m	~ dgamma(a.sigma, b.sigma)
  sigma.iota.s    ~ dgamma(a.sigma, b.sigma)
  sigma.chi.m     = 0.075
  sigma.chi.s     = 0.075
  
  
  for(i in 1:n){
  ## linear transformation of the parameters of tau and nu
  mu.tau[i]   = inprod(beta.tau, Xa[i,])
  mu.nu[i]    = inprod(beta.nu , Xw[i,])
  
  ## Priors
  ## ------        
  Zn[i] ~ dunif(0,1)
  Zn1[i] <- ifelse(Zn[i] <= pi1, 0, 1)
  Zn2[i] <- ifelse(Zn[i] > pi2, 1 , 0)
  Z[i]  <-  Zn1[i] + Zn2[i] + 1
  
  iota.m[i]  ~  dnorm(mu.iota.m , 1/pow(sigma.iota.m,2)) T(0,.9999)            ## truncated normal
  iota.s[i]  ~  dnorm(mu.iota.s , 1/pow(sigma.iota.s,2)) T(0,.9999)            ## truncated normal
  chi.m[i]   ~  dnorm(mu.chi.m  , 1/pow(sigma.chi.m,2) ) T(0,.9999)
  chi.s[i]   ~  dnorm(mu.chi.s  , 1/pow(sigma.chi.s,2) ) T(0,.9999)
  
  ## Data model
  ## ----------
  tau[i] <-   (Z[i] == 1) * (1 -  min(a[i],.9999) ) +           
  (Z[i] == 2) * (1 - (min(a[i], .9999)/(1 - iota.m[i])) ) +
  (Z[i] == 3) * (1 - (min(a[i], .9999)/(1 - chi.m[i] )) ) 
  nu[i]  <-   (Z[i] == 1) * ( w[i]*(1/(1 - min(a[i], .9999))) ) +
  (Z[i] == 2) * ( w[i]*(1/(1 - iota.m[i] - min(a[i], .9999)))*((1 - iota.m[i])/(1-iota.s[i])) - iota.s[i]/(1 - iota.s[i]) - (a[i]*iota.m[i]) / ((1 - iota.m[i] - min(a[i], .9999))*(1-iota.s[i])) ) +
  (Z[i] == 3) * ( w[i]*(1/(1 - chi.m[i]  - min(a[i], .9999)))*((1 - chi.m[i]) /(1-chi.s[i]) ) - chi.s[i] /(1 - chi.s[i])  - (a[i]*chi.m[i])  / ((1 - chi.m[i]  - min(a[i], .9999))*(1-chi.s[i])) )
  
  ## k.tau[i] <-  pnorm(1, mu.tau[i], 1/(sigma.tau^2)) - pnorm(0, mu.tau[i], 1/(sigma.tau^2)) 
  ## k.nu[i]  <-  pnorm(1, mu.nu[i],  1/(sigma.nu^2))  - pnorm(0, mu.nu[i],  1/(sigma.nu^2))  
  
  ## p.tau[i] <- (0 <= tau[i] && tau[i] <=1) * dnorm(tau[i], mu.tau[i], 1/(sigma.tau^2)) / k.tau[i]
  ## p.nu[i]  <- (0 <= nu[i]  && nu[i]  <=1) * dnorm(nu[i] , mu.nu[i] , 1/(sigma.nu^2) ) / k.nu[i]
  ## p[i]     <- ( p.tau[i] * p.nu[i] ) / M
  
  tau.scaled[i] = (tau[i] - mu.tau[i]) / sigma.tau
  b.tau[i]      = (  1    - mu.tau[i]) / sigma.tau
  a.tau[i]      = (  0    - mu.tau[i]) / sigma.tau
  nu.scaled[i]  = (nu[i]  - mu.nu[i]) / sigma.nu
  b.nu[i]       = (  1    - mu.nu[i]) / sigma.nu
  a.nu[i]       = (  0    - mu.nu[i]) / sigma.nu
  
  k.tau[i] <-  pnorm(b.tau[i], 0, 1) - pnorm(a.tau[i], 0, 1) 
  k.nu[i]  <-  pnorm(b.nu[i] , 0, 1) - pnorm(a.nu[i] , 0, 1)  
  
  p.tau[i] <- ( (0 <= tau[i] && tau[i] <=1) * (  dnorm(tau.scaled[i], 0, 1) / ( sigma.tau * k.tau[i] )  ) ) / M
  p.nu[i]  <- ( (0 <= nu[i]  && nu[i]  <=1) * (  dnorm(nu.scaled[i] , 0, 1) / ( sigma.nu  * k.nu[i]  )  ) ) / M
  
  ones.w[i] ~ dbern(p.tau[i])
  ones.a[i] ~ dbern(p.nu[i])
  
  
  }
  
  
  } "
   invisible(model)
}



