Overview
========

Instalation
===========

**NOTE:** [JAGS](https://sourceforge.net/projects/mcmc-jags/) is
required to use the current version of this package

``` {.r .rundoc-block rundoc-language="R" rundoc-exports="code"}
# Install the development version (requires the package "devtools", so install it first if it is not installed already)
devtools::install_github("eforensics")
```

Usage
=====

In R, check `help(eforensics)`

``` {.r .rundoc-block rundoc-language="R" rundoc-exports="code"}
library(eforensics)
## model
## -----
model    = 'bl'

## simulating data and parameters
## ------------------------------
help(ef_simulateData)
sim_data = ef_simulateData(n=7000, nCov=1, model=model)
data     = sim_data$data

## mcmc parameters
## ---------------
mcmc    = list(burnin=1, n.adapt=10, n.iter=10)

## samples
## -------
help(eforensics)
samples    = eforensics(w ~.-a, a ~ .-w,  data=data, model=model, mcmc=mcmc)

summary(samples)

```
