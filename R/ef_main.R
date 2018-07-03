


get_Z <- function(samples)
{
    ## replace matrix Z (sample size) x (number of interation) to a single column matrix with z.hat
    ## the estimated cluster of Zi. Zi is classified in the cluster it has highest estimated
    ## posterior probability to belong to
    for (i in 1:length(samples))
    {
        z            = samples[[i]][,base::grepl(pattern='Z.[0-9]*.', x=colnames(samples[[i]]))]
        k.hat        =   base::apply(z, 2, function(zi) which.max(c(sum(zi==1), sum(zi==2), sum(zi==3)) ) )
        piZi         = t(base::apply(z, 2, function(zi) c(sum(zi==1), sum(zi==2), sum(zi==3))/length(zi))) 
        colnames(piZi )= c("pi[Zi1]", "pi[Zi2]", "pi[Zi3]")
        samp         = samples[[i]][,!base::grepl(pattern='Z.[0-9]*.', x=colnames(samples[[i]]))]
        samp         = list(parameters=coda::as.mcmc(samp), k.hat=k.hat, piZi=piZi)
        samples[[i]] = samp
    }
    return(samples)
}

create_list <- function(samples)
{
    for (i in 1:length(samples))
    {
        samples[[i]] = list(parameters=samples[[i]])
    }
    return(samples)
}

ef_get_parameters_to_monitor <- function(model, all=FALSE)
{
    if(model == 'rn')           parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha")
    if(model == 'rn_no_scaled') parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha")
    if(model == 'normal')       parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha")
    if(model == 'rn_no_alpha')  parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "mu.iota.s", "mu.chi.s", "sigma.iota.s")

    if(model == 'rn_sep')           parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha")
    if(model == 'rn_no_scaled_sep') parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha")
    if(model == 'normal_sep')       parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha")
    if(model == 'rn_no_alpha_sep')  parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "mu.iota.s", "mu.chi.s", "sigma.iota.s")

    if(model == 'bl')           parameters = c("pi", 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "mu.iota.s", "mu.chi.s")
    
    if(model == 'bl.vd')        parameters = c("pi", 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "mu.iota.s", "mu.chi.s","psi.i")
    
    if(model == 'rn.vd')           parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha", "psi.i")
    if(model == 'rn_no_scaled.vd') parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha", "psi.i")
    if(model == 'normal.vd')       parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha", "psi.i")
    if(model == 'rn_no_alpha.vd')  parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "mu.iota.s", "mu.chi.s", "sigma.iota.s", "psi.i")
    
    if(model == 'rn_sep.vd')           parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha", "psi.i")
    if(model == 'rn_no_scaled_sep.vd') parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha", "psi.i")
    if(model == 'normal_sep.vd')       parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha", "psi.i")
    if(model == 'rn_no_alpha_sep.vd')  parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "mu.iota.s", "mu.chi.s", "sigma.iota.s", "psi.i")

    if(all) parameters = c(parameters, 'Z')

    return(parameters)
}

get_model <- function(model)
{
    if (model == 'rn')           return(rn())
    if (model == 'rn_no_scaled') return(rn_no_scaled())
    if (model == 'normal')       return(normal())
    if (model == 'rn_no_alpha')  return(rn_no_alpha())

    if (model == 'rn_sep')           return(rn_sep())
    if (model == 'rn_no_scaled_sep') return(rn_no_scaled_sep())
    if (model == 'normal_sep')       return(normal_sep())
    if (model == 'rn_no_alpha_sep')  return(rn_no_alpha_sep())

    if (model == 'bl')          return(bl())
  
    if (model == 'bl.vd')       return(bl.vd())
  
  if (model == 'rn.vd')           return(rn.vd())
  if (model == 'rn_no_scaled.vd') return(rn_no_scaled.vd())
  if (model == 'normal.vd')       return(normal.vd())
  if (model == 'rn_no_alpha.vd')  return(rn_no_alpha.vd())
  
  if (model == 'rn_sep.vd')           return(rn_sep.vd())
  if (model == 'rn_no_scaled_sep.vd') return(rn_no_scaled_sep.vd())
  if (model == 'normal_sep.vd')       return(normal_sep.vd())
  if (model == 'rn_no_alpha_sep.vd')  return(rn_no_alpha_sep.vd())

}

getRegMatrix <- function(func.call, data, weights, formula_number=1)
{
    args <- names(func.call)
    ## creating the dependent variable and the covariates matrix from the fomula 1
    f = paste0('formula', formula_number, sep='')
    idx.args  <- match(c(f,  "data", "weights"), args, 0L)
    func.call <- func.call[c(1L, idx.args)]
    names(func.call)[names(func.call)==f] = "formula"
    func.call$drop.unused.levels <- TRUE
    func.call[[1L]] <- quote(stats::model.frame)
    func.call[[3]] = quote(data)
    reg.matrix <- eval(func.call, parent.frame())
    ## response variable
    y   <- stats::model.response(reg.matrix, "numeric")
    ## weights
    w   <- as.vector(stats::model.weights(reg.matrix))
    if (!is.null(w) && !is.numeric(w)) stop("'weights' must be a numeric vector")
    offset <- as.vector(stats::model.offset(func.call))
    if (!is.null(offset)) {
        if (length(offset) != NROW(y)) 
            stop(gettextf("number of offsets is %d, should equal %d (number of observations)", length(offset), NROW(y)), domain = NA)
    }
    ## covariates
    mt1    <- attr(reg.matrix, "terms")
    if (stats::is.empty.model(mt1)) {
        x <- matrix(1, ncol=1,nrow=nrow(y))
        results <- list(coefficients = if (is.matrix(y)) matrix(, 0, 3) else numeric(), residuals = y, fitted.values = 0 * y, weights = w, rank = 0L, df.residual = if (!is.null(w)) sum(w != 0) else if (is.matrix(y)) nrow(y) else length(y))
        if (!is.null(offset)) {
            results$fitted.values <- offset
            results$residuals <- y - offset
        }
    } else {
        x <- stats::model.matrix(mt1, reg.matrix, contrasts)
    }
    return(list(y=y, X=x, w=w))
}


## {{{ docs }}}
#' Election Forensics Finite Mixture Model
#'
#' This function estimates a finite mixture model of election fraud
#'
#'
#' @param formula1 an object of the class \code{formula} as used in \code{\link{lm}}. The independent variable must the votes for the winner party or candidate, whose format must be either integers (number of votes) if the model used is the binomial model, or proportions (of votes) if the model used is the restricted normal (see \code{model})
#' @param formula2 an object of the class \code{formula} as used in \code{\link{lm}}. The independent variable must the abstention , whose format must be same as the independent variable in \code{formula1}
#' @param data a dara.frame with the independent variables (voters for the winner and abstention) and the covariates. If the independent variables are counts, it must contain also a variable \code{N} which must be an integer with the number of elegible voters
#' @param weights (not used)
#' @param mcmc a list containing \code{n.iter}, which is the number of iterations for the MCMC, \code{burn.in} for the burn-in period of the MCMC chain, \code{n.adapt} indicating the number of adaptative steps before the estimation (see \code{\link{rjags}})
#' @param model a string with either "bl", "rn", "rn_no_alpha". It indicates the forensics model of the data. "rn" stands for restricted normal model. "bl" stands binomial with logistic transformation of the probability parameter. These options indicate which model should be used as the distribution for the voters for the winner and abstention.
#' @param parameters a string vector with the names of the parameters to monitor. When \code{NULL} (default), it will monitor all the parameters, except the Z's. When \code{parameters='all'}, it will monitor all parameters, including Z, which is necessary to classify the observations as fraudulent cases or not.
#' @param na.action (not used)
#'
#' @return The function returns a nested list. The first element of the list is a \code{mcmc} object with the samples from the posterior distribution. The second element of the list is a list of summaries (HPD, Mean, etc)
#'
#' @export
## }}}
eforensics   <- function(formula1, formula2, data, weights, mcmc, model, parameters=NULL, na.action="exclude")
{

    ## error handling
    check_mcmc(mcmc)

    options(warn=-1)
    on.exit(options(warn=0))
    ## check if JAGS is installed
    ef_check_jags()
    
    ## ## construct the regression matrices (data.frames) based on the formula provided
    ## ## -----------------------------------------------------------------------------
    func.call <- match.call(expand.dots = FALSE)
    mat     = getRegMatrix(func.call, data, weights, formula_number=1)
    w       = mat$y
    Xw      = mat$X
    weightw = mat$w
    mat     = getRegMatrix(func.call, data, weights, formula_number=2)
    a       = mat$y
    Xa      = mat$X
    weighta = mat$w
    if(model == 'bl' | model == 'bl.vd'){
        data    = list(w = w, a = a, Xa = as.matrix(Xa), Xw = as.matrix(Xw), dxw = ncol(Xw), dxa = ncol(Xa), n = length(w), N = data$N)
    }else{
        data    = list(w = w, a = a, Xa = as.matrix(Xa), Xw = as.matrix(Xw), dxw = ncol(Xw), dxa = ncol(Xa), n = length(w))
    }

    ## get parameters to monitor
    ## -------------------------
    if(is.null(parameters)) parameters = ef_get_parameters_to_monitor(model)
    if(parameters[1] == 'all') parameters = ef_get_parameters_to_monitor(model, all=TRUE)

    ## get model
    ## ---------
    model.name = model
    model      = get_model(model.name)

    ## Debug/Monitoring message --------------------------
    msg <- paste0('\n','Burn-in: ', mcmc$burnin, '\n'); cat(msg)
    ## msg <- paste0('\n','Chains: ', mcmc$n.chains, '\n'); cat(msg)
    msg <- paste0('\n','Number of MCMC samples per chain: ', mcmc$n.iter, '\n'); cat(msg)
    msg <- paste0('\n','MCMC in progress ....', '\n'); cat(msg)
    ## ---------------------------------------------------

    ## MCMC
    ## ----
    time.init    = Sys.time()
    cat('\nCompiling the model...\n')     ; sim = rjags::jags.model(file=textConnection(model), data = data, n.adapt=mcmc$n.adapt, n.chain=mcmc$n.chains)
    cat('\nUpdating MCMC (burn-in) ...\n'); stats::update(sim, n.iter = mcmc$burn.in)
    cat('\nDrawing the samples...\n')     ; samples = rjags::coda.samples(model=sim, variable.names=parameters, n.iter=mcmc$n.iter)
    T.mcmc = Sys.time() - time.init

    if(!is.null(parameters) & "Z" %in% parameters)
        samples = get_Z(samples)
    else
        samples = create_list(samples)
    class(samples) = "eforensics"

    attr(samples, "formula.w") = formula1
    attr(samples, "formula.a") = formula2
    attr(samples, "model")     = model.name

    cat("\n\nEstimation Completed\n\n")
    return(samples)
}

## eforensics_par   <- function(formula1, formula2, data, weights, mcmc, model, parameters=NULL, na.action="exclude")
## {

##     ## check if JAGS is installed
##     ef_check_jags()
    
##     ## ## construct the regression matrices (data.frames) based on the formula provided
##     ## ## -----------------------------------------------------------------------------
##     func.call <- match.call(expand.dots = FALSE)
##     mat     = getRegMatrix(func.call, data, weights, formula_number=1)
##     w       = mat$y
##     Xw      = mat$X
##     weightw = mat$w
##     mat     = getRegMatrix(func.call, data, weights, formula_number=2)
##     a       = mat$y
##     Xa      = mat$X
##     weighta = mat$w
##     if(model == 'bl'){
##         data    = list(w = w, a = a, Xa = as.matrix(Xa), Xw = as.matrix(Xw), dxw = ncol(Xw), dxa = ncol(Xa), n = length(w), N = data$N)
##     }else{
##         data    = list(w = w, a = a, Xa = as.matrix(Xa), Xw = as.matrix(Xw), dxw = ncol(Xw), dxa = ncol(Xa), n = length(w))
##     }


##     ## get parameters to monitor
##     ## -------------------------
##     if(is.null(parameters)) parameters = ef_get_parameters_to_monitor(model)

##     ## get model
##     ## ---------
##     model = get_model(model)

##     ## Debug/Monitoring message --------------------------
##     msg <- paste0('\n','Burn-in: ', mcmc$burnin, '\n'); cat(msg)
##     ## msg <- paste0('\n','Chains: ', mcmc$n.chains, '\n'); cat(msg)
##     msg <- paste0('\n','Number of MCMC samples per chain: ', mcmc$n.iter, '\n'); cat(msg)
##     msg <- paste0('\n','MCMC in progress ....', '\n'); cat(msg)
##     ## ---------------------------------------------------

##     ## MCMC (parallel)
##     ## ----
##     cl <- parallel::makePSOCKcluster(min(parallel::detectCores()-1, mcmc$n.chains))
##     samples = dclone::jags.parfit(cl       = cl,
##                                   data     = data,
##                                   params   = parameters,
##                                   model    = textConnection(model),
##                                   n.adapt  = mcmc$n.adapt,
##                                   n.chains = mcmc$n.chains,
##                                   n.update = mcmc$burn.in,
##                                   n.iter   = mcmc$n.iter,
##                                   thin     = 1)
##     parallel::stopCluster(cl)


##     ## computing summary
##     ## -----------------
##     ## summary <- list(summary = summary(samples), HPD = coda::HPDinterval(samples))
##     ## results = list(samples=samples, stat=summary, time.elapsed=T.mcmc)
##     results = samples

##     cat("\n\nEstimation Completed\n\n")
##     return(results)
## }
