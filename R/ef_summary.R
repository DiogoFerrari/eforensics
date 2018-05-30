
## {{{ docs }}}
#' Classify observations into fraud distributions
#'
#' This function classifies the data points into one of three distributions: incremental fraud, extreme fraud, and no fraud.
#'
#'
#' @param data a data frame with the data
#' @param samples the output of the function \code{eforensics}  
#'
#'
#' @export
## }}}
ef_classify <- function(data, samples){
    if (!"k.hat" %in% names(samples[[1]])) 
        stop("You need to monitor parameter Z as well in order to classify the data.")

    k.list = list()
    for (i in 1:length(samples))
    {
        k.list[[i]] = samples[[i]]$k.hat
    }
    k.list = k.list %>%
        do.call(cbind,.) %>%
        tibble::as_data_frame(.)  %>% 
        data.table::setnames(., paste0("chain.", 1:length(samples), sep='')) %>%
        base::apply(., 1, function(zi) which.max(c(sum(zi==1), sum(zi==2), sum(zi==3)) ) )
    return(cbind(data, Estimated.Distribution=k.list))
}

#' @export
summary.ef <- function(object, ...){
    samp = list()
    for (i in 1:length(object))
    {
        samp[[i]] = coda::as.mcmc(object[[i]]$parameters)
    }
    samp = coda::as.mcmc(samp)
    return(samp %>% summary(.))
}

## {{{ docs }}}
#' Get parameters
#'
#' This function returns the parameters of the model
#'
#'
#' @param samples the output of the function \code{eforensics} 
#'
#' @export
## }}}
ef_get_parameters <- function(samples)
{
    samp = list()
    for (i in 1:length(samples))
    {
        samp[[i]] = samples[[i]]$parameters
    }
    return(samp)   
}
