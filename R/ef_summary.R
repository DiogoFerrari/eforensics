
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

## {{{ docs }}}
#' Summary 
#'
#' This function return summaries of the posterior distribution estimated by the function \code{eforensics}
#'
#'
#' @param object the output of the function \code{eforensics}
#' @param ... join.chaing=TRUE can be used to provide summaries of chains after they are combined together 
#'
#' @export
## }}}
summary.eforensics <- function(object, ...)
{
    args = as.list(match.call())
    if('join.chains' %in% names(args)) {
        join.chains = as.logical(args$join.chains)
    }else{
        join.chains = FALSE
    }
        

    x = object
    if (join.chains) {
        samp = x %>% purrr::map(.x=., ~.x['parameters'][[1]])  %>% do.call(rbind,.)
        HPD = samp %>%
            coda::as.mcmc(.) %>%
            coda::HPDinterval(.) %>%
            data.frame(Parameter = row.names(.), ., row.names = 1:nrow(.)) %>%
            dplyr::rename(HPD.lower = lower, HPD.upper = upper)  

        samp = samp %>%
            coda::as.mcmc(.) %>%
            summary(.) %>%
            .[[1]] %>%
            data.frame(Parameter = row.names(.), ., row.names = 1:nrow(.))  %>%
            dplyr::select(Parameter, Mean, SD) %>%
            dplyr::full_join(., HPD , by=c("Parameter"))  
    }else{
        samp = list()
        for (i in 1:length(x))
        {
            HPD = x[[i]]$parameters %>%
                coda::HPDinterval(.) %>%
                data.frame(Parameter = row.names(.), ., row.names = 1:nrow(.)) %>%
                dplyr::rename(HPD.lower = lower, HPD.upper = upper)  
            tab = x[[i]]$parameters %>%
                summary(.) %>%
                .[[1]] %>%
                data.frame(Parameter = row.names(.), ., row.names = 1:nrow(.))  %>%
                dplyr::select(Parameter, Mean, SD) %>%
                dplyr::full_join(., HPD , by=c("Parameter"))  
            samp[[i]] = tab
        }
        names(samp) = paste0('Chain ', 1:length(x)) 
    }
    return(samp)
}

#' @export
print.eforensics <- function(x, ...)
{
    summary(x)
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
