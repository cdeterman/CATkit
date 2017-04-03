
#' @export
pmc = function(x , alpha, perd, ...)
{
    UseMethod("pmc")
}

#' @method default
#' @export
pmc.default <- function(data, alpha = 0.05, perd = 24){
    
    if(!all(c("PR","Mesor", "Amp", "Phi") %in% colnames(data))){
        stop("data does not contain required variables ('PR','Mesor', 'Amp', 'Phi')")
    }
    
    rytpar <- matrix(0, nrow = 6, ncol = 13)
    sig <- matrix(0, nrow = 3, ncol = 3)   
    
    k <- nrow(data)
    
    out <- cpp_bmpop(data, sig, rytpar, 
                     alpha, k, perd)
    
    result <- out$rytpar[1,]
    result <- result[-8]
    names(result) <- c(
        "Period","P.R.","P",
        "Amplitude", "Acrophase", "Amplitude CI Lower","Amplitude CI Upper",
        "Acrophase CI Lower","Acrophase CI Upper","Mesor","Mesor CI","k")
    
    return(result)   
}


#' @method pmc list
#' @export
pmc.list <- function(data, alpha = 0.05, perd = 24, GrpID = NA){
    
    if(is.na(GrpID)){
        warning("No group defined")
    }
    
    if(!is.na(GrpID)){
        result <- lapply(data, function(x) {
            lapply(split(x, x[,GrpID]), function(y){
                pmc.default(as.matrix(y[,c("PR","Mesor", "Amp", "Phi")]), alpha, perd)
            })
        })
        
        result <- lapply(result, function(x) do.call('rbind', x))
        
    }else{
        result <- lapply(data, function(x) {
            pmc.default(as.matrix(x[,c("PR","Mesor", "Amp", "Phi")]), alpha, perd)
        }) 
    }
    
    
    return(result)   
}

