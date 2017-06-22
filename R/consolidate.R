
# function to collect output files from CATCosinor

#' @title Consolidate CATCosinor Results
#' @description This function takes the path originally passed to
#' the CATCosinor function and collects all the relevant results
#' for subsequent analysis in a \code{data.frame} object.
#' @param path The path where the CATCosinor .dat output files reside
#' @return A data.frame object
#' @export
consolidate <- function(path){
    files <- list.files(path, full.names = TRUE)
    files <- files[grepl('\\.dat', files)]
    
    cosdat <- as.data.frame(do.call('rbind', lapply(files, function(x){
        lines <- readLines(x)
        lines <- lines[c(1, seq(2, length(lines), 4))]
        tmp <- do.call('rbind', strsplit(lines[2:length(lines)], "\t"))
        colnames(tmp) <- strsplit(lines[1], "\t")[[1]]
        return(tmp)
    })))
    
    # format the numeric columns
    cosdat[,c(9:10,12:ncol(cosdat))] <- sapply(cosdat[,c(9:10,12:ncol(cosdat))], function(x) as.numeric(as.character(x)))
    
    return(cosdat)
}