
#' @export
pmc = function(x , alpha, perd, ...)
{
    UseMethod("pmc")
}

pmc.test <- function(data, alpha = 0.05, perd = 24){
    
    if(!all(c("PR","Mesor", "Amp", "PHI") %in% colnames(data))){
        stop("data does not contain required variables ('PR','Mesor', 'Amp', 'PHI')")
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

#' @importFrom assertive assert_is_character
#' @export
pmc.default <- function(data, alpha = 0.05, perd = 24, GrpID = NA){
    
    if(!all(c("PR","Mesor", "Amp", "PHI") %in% colnames(data))){
        stop("data does not contain required variables ('PR','Mesor', 'Amp', 'PHI')")
    }
    
    if(is.null(GrpID) | any(is.na(GrpID))){
        out <- pmc_internal(data, alpha, perd)    
    }else{
        
        # some quick input checks and define splitting factor
        if(is.character(GrpID)){
            split_factor <- data[,GrpID]
            rownames <- unique(data[,GrpID])
        }else{
            assert_is_numeric(GrpID)
            if(sum(GrpID) > nrow(data)){
                stop("Group Indices exceed data rows")
            }
            if(sum(GrpID) < nrow(data)){
                warning("Group Indices are less than number of rows. \nAnalysis will omit some rows")
            }
            split_factor <- rep(1:length(GrpID), GrpID)
            if(!is.null(names(idx))){
                rownames <- names(idx)
            }else{
                rownames <- NULL
            }
        }
        
        # loop over each group and bind results together
        out <- do.call('rbind',
                       lapply(split(data, split_factor), function(x){
                           pmc_internal(x, alpha, perd)
                       })
        )
        row.names(out) <- rownames
    }
    
    
    return(out)
}




pmc_internal <- function(X, alpha = 0.05, perd = 24){
    
    # create output data.frame
    out <- data.frame(
        k = integer(),
        P.R. = integer(),
        P = numeric(),
        Mesor = numeric(),
        Mesor.CI = numeric(),
        Amplitude =  numeric(),
        Amplitude.CI.Lower =  numeric(),
        Amplitude.CI.Upper =  numeric(),
        Acrophase =  numeric(),
        Acrophase.CI.Lower =  numeric(),
        Acrophase.CI.Upper =  numeric(),
        Period =  numeric()
    )
    
    sigma_matrix <- matrix(0, nrow = 3, ncol = 3)
    
    k <- nrow(X)
    df <- k - 1
    
    beta= X$Amp*cos(X$PHI * pi/180)
    gamma=-X$Amp*sin(X$PHI * pi/180)
    
    beta_hat <- mean(beta)
    gamma_hat <- mean(gamma)
    
    # fill covariance matrix (triangular)
    sigma_matrix[1,1] <- sum((X$Mesor - mean(X$Mesor))^2)/df
    sigma_matrix[1,2] <- sum((X$Mesor - mean(X$Mesor)) %*% (beta - mean(beta)))/df
    sigma_matrix[1,3] <- sum((X$Mesor - mean(X$Mesor)) %*% (gamma - mean(gamma)))/df
    sigma_matrix[2,2] <- sum((beta - mean(beta))^2)/df
    sigma_matrix[2,3] <- sum((beta - mean(beta)) %*% (gamma - mean(gamma)))/df
    sigma_matrix[3,3] <- sum((gamma - mean(gamma))^2)/df
    
    # mirror lower matrix
    sigma_matrix[lower.tri(sigma_matrix)] <- sigma_matrix[upper.tri(sigma_matrix)]
    
    # Percentage Rhythm
    PR <- mean(X$PR)
    
    # Population Mesor
    Mesor <- mean(X$Mesor)
    
    # Population Acrophase
    Phi <- -phsrd(beta, gamma)
    
    # Population Amplitude
    Amp <- module(beta, gamma)
    
    # t-value
    tval <- qt(1-alpha/2, df)
    
    # Mesor CI
    cim <- tval * sqrt(sigma_matrix[1,1]/k)
    
    # Amplitude CI
    c22 <- (sigma_matrix[2,2] * beta_hat^2 + 2*sigma_matrix[2,3]*beta_hat*gamma_hat + sigma_matrix[3,3] * gamma_hat^2) / (k * Amp^2)
    cia <- tval*sqrt(c22)
    
    # Acrophase CI
    c23 <- (-(sigma_matrix[2,2] - sigma_matrix[3,3]) * (beta_hat*gamma_hat) + sigma_matrix[2,3]*(beta_hat^2 - gamma_hat^2)) / (k * Amp^2)
    c33 <- (sigma_matrix[2,2] * gamma_hat^2 - 2 * sigma_matrix[2,3] * beta_hat * gamma_hat + sigma_matrix[3,3]*beta_hat^2) / (k * Amp^2)
    
    an1 <- Amp^2 - (c22*c33 - c23^2) * (tval^2)/c33
    
    if(an1 < 0){
        Phi1 <- 0
        Phi2 <- 0
    }else{
        Phi1 <- Phi + atan((c23 * tval^2 + tval*sqrt(c33) * sqrt(an1))/(Amp^2 - c22*tval^2)) * 180/pi
        Phi2 <- Phi + atan((c23 * tval^2 - tval*sqrt(c33) * sqrt(an1))/(Amp^2 - c22*tval^2)) * 180/pi
        Phi1 <- phase(Phi1)
        Phi2 <- phase(Phi2)
    }
    
    
    # p-value calculation
    r <- sigma_matrix[2,3]/sqrt(sigma_matrix[2,2]*sigma_matrix[3,3])
    fval <- k*(k-2)/(2*(k-1) * (1-r^2)) *
        (beta_hat^2/sigma_matrix[2,2]
         -2*r*beta_hat*gamma_hat/sqrt(sigma_matrix[2,2]*sigma_matrix[3,3])
         +gamma_hat^2/sigma_matrix[3,3]
        )
    p <- pf(fval, df1 = 2, df2 = k - 2, lower.tail = FALSE)
    
    
    
    # This could be much more concise but doing it like this to be explicit
    out[1,"k"] <- as.integer(k)
    out[1,"P.R."] <- as.integer(round(mean(X$PR)))
    out[1,"P"] <- p
    out[1,"Mesor"] <- round(Mesor, 3)
    out[1,"Mesor.CI"] <- round(cim, 3)
    out[1,"Amplitude"] <- round(Amp, 3)
    out[1,"Amplitude.CI.Lower"] <- ifelse(p < alpha, round(Amp - cia, 3), 0)
    out[1,"Amplitude.CI.Upper"] <- ifelse(p < alpha, round(Amp + cia, 3), 0)
    out[1,"Acrophase"] <- round(Phi, 3)
    out[1,"Acrophase.CI.Lower"] <- ifelse(p < alpha, round(Phi1 - 0.5, 3), 0)
    out[1,"Acrophase.CI.Upper"] <- ifelse(p < alpha, round(Phi2 - 0.5, 3), 0)
    out[1,"Period"] <- perd
    
    return(out)
    
}


# get polar coordinates from rectangular
phsrd <- function(beta, gamma){
    
    m_beta <- mean(beta)
    m_gamma <- mean(gamma)
    
    tmp <- 180/pi * atan(m_gamma/m_beta)
    
    # correction if in different quadrants
    if(m_beta < 0 | m_gamma < 0){
        if(m_beta > 0){
            # If in Quadrant IV
            tmp <- tmp + 360
        }else{
            # If in Quandrant II or III
            tmp <- tmp + 180
        }
    }
    return(tmp)
}

module <- function(beta, gamma){
    c=max(c(abs(mean(beta)),abs(mean(gamma))))
    d=min(c(abs(mean(beta)),abs(mean(gamma))))
    q=d/c
    amp=c*sqrt(1+q*q)
    return(amp)
}

phase <- function(Phi){
    if(Phi > 0){
        Phi <- Phi - 360
    }
    if(Phi < -360){
        Phi <- Phi + 360
    }
    return(Phi)
}






