
get_split_factors <- function(data, GrpID){
    
    # some quick input checks and define splitting factor
    if(is.character(GrpID)){
        if(!is.factor(data[,GrpID])){
            data[,GrpID] <- factor(data[,GrpID])
        }
        split_factor <- data[,GrpID]
    }else{
        assert_is_numeric(GrpID)
        if(sum(GrpID) > nrow(data)){
            stop("Group Indices exceed data rows")
        }
        if(sum(GrpID) < nrow(data)){
            warning("Group Indices are less than number of rows. \nAnalysis will omit some rows")
        }
        split_factor <- rep(1:length(GrpID), GrpID)
    }
    
    return(split_factor)
}

generate_sigma_matrix <- function(X, beta, gamma){
    
    sigma_matrix <- matrix(0, nrow = 3, ncol = 3)
    
    k <- nrow(X)
    df <- k - 1
    
    beta_hat <- mean(beta)
    gamma_hat <- mean(gamma)
    
    # fill covariance matrix (triangular)
    sigma_matrix[1,1] <- sum((X$mesor - mean(X$mesor))^2)/df
    sigma_matrix[1,2] <- sum((X$mesor - mean(X$mesor)) %*% (beta - mean(beta)))/df
    sigma_matrix[1,3] <- sum((X$mesor - mean(X$mesor)) %*% (gamma - mean(gamma)))/df
    sigma_matrix[2,2] <- sum((beta - mean(beta))^2)/df
    sigma_matrix[2,3] <- sum((beta - mean(beta)) %*% (gamma - mean(gamma)))/df
    sigma_matrix[3,3] <- sum((gamma - mean(gamma))^2)/df
    
    # mirror lower matrix
    sigma_matrix[lower.tri(sigma_matrix)] <- sigma_matrix[upper.tri(sigma_matrix)]
    
    return(sigma_matrix)
}