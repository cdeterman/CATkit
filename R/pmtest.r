

#' @export
pmtest <- function(X, alpha = 0.05, GrpID = NA, SubGrp = NA){

    if(is.na(GrpID) | is.null(GrpID)){
        stop("GrpID must be defined")
    }
    
    if(!is.na(SubGrp)){
        if(!all(SubGrp %in% X[,GrpID])){
            stop("Some SubGrp ids not found in GrpID column")
        }
    }
    
    # filter to any subgroups
    if(!is.na(SubGrp)){
        X <- X[X[,GrpID] %in% SubGrp,]    
    }
    
    out <- pmtest_internal(X, alpha = alpha, GrpID = GrpID)
    
    return(out)
}

# k - vector of group sizes
# K - total sample sizes
# M - list of group Mesors
# A - list of group Amplitudes
# m - number of groups
# sigma_m - sigma_matrix[1,1]

pmtest_internal <- function(X, alpha = 0.05, GrpID = NA){
    
    vars <- tolower(colnames(X))
    idx <- which(vars %in% c("pr", "mesor", "amp", "phi"))
    colnames(X)[idx] <- vars[idx]

	# total sample size
	K <- nrow(X)
	m <- length(unique(X[,GrpID]))
	k <- table(X[,GrpID])
	
	# split populations
	split_factor <- get_split_factors(X, GrpID)
	rownames <- names(split(X, split_factor))   

	X_list <- split(X, split_factor)
	
	beta_star <- X$amp*cos(X$phi * pi/180)
	gamma_star <- -X$amp*sin(X$phi * pi/180)
	
    beta <- sapply(X_list, function(x) x$amp*cos(x$phi * pi/180))
    beta_hat <- sapply(beta, mean)
    gamma <- sapply(X_list, function(x) -x$amp*sin(x$phi * pi/180))
    gamma_hat <- sapply(gamma, mean)
    
    # Population Mesors
    M <- sapply(X_list, function(x) mean(x$mesor))
	
    # Population Acrophases
    phi <- mapply(function(b, g) -phsrd(b, g), beta, gamma)
    phi_star <- -phsrd(beta_star, gamma_star)
	
    # Population Amplitudes
    A <- mapply(function(b, g) module(b, g), beta, gamma)
	A_star <- module(beta_star, gamma_star)
    
	# get sigma_matrix for each population
    sigma_matrices <- lapply(names(X_list), function(i) {generate_sigma_matrix(X_list[[i]], beta[[i]], gamma[[i]])})
	
	# Eq. 66
	#sigma_matrix_hat <- sum((k - 1) * sigma_matrices / (K - m))
    sigma_matrix_hat <- Reduce(`+`, lapply(1:m, function(x) (k[x] - 1) * sigma_matrices[[x]]/(K - m)))
	
	M_bar <- sum(k * M)/K
	beta_bar <- sum(k * beta_hat)/K
	gamma_bar <- sum(k * gamma_hat)/K
	
	print(mean(X$mesor))
	print(A_star)
	print(phi_star)

	# Eq. 67
	t_mat <- matrix(0, nrow = 3, ncol = 3)
	t_mat[1,1] <- sum(k * (M - M_bar)^2)
	t_mat[1,2] <- sum(k * (M - M_bar) * (beta_hat - beta_bar))
	t_mat[1,3] <- sum(k * (M - M_bar) * (gamma_hat - gamma_bar))
	t_mat[2,2] <- sum(k * (beta_hat - beta_bar)^2)
	t_mat[2,3] <- sum(k * (beta_hat - beta_bar) * (gamma_hat - gamma_bar))
	t_mat[3,3] <- sum(k * (gamma_hat - gamma_bar)^2)
	
	t_mat[lower.tri(t_mat)] <- t_mat[upper.tri(t_mat)]

	t1 <- t_mat[2:3, 2:3]

	# test for equal population mesors
	# Eq. 68
	m_fval <- t_mat[1,1]/((m - 1) * sigma_matrix_hat[1,1])
	m_p <- pf(m_fval, df1 = m-1, df2 = K-m, lower.tail = FALSE)


	# multivariate test of rhythm parameters
	# Eq. 69
	J <- (sigma_matrix_hat[2:3, 2:3]^-1) %*% t1
	print(J)
	# Eq. 70
	D <- det(diag(2) + J/(K-m))
	print("D")
	print(D)

	if(m > 2){
		# Eq. 71
	    rhythm_fval <- (K-m-1)/(m-1) * (sqrt(D) - 1)
		rhythm_p <- pf(rhythm_fval, df1 = 2*m - 2, 2*K - 2*m - 2, lower.tail = FALSE)
	}else{
		# Eq. 72
		rhythm_fval <- (sum(k) - 3)/2 * (D - 1)
		rhythm_p <- pf(rhythm_fval, df1 = 2, df2 = sum(k) - 3, lower.tail = FALSE)
	}

	# if multivariate test significant?
	# test for equivalent population acrophases (approximate)

	# tan(2*phi_tilda) = sum(mapply(`*`, k, A^2) * sin(2*phi * pi/180)) / sum(mapply(`*`, k, A^2) * cos(2*phi * pi/180))
	num <- sum(k * A^2 * sin(2*phi * pi/180))
	den <- sum(k * A^2 * cos(2*phi * pi/180))
	#num <- Reduce(`+`, Reduce(function(x, y) Map(`*`, x, y), list(k, A^2, sin(2*phi * pi/180))))
	#den <- Reduce(`+`, Reduce(function(x, y) Map(`*`, x, y), list(k, A^2, cos(2*phi * pi/180))))
	phi_tilda <- atan(num/den)/2 * pi/180

	# Eq. 74
	num <- sum(k * A^2 * sin((phi * pi/180) - phi_tilda)^2)/(m - 1)
	den <- sigma_matrix_hat[2,2] * sin(phi_tilda)^2 + 2 * sigma_matrix_hat[2,3] * cos(phi_tilda) * sin(phi_tilda) + sigma_matrix_hat[3,3] * cos(phi_tilda)^2
	phi_fval <- num/den
	phi_p <- pf(phi_fval, df1 = m-1, df2 = K-m, lower.tail = FALSE)
	
	# test for equivalent population amplitudes
	# NO KNOWN SOLUTION
	# ASSUMPTION - acrophases do not differ greatly
	# Eq. 75
	# num <- sum(mapply(`*`, k, (A - A_star)^2))/(m-1)
	# phi_star <- as.integer(phi_star)
	num <- sum(k * (A - A_star)^2)/(m-1)
	den <- (sigma_matrix_hat[2,2] * (cos(phi_star * pi/180)^2)) - (2 * sigma_matrix_hat[2,3] * cos(phi_star * pi/180)*sin(phi_star * pi/180)) + (sigma_matrix_hat[3,3]*(sin(phi_star * pi/180)^2))
	amp_fval <- num/den
	amp_p <- pf(amp_fval, df1 = m-1, df2 = K-m, lower.tail = FALSE)
	
	# combined params
	print(M_bar)
	print(A_star)
	print(phi_star)
	
	out <- data.frame(F = c(m_fval, amp_fval, phi_fval, rhythm_fval),
	                  P = c(m_p, amp_p, phi_p, rhythm_p))
	row.names(out) <- c("Mesor", "Amplitude", "Acrophase", "(A, phi)")
	# out <- data.frame(Mesor = m_p, Rhythm = rhythm_p, Acrophase = phi_p, Amplitude = amp_p)
	return(out)
}

