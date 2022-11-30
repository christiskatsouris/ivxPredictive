###########################################
# R Script Details:
###########################################

# Script name: supWaldIVX.R

# Program aim: This R program implements the sup Wald IVX test for structural break detection 
# on the slopes of a linear predictive regression model.  

# written by: 

# Christis G. Katsouris (May 2021)
# Department of Economics
# University of Southampton
# Southampton, United Kingdom

#################################################################################################################
### Function 1: Simulate data pair (X,Y) under the null hypothesis of no breaks
#################################################################################################################

simulate_data_null_function <- function( N = N_size, beta1=beta1, beta2=beta2, c1=c1, c2=c2, rho = rho, gamma.x = gamma.x )
{# begin of function
  
  N   <- N_size
  pi0 <- pi0
  # Here we only check the model without intercept
  beta1   <- beta1
  beta2   <- beta2
  c1      <- c1
  c2      <- c2
  rho     <- rho
  gamma.x <- gamma.x
  
  p   <- 2
  mu.vector <- matrix(0, nrow = p+1, ncol = 1 ) 
  Sigma     <- matrix(0, nrow = p+1, ncol = p+1 )
  
  sigma_v1v2 <- (rho)*sqrt(1)*sqrt(1)
  
  Sigma[1,1] <- 1
  Sigma[1,2] <- 0.10
  Sigma[1,3] <- -0.29  
  Sigma[2,1] <- 0.10
  Sigma[2,2] <- 1
  Sigma[2,3] <- sigma_v1v2
  Sigma[3,1] <- -0.29
  Sigma[3,2] <- sigma_v1v2
  Sigma[3,3] <- 1
  
  # generate random error sequence from Multivariate Normal Distribution
  innov.e  <- rmvnorm( n = N, mean = mu.vector, sigma = Sigma ) 
  innov.u  <- as.matrix( innov.e[ ,1] )
  innov.v  <- as.matrix( innov.e[ ,2:3] )
  
  C      <- diag(p)
  C[1,1] <- c1
  C[2,2] <- c2
  
  Rn <- diag(p) - C/(N^gamma.x)
  x  <- matrix(0,N,p)
  
  for (j in 1:p)
  {
    for(t in 2:N) 
    {  
      x[t,j] <- Rn[j,j]*x[t-1,j] + innov.v[t,j]
    }
  }
  
  y <- matrix(0, nrow = N, ncol = 1)
  for (t in 2:N)
  {
    y[t,1] <- as.numeric(beta1)*x[t-1,1] + as.numeric(beta2)*x[t-1,2] + innov.u[t,1]
  }
  
  simulated.data <- structure( list( y = y, x  = x ) ) 
  return( simulated.data )
  
}# end of function

########################################################################################
### Function 2: Estimate Wald IVX statistic 
########################################################################################

sup_Wald_IVX_function <- function( Yt = Yt, Xt = Xt, Xlag = Xlag, delta = delta, cz = cz, pi0 = pi0 )
{# begin of function 
  
  # Insert data
  pi0   <- pi0
  Yt    <- as.matrix( Yt )
  Xt    <- as.matrix( Xt )
  Xlag  <- as.matrix( Xlag )
  delta <- delta 
  cz    <- cz
  
  ####################################################################
  ## IVX Estimation Step ## 
  ####################################################################
  
  # length size of the time series
  n <- NROW(Xlag)
  N <- (n+1)
  # number of predictors
  p  <- NCOL(Xlag)
  h  <- 1
  
  rn <- matrix(0, p, p)
  for (i in 1:p) 
  {
    rn[i, i] <- lm(Xt[, i] ~ 0 + Xlag[, i])$coefficients
  }
  
  # autoregressive residual estimation
  u.hat <- Xt - Xlag %*% rn
  
  # Rz assigns a common mildly integrated root to all regressors 
  Cz      <- diag(p)
  Cz[1,1] <- cz
  Cz[2,2] <- cz
  
  Rz <- ( diag(p) - Cz / ( N^delta ) ) 
  
  zt      <- matrix(0, n, p)
  zt[1, ] <- Xlag[1, ]
  
  diff <-  Xt -  Xlag
  
  for (t in 2:n) 
  {
    zt[t, ] <- Rz %*% zt[t - 1, ] + ( diff[t ,] )
  }
  
  # bandwith parameter
  m  <- floor( N^(1 / 3) ) 
  
  ## Use the above estimated Xt, Zt and Yt to construct the Wald IVX statistic
  sup_Wald_IVX <- 0
  
  # Step 1: Estimate OLS regression
  n <- NROW(Xlag)
  # number of predictors
  p <- NCOL(Xlag)
  
  # Step 2: Estimate the Wald IVX statistic
  k_seq <- as.matrix( floor(pi0 * N):(N - floor(pi0 * N) ) )
  dim   <- nrow( k_seq )    
  
  ones            <- matrix( 1, nrow = n, ncol = 1 )
  Wald_IVX_vector <- matrix( 0, nrow = dim, ncol = 1 )
  
  s <- 1
  
  for (s in 1:dim)
  {# begin of for-loop
    
    k <- k_seq[ s,1 ]
    
    ###################################################################
    # Step 1: Obtain the estimate of the variance of OLS regression 
    ###################################################################
    
    Xlag1             <- as.matrix( Xlag )
    Xlag1[(k+1): n, ] <- 0 
    
    Xlag2             <- as.matrix( Xlag )
    Xlag2[1:k, ]      <- 0 
    
    regressors <- cbind( Xlag1, Xlag2 )
    regressors <- as.matrix( regressors )
    
    model_OLS       <- lm( Yt ~ regressors - 1 )
    epsilon_hat     <- matrix( residuals( model_OLS ) )
    cov_epsilon_hat <- as.numeric( crossprod( epsilon_hat ) / n )
    
    ###################################################################
    # Step 2: Obtain the estimates of the IVX estimators
    ###################################################################
    
    # Estimate the Z1 and Z2 matrices corresponding to before and after the structural break
    Zt <- as.matrix( zt )
    
    Z1 <- as.matrix( Zt )
    for (i in (k+1): n)
    {
      Z1[i, ] <- 0 
    }
    
    Z2 <- as.matrix( Zt )
    for (i in 1: k)
    {
      Z2[i, ] <- 0 
    }
    
    # B1_hat_ivx    <- ( t( Yt ) %*% (Z1) ) %*% ( pracma::pinv( t(Xlag1) %*% Z1 ) )
    # B2_hat_ivx    <- ( t( Yt ) %*% (Z2) ) %*% ( pracma::pinv( t(Xlag2) %*% Z2 ) )
    
    B1_hat_ivx    <- ( crossprod (Yt, Z1) ) %*% ( pracma::pinv( crossprod (Xlag1, Z1) ) )
    B2_hat_ivx    <- ( crossprod (Yt, Z2) ) %*% ( pracma::pinv( crossprod (Xlag2, Z2) ) )
    Bivx_distance <- as.matrix( B2_hat_ivx  - B1_hat_ivx )
    
    ##########################################################################
    # Covariance matrix estimation
    
    Sigma.uu  <- matrix(0, p, p)
    
    for (i in 1:n)
    {
      Sigma.uu <- Sigma.uu + crossprod( u.hat[i, , drop = FALSE] )
    }
    
    Sigma.uu   <- Sigma.uu / n
    Sigma.ue   <- matrix(0, 1, p)
    
    for (i in 1:p) 
    {
      Sigma.ue[, i] <- sum( epsilon_hat * u.hat[, i] )
    }
    
    Sigma.ue  <- t( Sigma.ue ) / n
    # bandwith parameter
    m         <- floor( n^(1 / 3) ) 
    
    Lampda.uu <- matrix(0, p, p)
    
    for (i in 1:m) 
    {
      a <- matrix(0, p, p)
      for (j in (i + 1):n) 
      {
        a <- a + t( u.hat[j, , drop = F]) %*% u.hat[j - i, , drop = F]
      }
      Lampda.uu <- Lampda.uu + (1 - i / (m + 1)) * a
    }
    
    Lampda.uu <- Lampda.uu / n
    
    # Estimation of the Omegauu matrix 
    Omega.uu <- Sigma.uu + Lampda.uu + t(Lampda.uu)
    
    # Newey-West matrix estimators (Lampda.ue as defined in the paper)
    q <- matrix(0, m, p)
    
    for (i in 1:m) 
    {
      r <- matrix(0, n - i, p)
      
      for (j in (i + 1):n)
      {
        r[j - i, ] <- u.hat[j, , drop = F] * epsilon_hat[j - i] 
      }
      q[i, ] <- (1 - i / (1 + m)) * colSums(r)
    }
    
    # residue finds the column average of the q matrix
    residue  <- apply(q, 2, sum) / n
    Omega.ue <- Sigma.ue + as.matrix(residue) 
    sigma_FM <- cov_epsilon_hat - t(Omega.ue)%*% ( pracma::pinv(Omega.uu) ) %*% (Omega.ue)
    
    ##########################################################################
    
    # Estimate Q1 matrix 
    Z1.mean <- colMeans( as.matrix( Z1[1:k, ] ) )
    M1      <- crossprod( Z1 )* cov_epsilon_hat - kronecker( k * tcrossprod( Z1.mean ) , sigma_FM )
    Q_matrix_part1 <- ( pinv( t(Z1)%*%Xlag1 ) )%*%( M1 )%*%( pinv( t(Xlag1)%*%Z1 ) )
    
    # Estimate Q2 matrix 
    Z2.mean <- colMeans( as.matrix( Z2[(k+1):n, ] ) )
    M2      <- crossprod( Z2 )* cov_epsilon_hat - kronecker( (n - k) * tcrossprod( Z2.mean ) , sigma_FM )
    Q_matrix_part2 <- ( pinv( t(Z2)%*%Xlag2 ) )%*%( M2 )%*%( pinv( t(Xlag2)%*%Z2 ) )
    
    Q_matrix           <- ( Q_matrix_part1 + Q_matrix_part2 )
    Wald_IVX_statistic <- ( Bivx_distance )%*%( pracma::pinv( Q_matrix ) )%*%t( Bivx_distance ) 
    
    Wald_IVX_vector[s,1] <- Wald_IVX_statistic 
    
  }# end of for-loop
  
  # Obtain the sup-Wald statistic
  sup_Wald_IVX <- max( Wald_IVX_vector  )
  
  return( data.frame(sup_Wald_IVX)  )
  
}# end of function

########################################################################################
