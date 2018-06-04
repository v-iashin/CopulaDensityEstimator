simulation <- function(n, nmc, cop, n0, alpha, nu, seed) {
    # Runs Monte-Carlo simulations for three different copula density 
    # estimators. Returns a list with calculated log of mean squared errors 
    # between theoretical and estimated values for each of the three bandwidths.
    # Args:
    #   n: Int. Number of generated observations.
    #   nmc: Int. Number of Monte-Carlo 
    #   cop: Copula. Copula class.
    #   n0: Int. Number of observations for the initialization.
    #   alpha: Float. Parameter of `mu`.
    #   nu: Positive float (or Int): Parameter of `a`.
    # Returns:
    #   List. A log of mean squared errors for each of the three bandwidths.
    
    # placeholders
    mse1 <- matrix(0, nmc, n)
    mse2 <- matrix(0, nmc, n)
    mse3 <- matrix(0, nmc, n)
    
    # set seed
    set.seed(seed)
    
    # start the progress bar
    pb <- txtProgressBar(min = 1, max = nmc, style = 3)
    
    for (i in 1:nmc) {
        # generate data
        U <- rCopula(n, cop)
        X <- qnorm(U)
        
        # default (constant)
        estim1 <- online_copula_density_constant(u, X, n0, alpha, nu)
        err1 <- estim1 - cbind(dCopula(u, cop)) %*% rep(1, n)
        mse1[i, ] <- colMeans(err1^2)
        
        # Silverman's
        estim2 <- online_copula_density_silverman(u, X, n0, alpha, nu)
        err2 <- estim2 - cbind(dCopula(u, cop)) %*% rep(1, n)
        mse2[i, ] <- colMeans(err2^2)
        
        # Specialized
        estim3 <- online_copula_density_specialized(u, X, n0, alpha, nu)
        err3 <- estim3 - cbind(dCopula(u, cop)) %*% rep(1, n)
        mse3[i, ] <- colMeans(err3^2)
        
        # update progress bar
        setTxtProgressBar(pb, i)
    }
    return(list(mse1, mse2, mse3))
}