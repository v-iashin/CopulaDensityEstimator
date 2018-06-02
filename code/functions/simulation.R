simulation <- function(n, nmc, cop, U, n0, alpha, nu, seed) {
    mse1 <- matrix(0, nmc, n)
    mse2 <- matrix(0, nmc, n)
    mse3 <- matrix(0, nmc, n)
    
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
        err3 <- estim3[[1]] - cbind(dCopula(u, cop)) %*% rep(1, n)
        mse3[i, ] <- colMeans(err3^2)
        
        # update progress bar
        setTxtProgressBar(pb, i)
    }
    return(list(mse1, mse2, mse3))
}