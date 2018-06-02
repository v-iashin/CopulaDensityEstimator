online_copula_density_constant <- function(u, X, n0, alpha, nu) {
    
    h_constant <- function(n) {
        h <- 1.0 * matrix(1, ngrid, d) * n^(-1/10)
        return(h)
    }
    
    n <- dim(X)[1]
    d <- dim(X)[2]
    ngrid <- dim(u)[1]
    
    Q <- matrix(0, ngrid, d)
    QmX <- matrix(0, ngrid, d)
    f <- matrix(1, ngrid, d)
    a <- matrix(NA, ngrid, d)
    sn <- matrix(NA, ngrid, n)
    cn <- matrix(NA, ngrid, n)
    
    ### Initialization using the first n0 observations
    X0 <- X[1:n0, ]
    h0 <- h_constant(n0)
    mu0 <- 0.1 * n0 ^ (-alpha)
    
    for (j in 1:d) {
        Q[, j] <- quantile(X0[, j], u[, j])
        QmX <- replicate(n0, Q[, j]) - t(replicate(ngrid, X0[, j]))
        f[, j] <- rowMeans(Kh(QmX, h0[, j]))
        a[, j] <- pmax(mu0, pmin(f[, j], nu * (log(n0)+1)))
    }
    QmX <- Q %x% rep(1, n0) - apply(X0, 2, rep, ngrid)
    h0rep <- h0 %x% rep(1, n0)
    KQmX <- Kh(QmX, h0rep)
    
    sn[, n0] <- colMeans(matrix(apply(KQmX, 1, prod), n0, ngrid))
    cn[, n0] <- sn[, n0] / apply(a, 1, prod)
    
    ### Online estimation
    for (i in (n0+1):n) {
        h <- h_constant(i)
        mu <- 0.1 * i ^ (-alpha)
        
        Xrep <- t(replicate(ngrid, X[i, ]))
        QmX <- Q - Xrep
        KQmX <- apply(Kh(QmX, h), 1, prod)
        Q <- Q + (1 / (i*a)) * (u - as.integer(Xrep <= Q))
        f <- (1 - 1/i)*f + (1/i)*Kh(QmX, h)
        a <- pmax(pmin(f, nu*(log(i)+1)), mu)
        sn[, i] <- (1 - 1/i)*sn[, i-1] + (1/i)*KQmX 
        cn[, i] <- sn[, i] / apply(a, 1, prod)
    }
    return(cn)
}

online_copula_density_silverman <- function(u, X, n0, alpha, nu) {
    
    h_silverman <- function(X, n, ngrid, mu_est, sigma2_est) {
        # update mu
        mu_est_prev <- mu_est
        mu_est_next <- mu_est_prev + (X - mu_est_prev) / n
        # update variance
        sigma2_est <- (n-2)/(n-1) * sigma2_est + (X - mu_est_prev)^2 / n
        # calculate h
        h <- (4/3)^(1/5) * sqrt(sigma2_est) * n^(-1/5)
        h_rep <- t(h) %x% rep(1, ngrid)
        return(list(h_rep, mu_est, sigma2_est))
    }
    
    n <- dim(X)[1]
    d <- dim(X)[2]
    ngrid <- dim(u)[1]

    Q <- matrix(0, ngrid, d)
    QmX <- matrix(0, ngrid, d)
    f <- matrix(1, ngrid, d)
    a <- matrix(NA, ngrid, d)
    sn <- matrix(NA, ngrid, n)
    cn <- matrix(NA, ngrid, n)
    mu_est <- matrix(0, d)
    sigma2_est <- matrix(1, d)
    
    ### Initialization using the first n0 observations
    X0 <- X[1:n0, ]
    h0 <- 1.0 * matrix(1, ngrid, d) * n0^(-1/10)
    mu0 <- 0.1 * n0 ^ (-alpha)
    
    for (j in 1:d) {
        Q[, j] <- quantile(X0[, j], u[, j])
        QmX <- replicate(n0, Q[, j]) - t(replicate(ngrid, X0[, j]))
        mu_est[j] <- mean(QmX)
        sigma2_est[j] <- sd(QmX)^2
        f[, j] <- rowMeans(Kh(QmX, h0[, j]))
        a[, j] <- pmax(mu0, pmin(f[, j], nu * (log(n0)+1)))
    }
    QmX <- Q %x% rep(1, n0) - apply(X0, 2, rep, ngrid)
    h0rep <- h0 %x% rep(1, n0)
    KQmX <- Kh(QmX, h0rep)
    
    sn[, n0] <- colMeans(matrix(apply(KQmX, 1, prod), n0, ngrid))
    cn[, n0] <- sn[, n0] / apply(a, 1, prod)
    
    ### Online estimation
    for (i in (n0+1):n) {
        mu <- 0.1 * i ^ (-alpha)
        h_res <- h_silverman(colMeans(QmX), i, ngrid, mu_est, sigma2_est)
        h <- h_res[[1]]
        mu_est <- h_res[[2]]
        sigma2_est <- h_res[[3]]
        
        Xrep <- t(replicate(ngrid, X[i, ]))
        QmX <- Q - Xrep
        KQmX <- apply(Kh(QmX, h), 1, prod)
        Q <- Q + (1 / (i*a))*(u - as.integer(Xrep <= Q))
        f <- (1 - 1/i)*f + (1/i)*Kh(QmX, h)
        a <- pmax(pmin(f, nu*(log(i)+1)), mu)
        sn[, i] <- (1 - 1/i)*sn[, i-1] + (1/i)*KQmX
        cn[, i] <- sn[, i] / apply(a, 1, prod)
    }
    return(cn)
}

online_copula_density_specialized <- function(u, X, n0, alpha, nu) {
    n = dim(X)[1]; d = dim(X)[2]; ngrid = dim(u)[1]
    eta = 1 # int u^2K(u)du
    etat = 1/(2*sqrt(pi)) # int K^2(u)du
    
    # Initialization using the first n0 observations
    X0 = X[1:n0, ]
    x1 = colSums(X0); x2 = colSums(X0^2); 
    #hinit = 1.06*sqrt(x2/n0-(x1/n0)^2)*n0^(-1/10)
    hinit = c(1,1)*n0^(-1/(d+4))
    Q =  matrix(rep(0, d*ngrid), nrow = ngrid)
    f =  matrix(rep(1, d*ngrid), nrow = ngrid)
    fdd =  matrix(rep(1, d*ngrid), nrow = ngrid) # 2nd order derivative marginals
    fll =  matrix(rep(1, d*ngrid), nrow = ngrid) # 2nd order derivatives joint
    a =  matrix(rep(1, d*ngrid), nrow = ngrid) 
    for (id in 1:d){
        Q[, id] = quantile(X0[, id], u[, id]) 
        QmX = matrix(rep(Q[, id], n0), nrow = ngrid, byrow=FALSE) - matrix(rep(X0[, id], ngrid), nrow = ngrid, byrow=TRUE)
        f[, id] = rowMeans(Kh(QmX, hinit[id])) 
        fdd[, id] = rowMeans(Kddh(QmX, hinit[id]))
        a[, id] = pmax(0.1*n0^(-alpha), pmin(f[, id], nu*(log(n0)+1)))
    }
    QmX1 = matrix(rep(Q[, 1], n0), nrow = ngrid, byrow=FALSE) - matrix(rep(X0[, 2], ngrid), nrow = ngrid, byrow=TRUE)
    QmX2 = matrix(rep(Q[, 2], n0), nrow = ngrid, byrow=FALSE) - matrix(rep(X0[, 2], ngrid), nrow = ngrid, byrow=TRUE)
    fll[, 1] = rowMeans(Kddh(QmX1, hinit[1])*Kh(QmX2, hinit[2]))
    fll[, 2] = rowMeans(Kddh(QmX2, hinit[2])*Kh(QmX1, hinit[1]))
    sn = matrix(rep(NA, n*ngrid), nrow = ngrid)
    cn = matrix(rep(NA, n*ngrid), nrow = ngrid)   
    KQmX = Kh(matrix(as.vector(apply(Q,1,rep,n0)), ncol=d, byrow=TRUE) - apply(X0, 2, rep, ngrid), matrix(rep(hinit, n0*ngrid), nrow=n0*ngrid, byrow=TRUE))
    sn[, n0] = as.vector(tapply(apply(KQmX, 1, prod), gl(ngrid, n0), mean))
    cn[, n0] = sn[, n0]/apply(a, 1, prod)
    
    # Online estimation
    for (i in (n0+1):n) {
        x1 = x1 + X[i,]; x2 = x2 + X[i,]^2;
        #hi = .5*c(1,1)*i^(-1/10) 
        #hi = 0.609*colMeans((f/fdd^2)^(1/5))*i^(-1/5)     # h minimizing the EQI
        fu = sn[,i-1]
        fll2 = rowSums(fll)^2
        himarg = ((3*etat*f)/(10*(eta*fdd)^2))^(1/5)*i^(-1/5)
        #himarg = matrix(rep(.5*c(1,1)*i^(-1/5), ngrid), nrow=ngrid)
        himargder = ((3*etat*f)/(10*(eta*fdd)^2))^(1/9)*i^(-1/9)
        #himargder = matrix(rep(.5*c(1,1)*i^(-1/5), ngrid), nrow=ngrid)
        hijoin = ((d*(d+2)*etat^d*fu)/(2*(d+4)*eta^2*fll2))^(1/(d+4))*i^(-1/(d+4))
        #hijoin = rep(.5*i^(-1/10), ngrid)
        hijoinder = ((d*(d+2)*etat^d*fu)/(2*(d+4)*eta^2*fll2))^(1/(d+8))*i^(-1/(d+8))
        Xrepi = matrix(rep(X[i,], each=ngrid), nrow=ngrid)     
        sn[ ,i] = (1 - 1/i)*sn[ ,i-1] + (1/i)*apply(Kh(Q-Xrepi, hijoin), 1, prod)   
        f = (1-1/i)*f + (1/i)*Kh(Q-Xrepi, himarg)
        fdd = (1-1/i)*fdd + (1/i)*Kddh(Q-Xrepi, himargder)
        fll[, 1] = (1-1/i)*fll[,1] + (1/i)*Kddh(Q[,1]-Xrepi[,1], hijoinder)*Kh(Q[,2]-Xrepi[,2], hijoinder)
        fll[, 2] = (1 - 1/i)*fll[,2] + (1/i)*Kddh(Q[,2]-Xrepi[,2], hijoinder)*Kh(Q[,1]-Xrepi[,1], hijoinder)
        Q = Q + (u - (Xrepi <= Q)*1)/(a*i) 
        a = pmax(pmin(f, nu*(log(i)+1)), 0.1*i^(-alpha))
        cn[ ,i] = sn[ ,i]/apply(a, 1, prod)
    }
    return(list(cn, f, fdd, sn, fll, himarg, hijoin))
}