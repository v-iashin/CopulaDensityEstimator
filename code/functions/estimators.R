online_copula_density_constant <- function(u, X, n0, alpha, nu) {
    
    h_constant <- function(n) {
        h <- 1.0 * matrix(1, ngrid, d) * n^(-1/(d+4))
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
    h_constant <- function(n) {
        h <- 1.0 * matrix(1, ngrid, d) * n^(-1/(d+4))
        return(h)
    }
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
    f <- matrix(1, ngrid, d)
    a <- matrix(NA, ngrid, d)
    sn <- matrix(NA, ngrid, n)
    cn <- matrix(NA, ngrid, n)
    mu_est <- matrix(0, d)
    sigma2_est <- matrix(1, d)
    
    ### Initialization using the first n0 observations
    X0 <- X[1:n0, ]
    h0 <- h_constant(n0)
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
    # TODO: Add a remark regarding special form of h for second order derivative
    h_constant <- function(n) {
        h <- 1.0 * matrix(1, ngrid, d) * n^(-1/(d+4))
        return(h)
    }
    h_marginals <- function(f, fdd, n, a) {
        eta <- 1 # int u^2K(u)du
        etat <- 1/(2*sqrt(pi)) # int K^2(u)du
        num <- f * (2 - a) * etat
        den <- 4 * fdd^2 * eta^2
        coeff <- (num/den)^(1/5)
        return(coeff * n^(-1/5))
    }
    h_marginals_dd <- function(f, fdd, n, a) {
        eta <- 1 # int u^2K(u)du
        etat <- 1/(2*sqrt(pi)) # int K^2(u)du
        num <- f * (2 - a) * etat
        den <- 4 * fdd^2 * eta^2
        coeff <- (num/den)^(1/9)
        return(coeff * n^(-1/9))
    }
    h_joint <- function(s, fll, n) {
        eta <- 1 # int u^2K(u)du
        etat <- 1/(2*sqrt(pi)) # int K^2(u)du
        d <- 2
        num <- d * etat^d * s * (d + 2)
        den <- 2 * eta^2 * rowSums(fll)^2 * (d + 4)
        coeff <- (num / den)^(1/(d+4))
        return(coeff * n^(-1/(d+4)))
    }
    h_joint_ll <- function(s, fll, n) {
        eta <- 1 # int u^2K(u)du
        etat <- 1/(2*sqrt(pi)) # int K^2(u)du
        d <- 2
        num <- d * etat^d * s * (d + 2)
        den <- 2 * eta^2 * rowSums(fll)^2 * (d + 4)
        coeff <- (num / den)^(1/(d+8))
        return(coeff * n^(-1/(d+8)))
    }
    
    n <- dim(X)[1]
    d <- dim(X)[2]
    ngrid <- dim(u)[1]
    
    Q <- matrix(0, ngrid, d)
    f <- matrix(1, ngrid, d)
    fdd <- matrix(1, ngrid, d) # 2nd order derivative marginals
    fll <- matrix(1, ngrid, d) # 2nd order derivatives joint
    a <- matrix(NA, ngrid, d)
    sn <- matrix(NA, ngrid, n)
    cn <- matrix(NA, ngrid, n)
    
    # Initialization using the first n0 observations
    X0 <- X[1:n0, ]
    h0 <- h_constant(n0)
    mu <- 0.1 * n0 * (-alpha)
    
    for (j in 1:d){
        Q[, j] <- quantile(X0[, j], u[, j])
        QmX <- replicate(n0, Q[, j]) - t(replicate(ngrid, X0[, j]))
        f[, j] <- rowMeans(Kh(QmX, h0[, j]))
        fdd[, j] <- rowMeans(Kddh(QmX, h0[, j]))
        a[, j] <- pmax(mu, pmin(f[, j], nu * (log(n0)+1)))
    }
    h0rep <- h0 %x% rep(1, n0)
    QmX <- Q %x% rep(1, n0) - apply(X0, 2, rep, ngrid)
    KQmX <- Kh(QmX, h0rep)
    
    QmX1 <- matrix(QmX[, 1], ngrid, n0, byrow = TRUE)
    QmX2 <- matrix(QmX[, 2], ngrid, n0, byrow = TRUE)
    fll[, 1] <- rowMeans(Kddh(QmX1, h0[, 1]) * Kh(QmX2, h0[, 2]))
    fll[, 2] <- rowMeans(Kddh(QmX2, h0[, 2]) * Kh(QmX1, h0[, 1]))
    s <- colMeans(matrix(apply(KQmX, 1, prod), n0, ngrid))
    sn[, n0] <- s
    cn[, n0] <- s / apply(a, 1, prod)

    # Online estimation
    for (i in (n0+1):n) {
        h_m <- h_marginals(f, fdd, i, a=4/5)
        h_m_dd <- h_marginals_dd(f, fdd, i, a=4/5)
        h_j <- h_joint(s, fll, i)
        h_j_ll <- h_joint_ll(s, fll, i)
        
        Xrep <- t(replicate(ngrid, X[i, ]))
        QmX <- Q - Xrep
        W <- apply(Kh(QmX, h_j), 1, prod)
        f <- (1 - 1/i) * f + (1/i) * Kh(QmX, h_m)
        fdd <- (1 - 1/i) * fdd + (1/i) * Kddh(QmX, h_m_dd)
        KQmX <- Kh(QmX, h_j_ll)
        KddQmX <- Kddh(QmX, h_j_ll)
        fll[, 1] <- (1 - 1/i) * fll[, 1] + (1/i) * KddQmX[, 1] * KQmX[, 2]
        fll[, 2] <- (1 - 1/i) * fll[, 2] + (1/i) * KddQmX[, 2] * KQmX[, 1]
        Q <- Q + (1 / (a*i)) * (u - as.integer(Xrep <= Q))
        a <- pmax(pmin(f, nu*(log(i)+1)), mu)
        s <- (1 - 1/i) * sn[, i-1] + (1/i) * W
        sn[, i] <- s
        cn[, i] <- s / apply(a, 1, prod)
    }
    return(cn)
}