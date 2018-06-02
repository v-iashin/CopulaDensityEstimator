Kh <- function(x, h){
    return(dnorm(x/h)/h)
}

Kddh <- function(x, h){
    return(((x/h)^2 - 1)*dnorm(x/h)/h^3)
}