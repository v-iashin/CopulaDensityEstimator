Kh <- function(x, h){
    # Calculates Gaussian kernel for given x and bandwidths h.
    # Args:
    #   x: Array. Values
    #   h: Array (of the same shape as x). Bandwidths
    # Returns:
    #   Array (of the same shape as x). Values of the Gaussian kernel.
    return(dnorm(x/h)/h)
}

Kddh <- function(x, h){
    # Calculates second order derivative for Gaussian kernel 
    # for given x and bandwidths h.
    # Args:
    #   x: Array. Values
    #   h: Array (of the same shape as x). Bandwidths
    # Returns:
    #   Array (of the same shape as x). Second order derivative for Gaussian kernel.
    return(((x/h)^2 - 1)*dnorm(x/h)/h^3)
}