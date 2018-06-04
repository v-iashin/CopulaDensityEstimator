library(latex2exp)

plot_results <- function(mse1, mse2, mse3, title, theta, pos, save_plot, filename) {
    # A wrapper for the plotting function tha decided whether to save a plot and 
    # at what path. Plots the results (mse) for each of the three bandwidths. 
    # Args:
    #   mse1, mse2, mse3: Arrays (nmc, n). Results of simulations.
    #   title: TeX. A formatted title for a plot.
    #   theta: float. A parameter of the copula that is used in the experiment.
    #   pos: String (for ex: 'topright' or 'bottomleft'): Position of the legend.
    #   save_plot: Boolean. Whether to save plot.
    #   filename: String. If save_plot is TRUE save plot to the filename.
    # Returns: nothing
     
    plot_results_ <- function(mse1, mse2, mse3, title, theta, pos) {
        # Plots the results (mse) for each of the three bandwidths.
        # Args:
        #   See plot_results
        # Returns: nothing
        
        # Dimensions and limits for the plot
        n0 <- sum(is.na(colMeans(mse1)))
        nmc <- dim(mse1)[1]
        nmax <- dim(mse1)[2]
        max_mse <- quantile(c(mse1, mse2, mse3), 0.97, na.rm = T)
        min_mse <- 0
        
        # COLORS
        # Red, Green, & Blue 
        # ColMain1 <- rgb(0.3, 0.6, 0.7)
        # ColCE1 <- rgb(0.3, 0.6, 0.7, alpha = 0.3)
        # ColMain2 <- rgb(0.7, 0.3, 0.3)
        # ColCE2 <- rgb(0.7, 0.3, 0.3, alpha = 0.3)
        # ColMain3 <- rgb(0.0, 0.8, 0.4)
        # ColCE3 <- rgb(0.0, 0.8, 0.4, alpha = 0.3)
        # Black & White
        ColMain1 <- rgb(0.1, 0.1, 0.1)
        ColCE1 <- rgb(0.1, 0.1, 0.1, alpha = 0.3)
        ColMain2 <- ColMain1
        ColCE2 <- ColCE1
        ColMain3 <- ColMain1
        ColCE3 <- ColCE1
        
        # CONFIDENCE INTERVALS (95 %)
        means1 <- colMeans(mse1)
        sd1 <- apply(mse1, 2, sd)
        confintUpper1 <- means1 + qnorm(0.975) * sd1 / sqrt(nmc)
        confintLower1 <- means1 - qnorm(0.975) * sd1 / sqrt(nmc)
        
        means2 <- colMeans(mse2)
        sd2 <- apply(mse2, 2, sd)
        confintUpper2 <- means2 + qnorm(0.975) * sd2 / sqrt(nmc)
        confintLower2 <- means2 - qnorm(0.975) * sd2 / sqrt(nmc)
        
        means3 <- colMeans(mse3)
        sd3 <- apply(mse3, 2, sd)
        confintUpper3 <- means3 + qnorm(0.975) * sd3 / sqrt(nmc)
        confintLower3 <- means3 - qnorm(0.975) * sd3 / sqrt(nmc)
        
        # PLOTS
        # 1
        plot(means1, type = 'l', col = ColMain1, lty = 'solid',
             main = TeX(sprintf("%s, $\\theta = %.1f$", title, theta)), 
             ylim = c(min_mse, max_mse), xlim = c(n0, floor(nmax*1.00)), 
             xlab = 'n', ylab = 'Mean Squared Error', lwd = 2)
        lines(means1, lwd = 2, col = ColMain1)
        # confidence interval
        polygon(c(1:nmax, rev(1:nmax)), c(confintLower1, rev(confintUpper1)), 
                col = ColCE1, border = FALSE)
        # 2
        lines(means2, lwd = 2, col = ColMain2, lty = 'dotted')
        # confidence interval
        polygon(c(1:nmax, rev(1:nmax)), c(confintLower2, rev(confintUpper2)), 
                col = ColCE2, border = FALSE)
        # 3
        lines(means3, lwd = 2, col = ColMain3, lty = 'longdash')
        # confidence interval
        polygon(c(1:nmax, rev(1:nmax)), c(confintLower3, rev(confintUpper3)), 
                col = ColCE3, border = FALSE)
        
        # LEGEND
        legend_box <- c(sprintf('Constant (%.4f)', means1[nmax]), 
                        sprintf('Silverman (%.4f)', means2[nmax]), 
                        sprintf('Specialized (%.4f)', means3[nmax]))
        legend(pos, legend = legend_box, title = "Bandwidth (lastest MSE)",
               col = c(ColMain1, ColMain2, ColMain3), lty = c(1, 3, 5), lwd = 2, 
               cex = 0.8)
        
        # GRID
        grid()
    }
    
    if (save_plot) {
        pdf(filename, width = 5, height = 4)
        plot_results_(mse1, mse2, mse3, title, theta, pos)
        dev.off()
    }
    
    else {
        plot_results_(mse1, mse2, mse3, title, theta, pos)
    }
}