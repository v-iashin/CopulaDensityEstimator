# Do not forget to set your working directory

library(ggplot2)
library(ggExtra)
library(gridExtra)
library(grid)
library(latex2exp) # TeX("$P_a$")

plot_independent_sample <- function(size, save_plot=FALSE, filename) {
    set.seed(13)
    a <- rnorm(size)
    w <- rlnorm(size)
    X <- data.frame(a, w)
    Q1 <- pnorm(a)
    Q2 <- plnorm(w)
    Q <- data.frame(Q1, Q2)
    
    sc1 <- ggplot(X, aes(a, w))
    sc1 <- sc1 + geom_point(size = 0.3, color = 'firebrick') 
    sc1 <- sc1 + theme_bw(9)
    sc1 <- sc1 + ggtitle('Joint PDF of random variables')
    sc1 <- sc1 + theme(plot.title = element_text(hjust = 0.5))
    sc1 <- ggMarginal(sc1, colour = 'firebrick', 
                      fill = 'firebrick', alpha = 0.3)
    
    sc2 <- ggplot(Q, aes(Q1, Q2))
    sc2 <- sc2 + geom_point(size = 0.3, color = 'lightseagreen') 
    sc2 <- sc2 + theme_bw(9)
    sc2 <- sc2 + ggtitle('Joint PDF of marginal CDFs')
    sc2 <- sc2 + theme(plot.title = element_text(hjust = 0.5))
    sc2 <- sc2 + xlab(TeX('$F_a(a)$')) + ylab(TeX('$F_w(w)$'))
    sc2 <- ggMarginal(sc2, colour = 'lightseagreen', 
                      fill = 'lightseagreen', alpha = 0.3)
    
    if (save_plot) {
        pdf(filename, width = 8, height = 4)
        grid.arrange(sc1, sc2, nrow = 1, top = ' ')
        dev.off()
    }
    else {
        grid.arrange(sc1, sc2, nrow = 1, top = ' ')
    }
}

plot_independent_sample(size = 10000, filename = 'indep.pdf', save_plot = TRUE)

plot_dependent_sample <- function(size, save_plot=FALSE, filename) {
    set.seed(13)
    rho <- 0.8
    a <- rnorm(size)
    w <- exp(rho * a + sqrt(1-rho^2) * rnorm(size))
    X <- data.frame(a, w)
    Q1 <- pnorm(a)
    Q2 <- plnorm(w)
    Q <- data.frame(Q1, Q2)
    
    sc1 <- ggplot(X, aes(a, w))
    sc1 <- sc1 + geom_point(size = 0.3, color = 'firebrick') 
    sc1 <- sc1 + theme_bw(9)
    sc1 <- sc1 + ggtitle('Joint PDF of random variables')
    sc1 <- sc1 + theme(plot.title = element_text(hjust = 0.5))
    sc1 <- ggMarginal(sc1, colour = 'firebrick', 
                      fill = 'firebrick', alpha = 0.3)
    
    sc2 <- ggplot(Q, aes(Q1, Q2))
    sc2 <- sc2 + geom_point(size = 0.3, color = 'lightseagreen') 
    sc2 <- sc2 + theme_bw(9)
    sc2 <- sc2 + ggtitle('Joint PDF of marginal CDFs')
    sc2 <- sc2 + theme(plot.title = element_text(hjust = 0.5))
    sc2 <- sc2 + xlab(TeX('$F_a(a)$')) + ylab(TeX('$F_w(w)$'))
    sc2 <- ggMarginal(sc2, colour = 'lightseagreen', 
                      fill = 'lightseagreen', alpha = 0.3)
    
    if (save_plot) {
        pdf(filename, width = 8, height = 4)
        grid.arrange(sc1, sc2, nrow = 1, top = ' ')
        dev.off()
    }
    else {
        grid.arrange(sc1, sc2, nrow = 1, top = ' ')
    }
}

plot_dependent_sample(size = 10000, filename = 'dep.pdf', save_plot = TRUE)

plot_copulas <- function(size, save_plot=FALSE, filename) {
    set.seed(13)
    a <- rnorm(size)
    w <- rlnorm(size)
    X <- data.frame(a, w)
    Q1 <- pnorm(a)
    Q2 <- plnorm(w)
    Q <- data.frame(Q1, Q2)
    
    sc1 <- ggplot(Q, aes(Q1, Q2))
    sc1 <- sc1 + geom_point(size = 0.3, color = 'lightseagreen') 
    sc1 <- sc1 + theme_bw(9)
    sc1 <- sc1 + ggtitle('Copula density (independent sample)')
    sc1 <- sc1 + theme(plot.title = element_text(hjust = 0.5))
    sc1 <- sc1 + xlab(TeX('$F_a(a)$')) + ylab(TeX('$F_w(w)$'))
    sc1 <- ggMarginal(sc1, colour = 'lightseagreen', 
                      fill = 'lightseagreen', alpha = 0.3)
    
    set.seed(13)
    rho <- 0.8
    a <- rnorm(size)
    w <- exp(rho * a + sqrt(1-rho^2) * rnorm(size))
    X <- data.frame(a, w)
    Q1 <- pnorm(a)
    Q2 <- plnorm(w)
    Q <- data.frame(Q1, Q2)
    
    sc2 <- ggplot(Q, aes(Q1, Q2))
    sc2 <- sc2 + geom_point(size = 0.3, color = 'lightseagreen') 
    sc2 <- sc2 + theme_bw(9)
    sc2 <- sc2 + ggtitle('Copula density (dependent sample)')
    sc2 <- sc2 + theme(plot.title = element_text(hjust = 0.5))
    sc2 <- sc2 + xlab(TeX('$F_a(a)$')) + ylab(TeX('$F_w(w)$'))
    sc2 <- ggMarginal(sc2, colour = 'lightseagreen', 
                      fill = 'lightseagreen', alpha = 0.3)
    
    final_plot <- grid.arrange(sc1, sc2, nrow = 1, top = ' ')
    
    if (save_plot) {
        ggsave(filename = filename, device = 'pdf', width = 7.5, height = 3.6,
               plot = final_plot)
    }
    else {
        print(final_plot)
    }
}

plot_copulas(size = 10000, filename = 'copulas.pdf', save_plot = TRUE)

blank_plot <- function(filename) {
    final_plot <- ggplot() + theme_minimal()
    
    ggsave(filename = 'blank.png', plot = final_plot, 
           device = 'png', width = 7.5, height = 3.6)
}

blank_plot('blank.png')
