library(copula)

source('functions/estimators.R', chdir = TRUE)
source('functions/kernels.R', chdir = TRUE)
source('functions/simulation.R', chdir = TRUE)
source('functions/plotting.R', chdir = TRUE)

# create grid, 9 by 9 in this case
u <- as.matrix(expand.grid(1:9 / 10, 1:9 / 10))

# number of simulations
nmc <- 100

# set parameters
n <- 5000
n0 <- 10
d <- 2 # dim
alpha <- 0.25 # mu
nu <- 1 # a

# change to TRUE if you wish to save plots to .pdf files
save_plot <- FALSE

# FRANK COPULA
{
    theta <- -1
    cop <- frankCopula(theta)
    sim_results <- simulation(n, nmc, cop, n0, alpha, nu, seed = 13)
    plot_results(sim_results[[1]], sim_results[[2]], sim_results[[3]], 
                 title = 'Frank Copula', theta, pos = 'topright', 
                 save_plot, '../text/plots/experiment_results/frank_1.pdf')
    
    theta <- 0
    cop <- frankCopula(theta)
    sim_results <- simulation(n, nmc, cop, n0, alpha, nu, seed = 13)
    plot_results(sim_results[[1]], sim_results[[2]], sim_results[[3]], 
                 title = '(Independence) Frank Copula', theta, pos = 'topright', 
                 save_plot, '../text/plots/experiment_results/frank0.pdf')
    
    theta <- 1
    cop <- frankCopula(theta)
    sim_results <- simulation(n, nmc, cop, n0, alpha, nu, seed = 13)
    plot_results(sim_results[[1]], sim_results[[2]], sim_results[[3]], 
                 title = 'Frank Copula', theta, pos = 'topright', 
                 save_plot, '../text/plots/experiment_results/frank1.pdf')
}

# CLAYTON COPULA
{
    theta <- -0.5
    cop <- claytonCopula(theta)
    sim_results <- simulation(n, nmc, cop, n0, alpha, nu, seed = 13)
    plot_results(sim_results[[1]], sim_results[[2]], sim_results[[3]], 
                 title = 'Clayton Copula', theta, pos = 'bottomleft', 
                 save_plot, '../text/plots/experiment_results/clayton_05.pdf')
    
    theta <- 0.5
    cop <- claytonCopula(theta)
    sim_results <- simulation(n, nmc, cop, n0, alpha, nu, seed = 13)
    plot_results(sim_results[[1]], sim_results[[2]], sim_results[[3]], 
                 title = 'Clayton Copula', theta, pos = 'topright', 
                 save_plot, '../text/plots/experiment_results/clayton05.pdf')
    
    theta <- 10
    cop <- claytonCopula(theta)
    sim_results <- simulation(n, nmc, cop, n0, alpha, nu, seed = 13)
    plot_results(sim_results[[1]], sim_results[[2]], sim_results[[3]], 
                 title = 'Clayton Copula', theta, pos = 'bottomleft', 
                 save_plot, '../text/plots/experiment_results/clayton10.pdf')
}

# STUDENT COPULA
{
    theta <- -0.5
    cop <- tCopula(theta)
    sim_results <- simulation(n, nmc, cop, n0, alpha, nu, seed = 13)
    plot_results(sim_results[[1]], sim_results[[2]], sim_results[[3]], 
                 title = 'Student Copula', theta, pos = 'topright', 
                 save_plot, '../text/plots/experiment_results/student_05.pdf')
    
    theta <- 0
    cop <- tCopula(theta)
    sim_results <- simulation(n, nmc, cop, n0, alpha, nu, seed = 13)
    plot_results(sim_results[[1]], sim_results[[2]], sim_results[[3]], 
                 title = 'Student Copula', theta, pos = 'topright', 
                 save_plot, '../text/plots/experiment_results/student0.pdf')
    
    theta <- 0.5
    cop <- tCopula(theta)
    sim_results <- simulation(n, nmc, cop, n0, alpha, nu, seed = 13)
    plot_results(sim_results[[1]], sim_results[[2]], sim_results[[3]], 
                 title = 'Student Copula', theta, pos = 'topright', 
                 save_plot, '../text/plots/experiment_results/student05.pdf')
}
