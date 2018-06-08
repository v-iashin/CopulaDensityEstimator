# Don't forget to set working directory

library(copula)

mai = c(0.1, 0.1, 0.1, 0.1)
width = 8
height = 2.5

pdf('frank.pdf', width = width, height = height)
old.par <- par(mfrow=c(1, 3), mai = mai)
theta <- -1
cop <- frankCopula(theta)
persp(cop, dCopula, nticks=2, zlab='', main=bquote(theta == .(theta)))
theta <- 0.001
cop <- frankCopula(theta)
persp(cop, dCopula, nticks=2, zlab='', main=bquote(theta == .(theta)))
theta <- 1
cop <- frankCopula(theta)
persp(cop, dCopula, nticks=2, zlab='', main=bquote(theta == .(theta)))
par(old.par)
dev.off()

pdf('gumbel.pdf', width = width, height = height)
old.par <- par(mfrow=c(1, 3), mai = mai)
theta <- 2
cop <- gumbelCopula(theta)
persp(cop, dCopula, nticks=2, zlab='', main=bquote(theta == .(theta)))
theta <- 5
cop <- gumbelCopula(theta)
persp(cop, dCopula, nticks=2, zlab='', main=bquote(theta == .(theta)))
theta <- 10
cop <- gumbelCopula(theta)
persp(cop, dCopula, nticks=2, zlab='', main=bquote(theta == .(theta)))
par(old.par)
dev.off()

pdf('clayton.pdf', width = width, height = height)
old.par <- par(mfrow=c(1, 3), mai = mai)
theta <- -0.5
cop <- claytonCopula(theta)
persp(cop, dCopula, nticks=2, zlab='', main=bquote(theta == .(theta)))
theta <- 0.5
cop <- claytonCopula(theta)
persp(cop, dCopula, nticks=2, zlab='', main=bquote(theta == .(theta)))
theta <- 10
cop <- claytonCopula(theta)
persp(cop, dCopula, nticks=2, zlab='', main=bquote(theta == .(theta)))
par(old.par)
dev.off()

pdf('student.pdf', width = width, height = height)
old.par <- par(mfrow=c(1, 3), mai = mai)
theta <- -0.5
cop <- tCopula(theta)
persp(cop, dCopula, nticks=2, zlab='', main=bquote(theta == .(theta)))
theta <- 0
cop <- tCopula(theta)
persp(cop, dCopula, nticks=2, zlab='', main=bquote(theta == .(theta)))
theta <- 0.5
cop <- tCopula(theta)
persp(cop, dCopula, nticks=2, zlab='', main=bquote(theta == .(theta)))
par(old.par)
dev.off()

pdf('normal.pdf', width = width, height = height)
old.par <- par(mfrow=c(1, 3), mai = mai)
theta <- -0.5
cop <- normalCopula(theta)
persp(cop, dCopula, nticks=2, zlab='', main=bquote(theta == .(theta)))
theta <- 0
cop <- normalCopula(theta)
persp(cop, dCopula, nticks=2, zlab='', main=bquote(theta == .(theta)))
theta <- 0.5
cop <- normalCopula(theta)
persp(cop, dCopula, nticks=2, zlab='', main=bquote(theta == .(theta)))
par(old.par)
dev.off()