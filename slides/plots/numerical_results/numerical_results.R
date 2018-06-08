# Don't forget to set working directory

library(copula)
library(latex2exp)

width = 4
height = 4

pdf('frank0.pdf', width = width, height = height)
theta <- 0
cop <- frankCopula(theta)
title <- TeX('Independent Copula')
persp(cop, dCopula, nticks = 2, zlab = '', main=title)
dev.off()

pdf('student05.pdf', width = width, height = height)
theta <- 0.5
cop <- tCopula(theta)
title <- TeX('Student Copula $\\theta = 0.5$')
persp(cop, dCopula, nticks = 2, zlab = '', main=title)
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