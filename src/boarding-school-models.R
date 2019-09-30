#!/usr/bin/R
library(xtable)
library(deSolve)
library(data.table)
library(bbmle)

data <- read.csv("data/boarding-school.csv")

closed.sir <- function (t, states, params) {
    s0 <- states[1]
    i0 <- states[2]
    r0 <- states[3]

    beta <- params[1]
    gamma <- params[2]

    dS <- -beta * s0 * i0
    dI <- beta * s0 * i0 - gamma * i0
    dR <- gamma * i0
    list(c(dS, dI, dR))
}

mle.sir <- function(b, g) {
    t <- seq(0, 14)    
    beta <- exp(b)
    gamma <- exp(g)

    results <- as.data.frame(ode(y=c(S=762, I=1, R=0), 
                             times=t, 
                             closed.sir,
                             parms=c(beta, gamma)))
    nll <- -sum(dpois(x=data$cases, lambda=tail(results$I, 14),  log=TRUE))
    return(nll)
    #return(results)
}

initial <- list(b=-6, g=-0.7)
fit0 <- mle2(mle.sir, start=initial)

fit <- mle2(mle.sir, start=as.list(coef(fit0)))

pred <- as.data.frame(ode(c(S=762, I=1, R=0), times=seq(0, 14, 0.5), closed.sir, parms=c(exp(coef(fit)))))

plot(cases~day, data=data, type="b", ylab="Infected")
lines(pred$I~pred$time)

sir.beta <- as.numeric(exp(coef(fit)[1])*763)
sir.gamma <- as.numeric(exp(coef(fit)[2]))
sir.nll <- as.numeric(logLik(fit))
sir.aic <- AIC(fit)

closed.seir <- function(t, states, params) {
    s0 <- states[1]
    e0 <- states[2]
    i0 <- states[3]
    r0 <- states[4]

    beta <- params[1]
    gamma <- params[2]
    mu <- 0
    sigma <- params[3]

    dS <- 0 - (beta * i0 + 0) * s0
    dE <- beta * s0 * i0 - (0 + sigma) * e0
    dI <- sigma * e0 - (0 + gamma) * i0
    dR <- gamma * i0 - 0 * r0
    list(c(dS, dE, dI, dR))
}

mle.seir <- function(b, g, s) {
    t <- seq(0, 14)    
    beta <- exp(b)
    gamma <- exp(g)
    sigma <- exp(s)

    results <- as.data.frame(ode(y=c(S=760, E=2, I=1, R=0), 
                             times=t, 
                             closed.seir,
                             parms=c(beta, gamma, sigma)))
    nll <- -sum(dpois(x=data$cases, lambda=tail(results$I, 14),  log=TRUE))
    return(nll)
    #return(results)
}

initial <- list(b=-6, g=-0.7, s=-0.9)
fit0.seir <- mle2(mle.seir, start=initial)

fit.seir <- mle2(mle.seir, start=as.list(coef(fit0.seir)))

pred.seir <- as.data.frame(ode(c(S=762, E=2, I=1, R=0), times=seq(0, 14, 0.5), closed.seir, parms=c(exp(coef(fit.seir)))))

plot(cases~day, data=data, type="b", ylab="Proportion Infected")
lines(pred.seir$I~pred.seir$time)

seir.beta <- as.numeric(exp(coef(fit.seir)[1]))
seir.gamma <- as.numeric(exp(coef(fit.seir)[2]))
seir.sigma <- as.numeric(exp(coef(fit.seir)[3]))
seir.nll <- as.numeric(logLik(fit.seir))
seir.aic <- AIC(fit.seir)