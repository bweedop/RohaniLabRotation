#!/usr/bin/R

################################################################################
# boarding-school-models.R - This script is intended to fit a number of        #
#                            epidemiological models to real data (boarding     #
#                            school, London, 1978). The model coeficients will #
#                            be estimated using maximum likelihood as          #
#                            implemented by the 'bbmle' package.               #
################################################################################

packages <- c("xtable", "deSolve", "data.table", "bbmle")
packages.to.install <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(packages.to.install)) {
    install.packages(packages.to.install)
}
invisible(lapply(packages, require, character.only = TRUE))

# Load the boarding school data
flu.data <- read.csv("data/boarding-school.csv")

################################################################################
# SIR model with no vital dynamics                                             #
# Beta and gamma estimated using maximum likelihood                            #
################################################################################

# Closed (no vital dynamics) SIR model
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

# mle.sir - Maximum likelihood estimation function for closed SIR model
mle.sir <- function(b, g) {
    t <- seq(0, 14)    
    beta <- exp(b)
    gamma <- exp(g)

    results <- as.data.frame(ode(y=c(S=762, I=1, R=0), 
                             times=t, 
                             closed.sir,
                             parms=c(beta, gamma)))
    nll <- -sum(dpois(x=flu.data$cases, lambda=tail(results$I, 14), log=TRUE))
    return(nll)
    #return(results)
}

# initial - Initial estimates of beta and gamma
initial <- list(b=-6, g=-0.7)
# fit0 - preliminary fit of model to data using the initial estimates
fit0 <- mle2(mle.sir, start=initial)
# fit - subsequent fit using the estimation of parameters returned from fit0
fit <- mle2(mle.sir, start=as.list(coef(fit0)))
# pred - model estimates using the estimated beta and gamma values
pred <- as.data.frame(ode(c(S=762, I=1, R=0), 
                      times=seq(0, 14, 0.5), 
                      closed.sir, 
                      parms=c(exp(coef(fit)))))

# Plotting raw data
plot(cases~day, data=flu.data, type="b", ylab="Numbers of Infecteds")
# Plotting model predictions
lines(pred$I~pred$time, col="red")

# Final estimates of beta and gamma
sir.beta <- as.numeric(exp(coef(fit)[1])*763)
sir.gamma <- as.numeric(exp(coef(fit)[2]))
# Final R0 estimation
sir.r0 <- sir.beta/sir.gamma
# Final negative log-likelihood
sir.nll <- as.numeric(logLik(fit))
# AIC score for simple SIR model
sir.aic <- AIC(fit)

################################################################################
# SEIR Model with no vital dynamics (mu = 0)                                   #
# Beta, gamma, sigma estimated using maximum likelihood                        #
################################################################################

# closed.seir - Initializing SEIR model to be used with deSolve
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

# mle.seir - function to use for maximum likelihood estimation. Returns negative 
#            log-likelihood which is used to optimize coefficients 
mle.seir <- function(b, g, s) {
    t <- seq(0, 14)    
    beta <- exp(b)
    gamma <- exp(g)
    sigma <- exp(s)

    results <- as.data.frame(ode(y=c(S=761, E=1, I=1, R=0), 
                             times=t, 
                             closed.seir,
                             parms=c(beta, gamma, sigma)))
    nll <- -sum(dpois(x=flu.data$cases, lambda=tail(results$I, 14), log=TRUE))
    return(nll)
}

initial <- list(b=-5, g=-0.45, s=-0.45)
fit0.seir <- mle2(mle.seir, start=initial)

fit.seir <- mle2(mle.seir, start=as.list(coef(fit0.seir)))

pred.seir <- as.data.frame(ode(c(S=761, E=1, I=1, R=0), 
                               times=seq(0, 14, 0.5), 
                               closed.seir,
                               parms=c(exp(coef(fit.seir)))))

plot(cases~day, data=flu.data, type="b", ylab="Number of Infecteds")
lines(pred.seir$I~pred.seir$time, col="red")

seir.beta <- as.numeric(exp(coef(fit.seir)[1])*763)
seir.gamma <- as.numeric(exp(coef(fit.seir)[2]))
seir.sigma <- as.numeric(exp(coef(fit.seir)[3]))
seir.r0 <- seir.beta/seir.gamma
seir.nll <- as.numeric(logLik(fit.seir))
seir.aic <- AIC(fit.seir)

################################################################################
# SICR with no vital dynamics (mu = 0)                                         #
################################################################################

# closed.sicr - Initializing sicr model to be used with deSolve
closed.sicr <- function(t, states, params) {
    s0 <- states[1]
    i0 <- states[2]
    c0 <- states[3]
    r0 <- states[4]

    beta <- params[1]
    gamma <- params[2]
    zeta <- params[3]
    eta <- params[4]
    mu <- 0

    dS <- mu - beta * s0 * i0 - mu * s0
    dI <- beta * s0 * i0 - (mu + gamma + eta) * i0
    dC <- eta * i0 - zeta * c0
    dR <- gamma * i0 + zeta * c0 - mu * r0
    list(c(dS, dI, dC, dR))
}

# mle.sicr - function to use for maximum likelihood estimation. Returns negative 
#            log-likelihood which is used to optimize coefficients 
mle.sicr <- function(b, g, z, e) {
    t <- seq(0, 14)    
    beta <- exp(b)
    gamma <- exp(g)
    zeta <- exp(z)
    eta <- exp(e)

    results <- as.data.frame(ode(y=c(S=762, I=1, C=0, R=0), 
                             times=t, 
                             closed.sicr,
                             parms=c(beta, gamma, zeta, eta)))
    nll <- -sum(dpois(x=flu.data$cases, lambda=tail(results$C, 14), log=TRUE))
    return(nll)
}

initial <- list(b=-5.76, g=-12, z=-0.78, e=-0.0857)
fit0.sicr <- mle2(mle.sicr, start=initial)

fit.sicr <- mle2(mle.sicr, start=as.list(coef(fit0.sicr)))

pred.sicr <- as.data.frame(ode(c(S=762, I=1, C=0, R=0), 
                               times=seq(0, 14, 0.5), 
                               closed.sicr,
                               parms=c(exp(coef(fit.sicr)))))

png("sicr-model.png")
plot(cases~day, data=flu.data, type="b", ylab="Number of Infecteds")
lines(pred.sicr$C~pred.sicr$time, col="red")
dev.off()

sicr.beta <- as.numeric(exp(coef(fit.sicr)[1])*763)
sicr.gamma <- as.numeric(exp(coef(fit.sicr)[2]))
sicr.zeta <- as.numeric(exp(coef(fit.sicr)[3]))
sicr.eta <- as.numeric(exp(coef(fit.sicr)[4]))
sicr.r0 <- sicr.beta/(sicr.gamma + sicr.zeta)
sicr.nll <- as.numeric(logLik(fit.sicr))
sicr.aic <- AIC(fit.sicr)

################################################################################
# SICR with no vital dynamics (mu = 0)                                         #
################################################################################

closed.seicr <- function(t, states, params) {
    s0 <- states[1]
    e0 <- states[2]
    i0 <- states[3]
    c0 <- states[4]
    r0 <- states[5]

    beta <- params[1]
    gamma <- params[2]
    zeta <- params[3]
    eta <- params[4]
    sigma <- params[5]
    mu <- 0

    dS <- mu - beta * s0 * i0 - mu * s0
    dE <- beta * s0 * i0 - (mu + sigma) * e0
    dI <- sigma * e0 - (mu + gamma + eta) * i0
    dC <- eta * i0 - zeta * c0
    dR <- gamma * i0 + zeta * c0 - mu * r0
    list(c(dS, dE, dI, dC, dR))
}

# mle.seicr - function to use for maximum likelihood estimation. Returns negative 
#             log-likelihood which is used to optimize coefficients 
mle.seicr <- function(b, g, z, e, s) {
    t <- seq(0, 14)    
    beta <- exp(b)
    gamma <- exp(g)
    zeta <- exp(z)
    eta <- exp(e)
    sigma <- exp(s)

    results <- as.data.frame(ode(y=c(S=761, E=1, I=1, C=0, R=0), 
                             times=t, 
                             closed.seicr,
                             parms=c(beta, gamma, zeta, eta, sigma)))
    nll <- -sum(dpois(x=flu.data$cases, lambda=tail(results$C, 14), log=TRUE))
    return(nll)
}

initial <- list(b=-4.98, g=-10.71, z=-0.71, e=-0.34, s=-0.35)
fit0.seicr <- mle2(mle.seicr, start=initial)

fit.seicr <- mle2(mle.seicr, start=as.list(coef(fit0.seicr)))

pred.seicr <- as.data.frame(ode(c(S=761, E=1, I=1, C=0, R=0), 
                               times=seq(0, 14, 0.5), 
                               closed.seicr,
                               parms=c(exp(coef(fit.seicr)))))

png("seicr-model.png")
plot(cases~day, data=flu.data, type="b", ylab="Number of Infecteds")
lines(pred.seicr$C~pred.seicr$time, col="red")
dev.off()

seicr.beta <- as.numeric(exp(coef(fit.seicr)[1])*763)
seicr.gamma <- as.numeric(exp(coef(fit.seicr)[2]))
seicr.zeta <- as.numeric(exp(coef(fit.seicr)[3]))
seicr.eta <- as.numeric(exp(coef(fit.seicr)[4]))
seicr.r0 <- seicr.beta/(seicr.gamma + seicr.zeta)
seicr.nll <- as.numeric(logLik(fit.seicr))
seicr.aic <- AIC(fit.seicr)

################################################################################
# Create LateX table with xtable for all models and accompanying values        #
################################################################################

ml.df <- data.frame(Model=c("SIR", "SEIR", "SICR", "SEICR"), 
                    R0=c(sir.r0, seir.r0, sicr.r0, NA), 
                    nll=c(-sir.nll, -seir.nll, -sicr.nll, NA),
                    AIC=c(sir.aic, seir.aic, sicr.aic, NA),
                    beta=c(sir.beta, seir.beta, sicr.beta, NA),
                    gamma=c(sir.gamma, seir.gamma, sicr.gamma, NA),
                    zeta=c(NA, NA, sicr.zeta, NA),
                    sigma=c(NA, seir.sigma, NA, NA))

xtable(ml.df)