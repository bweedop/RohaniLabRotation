#!/usr/bin/R
packages <- c("deSolve")
packages.to.install <- packages[!(packages %in% installed.packages()[,'Package'])]
if (length(packages.to.install)) {
    install.packages(packages.to.install)
}
invisible(lapply(packages, require, character.only = TRUE))

################################################################################
# SIR with Sinusoidal Forcing (CH 5)                                           #
################################################################################
betafx <- function(beta0 = 1.65, beta1 = 0.45, t) {
    omega <- (2*pi)/365
    return(beta0 * (1 + beta1 * cos(omega * t)))
}

seasonal.sir <- function(t, states, params) {
    s0 <- states[1]
    i0 <- states[2]
    r0 <- states[3]

    beta <- betafx(beta0 = params[1], beta1 = params[2], t = t)
    mu <- 0

    dS <- mu - beta * s0 * i0 - mu * s0
    dI <- beta * s0 * i0 - gamma * i0 - mu * i0
    dR <- gamma * i0 - mu * r0
    return(list(c(dS, dI, dR)))
}

s0 <- 0.9
i0 <- 0.0001
r0 <- 1-s0-i0

t_range <- seq(0, 50*365, 1)

initial.states <- c(S = s0, I = i0, R = r0)

params <- c(beta0=1.4, beta1=0.1, gamma = 1/13)

out <- ode(y = initial.states, times = t_range, func = seasonal.sir, parms = params)