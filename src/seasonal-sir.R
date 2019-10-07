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

    gamma <- params[3]
    mu <- params[4]

    dS <- mu - betafx(beta0 = params[1], beta1 = params[2], t = t) * s0 * i0 - mu * s0
    dI <- betafx(beta0 = params[1], beta1 = params[2], t = t) * s0 * i0 - gamma * i0 - mu * i0
    dR <- gamma * i0 - mu * r0
    return(list(c(dS, dI, dR)))
}

s0 <- 0.06
i0 <- 0.001
r0 <- 1-s0-i0

t_range <- seq(0, 50*365, 1)

initial.states <- c(S = s0, I = i0, R = r0)

params <- c(beta0 = 1.6, beta1 = 0.05, gamma = 0.08, mu = 0.0000547945)

out <- ode(y = initial.states, times = t_range, func = seasonal.sir, parms = params)

png("figures/seasonal-sir.png", width = 14, height = 8, units = "in", res=300)
par(mfrow=c(3,1), pin = c(10, 1.5), mai = c(0, 0.6, 0, 0.1), mar=c(0.5, 4, 2, 0.5), oma = c(0.5, 0, 0, 0.1))
plot(out[, 2] ~ out[, 1], type="l", col = "blue", xaxt = "n", xlab = "", ylab = "S")
par(mar=c(0.5, 4, 0, 0.5))
plot(out[, 3] ~ out[, 1], type = "l", col = "red", xaxt = "n", xlab = "", ylab = "I")
par(mar=c(4, 4, 0.5, 0.5))
plot(out[, 4] ~ out[, 1], type = "l", col = "green", xlab = "time (days)", xaxt = "n", ylab = "R")
xticks <- seq(0, 50*365, 365)
axis(1, at=xticks)
dev.off()