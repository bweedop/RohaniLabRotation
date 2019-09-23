#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from sympy import *

####################
# Closed SIR Model #
####################

# Setting up the time parameters

# time zero
t0 = 0
# end time
tn = 120
# time step
t_step = 5

time_steps = range(t0, tn, t_step)

# Parameters for the model
# Instantaneous, constant params

# Transmission rate
beta = 0.3
# Recovery rate
gamma = 1/7

# State variables - initial

# Proportion of susceptibles
s0 = 9999/10000
# Proportion of infecteds
i0 = 1/10000
# Proportion of recovereds
r0 = 0.0

init_params = (s0, i0, r0)

def closed_sir (params, time):
    res = np.zeros((3))
    res[0] = -beta * params[0] * params[1]
    res[1] = beta * params[0] * params[1] - gamma * params[1]
    res[2] = gamma * params[1]
    return res

# Use scipy package to integrate solutions
sol = odeint(closed_sir, init_params, time_steps)

# Plot the results in one figure
plt.plot(sol)
plt.xlabel("time (days)")
plt.ylabel("Proportion of Population")
plt.legend(("Susceptibles", "Infecteds", "Recovered"),
           loc = "upper right")
plt.savefig("closed-sir.png")
