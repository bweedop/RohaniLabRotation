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
t_step = 1

time_steps = range(t0, tn, t_step)

# Parameters for the model
# Instantaneous, constant params

# Transmission rate
beta = 1.65
# Recovery rate
gamma = 0.4545

# State variables - initial

# Proportion of susceptibles
s0 = 0.9999
# Proportion of infecteds
i0 = 0.0001
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
plt.title("Simple SIR without demographics (Beta: {}, Gamma: {})".format(beta
                                                                         , gamma))
plt.xlabel("time (days)")
plt.ylabel("Proportion of Population")
plt.legend(("Susceptibles", "Infecteds", "Recovered"),
           loc = "upper right")
plt.savefig("closed-sir.png")

##########################################
# Multiple values of beta and gamma plot ##
##########################################

# Transmission rate
#beta_values = [2**x for x in range(5, 9)]
#Recovery rate
#gamma_values = [0.33, 0.11, 0.055, 0.037]

##fig, axs = plt.subplots(4, 4)
#fig.suptitle("SIR over multiple values")
#for beta in beta_values:
#    for gamma in gamma_values:
#        sol = odeint(closed_sir, init_params, time_steps)
#        axs[beta_values.index(beta), gamma_values.index(gamma)].plot(sol)
#        axs[beta_values.index(bet#a), gamma_values.index(gamma)].set_title('Beta:{} Gamma:{}'.format(beta, gamma))
#plt.savefig("mutliple-plots.png")
