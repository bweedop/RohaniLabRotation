#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from sympy import *
import pandas as pd

####################
# Closed SIR Model #
####################

# Setting up the time parameters

# time zero
t0 = 0
# end time
tn = 14
# time step
t_step = 1

time_steps = range(t0, tn, t_step)

# Parameters for the model
# Instantaneous, constant params

# Transmission rate
beta = 1.66
# Recovery rate
gamma = 0.4545

# State variables - initial

# Proportion of susceptibles
s0 = 752/753
# Proportion of infecteds
i0 = 1/753
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
data = pd.read_csv("data/boarding-school.csv")
plt.scatter(data['day'], data['cases']/753)
plt.title("Simple SIR without demographics\n" +r"($\beta$: {}, $\gamma$: {})".format(beta, gamma))
plt.xlabel("time (days)")
plt.ylabel("Proportion of Population")
plt.legend(("Susceptibles", "Infecteds", "Recovered"),
           loc = "upper right")
plt.savefig("figures/closed-sir.png")

##########################################
# Multiple values of beta and gamma plot #
##########################################

# Transmission rate
beta_values = [x for x in range(1, 5)]
# Recovery rate
gamma_values = [0.3, 0.1, 0.02, 0.03]

# time zero
t0 = 0
# end time
tn = 25
# time step
t_step = 1

time_steps = range(t0, tn, t_step)

widths = [5] * 4
heights = [5] * 4
gs_kw = dict(width_ratios=widths, height_ratios=heights)

fig, axs = plt.subplots(4, 4, gridspec_kw=gs_kw)
fig.subplots_adjust(wspace = 0.7, hspace = 0.7)
fig.suptitle("SIR over multiple values")
for beta in beta_values:
    for gamma in gamma_values:
        sol = odeint(closed_sir, init_params, time_steps)
        axs[beta_values.index(beta), gamma_values.index(gamma)].plot(sol)
        axs[beta_values.index(beta), gamma_values.index(gamma)].set_title(r'$\beta$:{} $\gamma$:{}'.format(beta, gamma))
plt.figure(num = 2, figsize=(9,6))
#plt.savefig("mutliple-plots.png", dpi = 300)