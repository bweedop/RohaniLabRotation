#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

# Constant parameters
beta = 1.3
gamma = 0.333
mu = 0.0003

# Setting up time range
t0 = 1
tn = 50*365
t_step = 1
t_span = range(t0, tn, t_step)
plot_time = range(1, 50)

# Initial state parameter values
s0 = 0.1
i0 = 0.0001
r0 = 1 - s0 - i0
init_params = (s0, i0, r0)

# Differential equations function
def demographics_sir (params, t):
    res = np.zeros((3))
    res[0] = mu - beta * params[0] * params[1] - mu * params[0]
    res[1] = beta * params[0] * params[1] - gamma * params[1] - mu * params[1]
    res[2] = gamma * params[1] - mu * params[2]
    return res

# Solution using the odeint function provided by scipy
sol = odeint(demographics_sir, init_params, t_span)

# Plotting results
fig, axs = plt.subplots(3, sharex=True)
fig.suptitle(r'SIR Model With Demographics ($\beta$: {}, $\gamma$: {})'.format(beta, gamma))
axs[0].title.set_text("Susceptible")
axs[0].plot(sol[:,0], color = "green")
axs[1].title.set_text("Infected")
axs[1].plot(sol[:,1], color = "red")
axs[2].title.set_text("Recovered")
axs[2].plot(sol[:,2], color = "orange")
plt.subplots_adjust(hspace = 0.4)
plt.savefig("demographics-sir.png")
