#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import least_squares, curve_fit
import pandas as pd

################################################################################
# Closed SEIR Model for Boarding School Data                                   #
################################################################################

pop_size = 763

data = pd.read_csv("./data/boarding-school.csv")

i_data = np.array(data['cases']/pop_size, dtype=float)
t_data = np.array(data['day'], dtype=float)

N = 1.0
I0 = i_data[0]
S0 = N - I0
E0 = I0
R0 = 0.0

def closed_sir (y, x, beta, gamma, sigma):
    S = -beta * y[0] * y[2]
    E = beta * y[0] * y[2] - sigma * y[1]
    I = sigma * y[1] - gamma * y[2]
    R = gamma * y[2]
    return S, E, I, R

def fit_odeint(x, beta, gamma, sigma):
    return odeint(closed_sir, (S0, E0, I0, R0), x, args=(beta, gamma, sigma))[:,2]

popt, pcov = curve_fit(fit_odeint, t_data, i_data)

beta = popt[0]
gamma = popt[1]
sigma = popt[2]

# Use scipy package to integrate solutions
sol = odeint(closed_sir, (S0, E0, I0, R0), t_data, args=(beta, gamma, sigma))

# Plot the results in one figure
plt.plot(sol)
plt.scatter(t_data, i_data)
plt.title("Simple SIR without demographics\n" +r"($\beta$: {}, $\gamma$: {})".format(round(beta, 2), round(gamma, 2)))
plt.xlabel("time (days)")
plt.ylabel("Proportion of Population")
plt.legend(("Susceptibles", "Exposed", "Infecteds", "Recovered"),
           loc = "upper right")
plt.savefig("figures/exposed-sir.png")
plt.clf()