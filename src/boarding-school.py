#!/usr/bin/python3

###################################################
# Closed SIR Model Fitted to Boarding School data #
###################################################

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import curve_fit, least_squares, leastsq  
from tabulate import tabulate
import pandas as pd
from numpy import exp, sin  

pop_size = 763

data = pd.read_csv("./data/boarding-school.csv")

i_data = np.array(data['cases']/pop_size, dtype=float)
t_data = np.array(data['day'], dtype=float)

N = 1.0
I0 = i_data[0]
S0 = N - I0
R0 = 0.0

def closed_sir (y, x, beta, gamma):
    S = -beta * y[0] * y[1]
    I = beta * y[0] * y[1] - gamma * y[1]
    R = gamma * y[1]
    return S, I, R

def fit_odeint(x, beta, gamma):
    return odeint(closed_sir, (S0, I0, R0), x, args=(beta, gamma))[:,1]

popt, pcov = curve_fit(fit_odeint, t_data, i_data)

beta = popt[0]
gamma = popt[1]

sol = odeint(closed_sir, (S0, I0, R0), t_data, args=(beta, gamma))

# Residual sum of squares
residuals = (data['cases']/763 - sol[:,1])**2
rss = np.sum(residuals)

# Plot the results in one figure
plt.plot(sol)
plt.scatter(t_data, i_data)
plt.title("Simple SIR without demographics\n" +r"($\beta$: {}, $\gamma$: {})".format(round(beta, 2), round(gamma, 2)))
plt.xlabel("time (days)")
plt.ylabel("Proportion of Population")
plt.legend(("Susceptibles",  "Infecteds", "Recovered"),
           loc = "upper right")
plt.savefig("figures/closed-sir.png")
plt.clf()