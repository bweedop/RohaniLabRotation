import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

# Establish range of time
t0 = 1
t_n = 50*365
t_step = 1
t_range = range(t0, t_n, t_step)

# Set up constant parameters
beta = 1.56
gamma = 0.1428
mu = 0.0000391
sigma = 0.071

# Set up initial proportions of classes
s0 = 0.1
e0 = 0.0001
i0 = 0.0001
init_params = (s0, e0, i0)

def latent_sir (params, t):
    res = np.zeros((3))
    res[0] = mu - (beta * params[2] + mu) * params[0]
    res[1] = beta * params[0] * params[2] - (mu + sigma) * params[1]
    res[2] = sigma * params[1] - (mu + gamma) * params[2]
    return res

sol = odeint(latent_sir, init_params, t_range)

recovered = 1.0 - (sol[:,0]+sol[:,1]+sol[:,2])

# Plotting results
fig, axs = plt.subplots(4, sharex=True)
fig.suptitle('SEIR Model With Demographics\n' + r'($\beta$: {}, $\gamma$: {}, $\sigma$: {}, $\mu$: {})'.format(beta, gamma, sigma, mu), fontsize = 10)
axs[0].title.set_text("Susceptible")
axs[0].plot(sol[:,0], color = "green")
axs[1].title.set_text("Exposed")
axs[1].plot(sol[:,1], color = "yellow")
axs[2].title.set_text("Infected")
axs[2].plot(sol[:,2], color = "red")
axs[3].title.set_text("Recovered")
axs[3].plot(recovered, color = "blue")
axs[3].set_xlabel('time (d)')
plt.subplots_adjust(hspace = 0.4)
plt.figure(num = 1, figsize=(6,9))
plt.savefig("latent-sir.png", dpi = 600)