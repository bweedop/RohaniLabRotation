import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import pandas as pd
import random as rd

rd.seed(1881)

def sir_onestep (x, params):
    num_sus = x[1]
    num_inf = x[2]
    