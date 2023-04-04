#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 15:26:51 2023

@author: carlosmh1
"""

import pandas as pd
import numpy as np
import os
from scipy.optimize import fsolve
from scipy.stats import norm

# FUNCTIONS DECLARATION

def clear_screen():
    if os.name == 'nt':  # For Windows
        os.system('cls')
    else:  # For macOS and Linux
        os.system('clear')

def fun(x, t, g, c):
    z = np.exp(-x*t)
    return np.dot(z, g) - c

# MODIFY PARAMETERS AS NEEDED
alpha = 0.05
fnam = "test_data.csv"

clear_screen()
Z = norm.ppf(1-alpha/2)



# IMPORT DATA
df = pd.read_csv(fnam, header=0, skip_blank_lines=True, usecols=[0, 1, 2])
df.dropna(inplace=True)

t = np.array(df['t'])
n = np.array(df['n'])
h = np.array(df['h'])

# CALCULATE VALUES
N = n[0]
K = np.sum(h)
g = h/K
f = (n[:-1] - n[1:])/N
f = np.append(f, 0)

# CREATE A NEW DATAFRAME
data = pd.DataFrame({"t": t, "n": n, "h": h, "g": g, "f": f})

# SAVE THE NEW DATAFRAME
fnam2 = fnam.split('.')[0] + '_added.csv'
data.to_csv(fnam2, index=False)

# CALCULATE R0
R0 = K / N

# CALCULATE LONGEVITY
EL = np.dot(t, f)
EL2 = np.dot(t * t, f)
VL = EL2 - EL * EL
CI_L = np.array([EL, VL, EL - Z * np.sqrt(VL / N), EL + Z * np.sqrt(VL / N)])

# CALCULATE GENERATION TIME
ET = np.dot(t, g)
ET2 = np.dot(t * t, g)
VT = ET2 - ET * ET
CI_mu = np.array([ET, VT, ET - Z * np.sqrt(VT / K), ET + Z * np.sqrt(VT / K)])

# CALCULATE POPULATION GROWTH RATE, r
x_initial_guess = 0.05
r_sol = fsolve(fun, x_initial_guess, args=(t, g, 1 / R0), xtol=1e-6)

mu = 1 / R0
z = np.exp(-(2 * r_sol) * t)
s2 = np.dot(z, g) - mu ** 2
s = np.sqrt(s2)

r_L_sol = fsolve(fun, r_sol, args=(t, g, mu + Z * s / np.sqrt(K)), xtol=1e-6)
r_U_sol = fsolve(fun, r_sol, args=(t, g, mu - Z * s / np.sqrt(K)), xtol=1e-6)

# GROWTH RATES
CI_r = [r_sol[0], r_L_sol[0], r_U_sol[0]]
CI_lam = [np.exp(i) for i in CI_r]

# DISPLAY RESULTS


print('Initial number of individuals N : ', N)
print('Offspring size K : ', K)
print('R0 : ', R0)
print('Longevity : ', CI_L)
print('Generation time : ', CI_mu)
print('r : ', CI_r)
print('lambda: ', CI_lam)
print('New data saved to: ', fnam2)

