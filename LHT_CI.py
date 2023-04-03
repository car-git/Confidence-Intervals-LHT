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

alpha = 0.05
fnam = "test_data.csv"
dir_path = "/Users/carlosmh1/Library/CloudStorage/Dropbox/Current_Research/Estimation of LHT/"


Z = norm.ppf(1-alpha/2)
# print(Z)

def fun(x, t, g, c):
    z = np.exp(-x*t)
    return np.dot(z, g) - c


# create an empty pandas dataframe with column names
data = pd.DataFrame(columns=['t', 'n', 'h'])

t = pd.Series(dtype=int)
n = pd.Series(dtype=int)
h = pd.Series(dtype=int)
g = pd.Series(dtype='float64')
f = pd.Series(dtype='float64')

file = os.path.join(dir_path, fnam)

df = pd.read_csv(fnam, header=0,skip_blank_lines=True, usecols=[0, 1, 2])
df.dropna(inplace=True)

data['t'] = df.iloc[:, 0]
data['n'] = df.iloc[:, 1]
data['h'] = df.iloc[:, 2]
# print(data)

t = np.array(df['t'])
n = np.array(df['n'])
h = np.array(df['h'])

N = n[0]
K = np.sum(h)
g = h/K
f = (n[:-1] - n[1:])/N
f = np.append(f, 0)

data['g'] = g
data['f'] = f
# print(data)

name_parts = fnam.split('.')
fnam2 = name_parts[0]+'_added.csv'
data.to_csv(fnam2, index=False)
print("")
print("")
print("")
print('Initial number of individuals N : ', N)
print('Offspring size K : ', K)
# R0
R0 = K/N
print('R0 : ',R0)

# Longevity 
EL = np.dot(t, f)
EL2 = np.dot(t*t, f)
VL = EL2-EL*EL
CI_L = np.array([EL , VL, EL-Z*np.sqrt(VL/N),EL+Z*np.sqrt(VL/N)])
print('Longevity : ',CI_L)

# Generation time
ET = np.dot(t, g)
ET2 = np.dot(t*t, g)
VT = ET2-ET*ET

CI_mu = np.array([ET , VT, ET-Z*np.sqrt(VT/K),ET+Z*np.sqrt(VT/K)])
print('Generation time : ',CI_mu)

# Population growth rate, r

x_initial_guess = 0.05  # initial guess for x
r_sol = fsolve(fun, x_initial_guess, args=(t, g, 1/R0),xtol = 1e-6)

mu = 1/R0
z = np.exp(-(2*r_sol)*t)
s2 = np.dot(z, g) - mu**2
s = np.sqrt(s2)

r_L_sol = fsolve(fun, r_sol, args=(t, g, mu+Z*s/np.sqrt(K)),xtol = 1e-6)
r_U_sol = fsolve(fun, r_sol, args=(t, g, mu-Z*s/np.sqrt(K)),xtol = 1e-6)

# Growth rates
CI_r = [r_sol[0], r_L_sol[0], r_U_sol[0]]
print('r : ',CI_r)
CI_lam = [np.exp(i) for i in CI_r]
print('lambda: ',CI_lam)
print('New data saved to: ',fnam2)
