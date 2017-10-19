from ctypes import *
import matplotlib.pyplot as plt
import numpy as np
from math import pi, sin, exp

def get_sol(x, entrada):

    if entrada == 0 or entrada  == 8:
        alpha = 1.0
        beta  = 1.0
        f = np.vectorize(lambda x: sin(2*pi*x))
    
    if entrada == 2:
        f = np.vectorize(lambda x: x)

    if entrada == 1 or entrada == 3:
        f = np.vectorize(lambda x: x - x*x)
    
    if entrada == 4 or entrada == 5 or entrada == 6 or entrada == 7:
        f = np.vectorize(lambda x: exp(2*x))
    
    return f(x)

lib = cdll.LoadLibrary('./libFem1d.dylib')


lib.Fem1dTest.restype = POINTER(c_float)

import argparse

parser = argparse.ArgumentParser(description='Metodo de Garlekin')
parser.add_argument('--entrada', type=int, default=0, help='entrada do problema')
parser.add_argument('--pontos', type=int, default=20, help='entrada do problema')
args = parser.parse_args()




n = args.pontos
entrada = args.entrada

calc_c = lib.Fem1dTest(n, entrada)
h = np.array([1.0/n]*n)
x = np.zeros(n+1)
x[0] = 0.0

for i in xrange(1,n+1):
    x[i] = x[i-1] + h[i-1]

calc = []
for i in range(n+1):
    calc.append(calc_c[i])

sol = get_sol(x, entrada)


fig, ax = plt.subplots(1,2,figsize=(15,8))
ax[0].plot(x,calc, label="Calculada")
ax[0].scatter(x, sol, label="Analitica", color="red")
ax[0].legend()
ax[1].plot(x, np.abs(calc-sol))
ax[0].set_title("Solucao")
ax[1].set_title("Erro")

plt.show()
