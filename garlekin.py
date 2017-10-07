import numpy as np
import matplotlib.pyplot as plt
from math import sin, pi
import argparse 




parser = argparse.ArgumentParser(description='Metodo de Garlekin')
parser.add_argument('--entrada', type=int, default=0, help='entrada do problema')
parser.add_argument('--pontos', type=int, default=20, help='entrada do problema')

args = parser.parse_args()

entrada = args.entrada
n = args.pontos

h = np.array([1.0/n]*n)

#h = [0.05, 0.05, 0.1, 0.1, 0.1, 0.2, 0.1, 0.1, 0.1, 0.05, 0.05]
#n = len(h)


x = np.zeros(n-1)
x[0] = h[0]

for i in xrange(1,n-1):
    x[i] = x[i-1] + h[i]

######################################################################
# Selecao da entrada e lado direito do problema
######################################################################

if entrada == 0:
    alpha = 1
    beta  = 1
    f = np.vectorize(lambda x: sin(2*pi*x))
    sol_analitica = f(x)

    f = lambda y: (alpha*4*pi*pi*sin(2*pi*y)+beta*sin(2*pi*y))
    f = np.vectorize(f)
    f = f(x)

    nome = "u(x) = sin(2 pi x) alpha = %f beta = %f" %( alpha, beta)

elif entrada == 1:
    alpha=1; beta=0; 

    sol_analitica = np.apply_along_axis(lambda x: x - x*x, 0, x)
    f = np.zeros(n-1) + 2*alpha
    nome = "u(x) = x - x*x alpha = %f beta = %f" %( alpha, beta)
    
elif entrada == 2:
    alpha=0; beta=1;
    sol_analitica = np.apply_along_axis(lambda x: x - x*x, 0, x)
    f = np.apply_along_axis(lambda x: x - x*x, 0, x)
    nome = "u(x) = x - x*x alpha = %f beta = %f" %( alpha, beta)


######################################################################
# Montagem da Matriz
######################################################################

K = np.zeros((n-1, n-1))

for i in xrange(n-1):
    K[i,i] = alpha*(1/h[i] + 1/h[i+1]) + beta*(h[i] + h[i+1])/3

    if i < n-2 :
        K[i,i+1] = -alpha/h[i+1] + beta*h[i+1]/6

    if i > 0:
        K[i,i-1] = -alpha/h[i] + beta*h[i]/6


######################################################################
# Montagem do lado direito
######################################################################
F = np.zeros(n-1)


for i in xrange(n-1):

    F[i] = f[i]*h[i]/3 + f[i]*h[i+1]/3

    if i > 0:
        F[i] += f[i-1]*h[i]/6 

    if i < n - 2:
        F[i] += f[i+1]*h[i+1]/6



######################################################################
# Solucao do Sistema Linear
######################################################################
y = np.linalg.solve(K, F)



######################################################################
# Plot da Solucao Analitica e da Solucao Calculada
######################################################################
fig, ax = plt.subplots(1,2,figsize=(15,8))
ax[0].plot(x,y, label="Calculada")
ax[0].scatter(x, sol_analitica, label="Analitica")
ax[0].legend()
ax[1].plot(x, np.abs(y-sol_analitica))

ax[0].set_title("Solucao")
ax[1].set_title("Erro")



plt.show()
