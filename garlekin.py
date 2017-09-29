import numpy as np
import matplotlib.pyplot as plt


n = 40

alpha = 1
beta  = 0
f = np.zeros(n-1) + 2*alpha




h = np.array([1.0/n]*n)
x = np.linspace(h[0], 1.0, num=n-1, endpoint=False)

sol_analitica = np.apply_along_axis(lambda x: x - x*x, 0, x)
f = sol_analitica.copy()
beta=1; alpha=0;

print f


K = np.zeros((n-1, n-1))

for i in xrange(n-1):
    K[i,i] = alpha*(1/h[i] + 1/h[i+1]) + beta*(h[i] + h[i+1])/3

    if i < n-2 :
        K[i,i+1] = -alpha/h[i+1] + beta*h[i+1]/6

    if i > 0:
        K[i,i-1] = -alpha/h[i] + beta*h[i]/6



F = np.zeros(n-1)


for i in xrange(n-1):

    F[i] = f[i]*h[i]/3 + f[i]*h[i+1]/3

    if i > 0:
        F[i] += f[i-1]*h[i]/6 

    if i < n - 2:
        F[i] += f[i+1]*h[i+1]/6




y = np.linalg.solve(K, F)


plt.plot(x,y)
plt.plot(x, sol_analitica)
plt.show()