from ctypes import *
import matplotlib.pyplot as plt

lib = cdll.LoadLibrary('./libFem1d.dylib')


lib.Fem1dTest.restype = POINTER(c_float)

n = 100
entrada = 1

x = lib.Fem1dTest(n, entrada)

xs = []
for i in range(n-1):
    xs.append(x[i])



plt.plot(xs)
plt.show()