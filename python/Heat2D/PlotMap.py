import pylab as plt

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import interpolate
from matplotlib import cm


def plot_surface(x, y, z, ax=None, n=20):

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

    finter = interpolate.interp2d(x, y, z)

    X = np.linspace(np.amin(x), np.amax(x), n, endpoint=True)
    Y = np.linspace(np.amin(y), np.amax(y), n, endpoint=True)
    Z = np.zeros(n * n)

    c = 0
    for j in range(n):
        for i in range(n):
            Z[c] = finter(X[i], Y[j])
            c += 1

    X, Y = np.meshgrid(X, Y)
    Z = np.reshape(Z, (n, n))

    ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)

    plt.show()


def plot_map(x, y, z, ax=None, fig=None, n=20):

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    xmin = np.amin(x)
    xmax = np.amax(x)

    ymin = np.amin(y)
    ymax = np.amax(y)

    X = np.linspace(xmin, xmax, n, endpoint=True)
    Y = np.linspace(ymin, ymax, n, endpoint=True)
    X, Y = np.meshgrid(X, Y)

    finter = interpolate.interp2d(x, y, z)

    X = np.linspace(np.amin(X), np.amax(X), n, endpoint=True)
    Y = np.linspace(np.amin(Y), np.amax(Y), n, endpoint=True)
    Z = np.zeros(n * n)

    c = 0
    for j in range(n):
        for i in range(n):
            Z[c] = finter(X[i], Y[j])
            c += 1

    X, Y = np.meshgrid(X, Y)
    zmin = np.amin(Z)
    zmax = np.amax(Z)

    Z = np.reshape(Z, (n, n))

    im = ax.imshow(Z, interpolation="gaussian", cmap='seismic', origin='lower',
                   extent=[xmin, xmax, ymin, ymax], vmin=zmin, vmax=zmax)

    fig.colorbar(im)
    plt.show()


if __name__ == '__main__':

    t = 1.0

    n = 50

    X = np.linspace(0.0, 2.0, n, endpoint=True)
    Y = np.linspace(0.0, 1.5, n, endpoint=True)

    x, y = np.meshgrid(X, Y)

    z = np.sin(np.pi * x) * np.sin(np.pi * y)

    plot_surface(x, y, z)

    plot_map(x, y, z)
