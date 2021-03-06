import pylab as plt

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import interpolate
from matplotlib import cm
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.mlab import griddata
from Element import *


def plot_solution(sol, context, uxAxis=None, uyAxis=None, fig=None):

    xArray = np.zeros(context.QtdNodes())
    yArray = np.zeros(context.QtdNodes())

    valuesUx = np.zeros_like(xArray)
    valuesUy = np.zeros_like(xArray)

    i = 0


    for node in context.nodesIter():
        x, y = node.coords

        xArray[i] = x
        yArray[i] = y

        eqx, eqy = context.getEq(node)

        px, py = context.getP(node)

        valuesUx[i] = px if eqx is None else sol[eqx]
        valuesUy[i] = py if eqy is None else sol[eqy]

        i+=1

    if uxAxis is not None:
        plot_map(xArray, yArray, valuesUx, ax=uxAxis, fig=fig)

    if uyAxis is not None:
        plot_map(xArray, yArray, valuesUy, ax=uyAxis, fig=fig)






def plot_surface(x, y, z, xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None, ax=None, fig=None, n=20, out='surface.png'):

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    xmin = np.amin(x) if xmin is None else xmin
    xmax = np.amax(x) if xmax is None else xmax

    ymin = np.amin(y) if ymin is None else ymin
    ymax = np.amax(y) if ymax is None else ymax

    # define grid.
    X = np.linspace(xmin, xmax, n, endpoint=True)
    Y = np.linspace(ymin, ymax, n, endpoint=True)

    X, Y = np.meshgrid(X, Y)
    # grid the data.
    Z = griddata(x, y, z, X, Y, interp='linear')

    zmin = np.amin(Z) if zmin is None else zmin
    zmax = np.amax(Z) if zmax is None else zmax

    Z = np.reshape(Z, (n, n))

    ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False, vmin=zmin, vmax=zmax)
    ax.set_zlim(zmin, zmax)

    plt.savefig(out)

def plot_map(x, y, z, xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None, ax=None, fig=None, nx=20, ny=20, out=None):

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    xmin = np.amin(x) if xmin is None else xmin
    xmax = np.amax(x) if xmax is None else xmax

    ymin = np.amin(y) if ymin is None else ymin
    ymax = np.amax(y) if ymax is None else ymax

    # define grid.
    X = np.linspace(xmin, xmax, nx, endpoint=True)
    Y = np.linspace(ymin, ymax, ny, endpoint=True)
    # grid the data.
    Z = griddata(x, y, z, X, Y, interp='linear')

    zmin = np.amin(Z) if zmin is None else zmin
    zmax = np.amax(Z) if zmax is None else zmax

    Z = np.reshape(Z, (ny, nx))

    im = ax.imshow(Z, interpolation="bilinear", cmap='rainbow', origin='lower',
                   extent=[xmin, xmax, ymin, ymax], vmin=zmin, vmax=zmax)

    if fig is not None:
        fig.colorbar(im, orientation ='horizontal')
    if out is not None:
        plt.savefig(out)

def plot_elements(elements, ax=None):

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    for element in elements:
        xs = []
        ys = []
        for node in element.nodes:
            xs.append(node.coords[0])
            ys.append(node.coords[1])

        ax.plot(xs, ys, 'ko-', markersize=0)


if __name__ == '__main__':

    t = 1.0

    n = 50

    X = np.linspace(0.0, 2.0, n, endpoint=True)
    Y = np.linspace(0.0, 1.5, n, endpoint=True)

    x, y = np.meshgrid(X, Y)

    z = np.sin(np.pi * x) * np.sin(np.pi * y)

    plot_surface(x, y, z)

    plot_map(x, y, z)
