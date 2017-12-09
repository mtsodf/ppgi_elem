import pylab as plt

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import interpolate
from matplotlib import cm
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.mlab import griddata


def plot_elements(elements, sol):

    import matplotlib.pyplot as plt

    vertices = []
    codes = []

    for element in elements:

        codes = [Path.MOVETO] + [Path.LINETO] * 3 + [Path.CLOSEPOLY]
        vertices = [(1, 1), (1, 2), (2, 2), (2, 1), (0, 0)]

    codes += [Path.MOVETO] + [Path.LINETO] * 2 + [Path.CLOSEPOLY]
    vertices += [(4, 4), (5, 5), (5, 4), (0, 0)]

    vertices = np.array(vertices, float)
    path = Path(vertices, codes)

    pathpatch = PathPatch(path, facecolor='None', edgecolor='green')

    fig, ax = plt.subplots()
    ax.add_patch(pathpatch)
    ax.set_title('A compound path')

    ax.dataLim.update_from_data_xy(vertices)
    ax.autoscale_view()

    plt.show()


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

    im = ax.imshow(Z, interpolation="nearest", cmap='seismic', origin='lower',
                   extent=[xmin, xmax, ymin, ymax], vmin=zmin, vmax=zmax)

    fig.colorbar(im)
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
        ax.plot(xs, ys, 'ko-', markersize=3)


if __name__ == '__main__':

    t = 1.0

    n = 50

    X = np.linspace(0.0, 2.0, n, endpoint=True)
    Y = np.linspace(0.0, 1.5, n, endpoint=True)

    x, y = np.meshgrid(X, Y)

    z = np.sin(np.pi * x) * np.sin(np.pi * y)

    plot_surface(x, y, z)

    plot_map(x, y, z)
