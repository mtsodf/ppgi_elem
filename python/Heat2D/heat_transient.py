#!/usr/bin/env python
# -*- coding: utf-8 -*-

from heat import ConstructCase
from FiniteElement import *
from PlotMap import *
import argparse
import glob
import os


def SolveSedo(M, K, F, d0, nsteps, alpha, dt):
    v0 = np.linalg.solve(M, F - np.dot(K, d0))
    
    dpreditor = d0 + np.dot(dt * (1 - alpha), v0)

    v0 = np.linalg.solve(M + np.dot(alpha * dt, K) , F - np.dot(K, dpreditor))
    d0 = dpreditor + np.dot(alpha * dt, v0)


    return d0.copy()


def SolveSedoNewtonImplicit(M, K, F, d0, nsteps, dt):

    K = K * dt
    F = F * dt

    dcurrent = d0

    dprevious = dcurrent

    dcurrent = np.linalg.solve(M + K, F + np.dot(M, dprevious))

    return dcurrent


def main():

    if not os.path.isdir("out"):
        os.mkdir("out")

    os.chdir("out")
    for f in glob.glob("step_*.png"):
        os.remove(f)

    parser = argparse.ArgumentParser(description='Transferencia de Calor 2D')
    parser.add_argument('--entrada', type=int, default=5,
                        help='quantidade de elementos na direção x')
    parser.add_argument('--nx', type=int, default=20,
                        help='quantidade de elementos na direção x')
    parser.add_argument('--ny', type=int, default=20,
                        help='quantidade de elementos na direção y')
    parser.add_argument('--dx', type=float, default=0.1,
                        help='tamanho da célula na direção x')
    parser.add_argument('--dy', type=float, default=0.1,
                        help='tamanho da célula na direção y')
    parser.add_argument('--alpha', type=float, default=0.5,
                        help='alpha para método de solução')
    parser.add_argument('--dt', type=float, default=0.002,
                        help='dt para método de solução')
    parser.add_argument('--nsteps', type=int, default=40,
                        help='numero de passos de tempo')

    parser.add_argument('-n', '--newton', action="store_true")
    parser.add_argument('-p', '--plot3D', action="store_true")

    args = parser.parse_args()

    nx = args.nx
    ny = args.ny
    entrada = args.entrada
    newton = args.newton
    alpha = args.alpha
    dt = args.dt
    nsteps = args.nsteps
    plot3D = args.plot3D

    elements, nodes, neq, sol = ConstructCase(entrada, nx, ny, verbose=False)


    nelem = len(elements)

    def initialCondition(x, y):
        return sol(x, y, t=0.0)

    # ***************************************************************
    #                   Calculo do Lado Direito
    # ***************************************************************
    F = CalcF(elements, neq, t=0.0)


    # ***************************************************************
    #                Construindo Matriz de Rigidez
    # ***************************************************************

    K = BuildStiffness(elements, neq)
    print "Tamanho K -> ", K.shape

    M = BuildM(elements, neq)
    print "Tamanho M -> ", M.shape

    d0 = np.zeros(neq)

    for node in nodes:
        if node.eq is not None:
            d0[node.eq] = initialCondition(node.coords[0], node.coords[1])

    X = np.array([node.coords[0] for node in nodes if node.eq is not None])
    Y = np.array([node.coords[1] for node in nodes if node.eq is not None])

    sols = [d0]

    t = 0.0
    for step in range(1, nsteps + 1):
        t = step*dt
        
        print "Calcundo time step %d. Tempo: %f." % (step, t)

        F = CalcF(elements, neq, t)
        if newton:
            print "Utilizando metodo de newton"
            sols.append(SolveSedoNewtonImplicit(M, K, F, d0, nsteps, dt))
        else:
            sols.append(SolveSedo(M, K, F, d0, nsteps, alpha, dt))

        d0 = sols[-1]

    zmax = np.amax(sols[0])
    zmin = np.amin(sols[0])

    for s in sols:
        aux = np.amax(s)
        zmax = aux if aux > zmax else zmax
        zmin = aux if aux < zmin else zmin

    for step in range(nsteps + 1):

        print "Plotando time step %d" % step

        fig = plt.figure(figsize=(16, 9))

        X, Y, Z = GetGridValues(nodes, sols[step])

        t = step*dt

        Zsol = np.zeros(len(Z))

        if plot3D:

            ax = fig.add_subplot(121, projection='3d')
            ax2 = fig.add_subplot(122, projection='3d')

            plot_surface(X, Y, Z, ax=ax, fig=fig, zmin=zmin, zmax=zmax + zmax / 1000)

            for i in xrange(len(X)):
                Zsol[i] = sol(X[i], Y[i], t)

            print "Norma da diferenca = ", np.linalg.norm(Z-Zsol)/np.linalg.norm(Zsol)

            plot_surface(X, Y, Zsol, ax=ax2, fig=fig, zmin=zmin, zmax=zmax + zmax / 1000)

        else:
            ax = fig.add_subplot(111)
            plot_map(X, Y, Z, ax=ax, fig=fig,
                     zmin=zmin, zmax=zmax + zmax / 1000, nx=nx, ny=ny)
            plot_elements(elements, ax=ax)


        plt.suptitle("t = %6.4f" % (step * dt))

        plt.savefig('step_%s.png' % str(step).zfill(2))
        plt.close()

    import imageio
    with imageio.get_writer('movie.gif', mode='I') as writer:
        filenames = glob.glob("step_*.png")
        filenames.sort()
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)


if __name__ == '__main__':
    main()
