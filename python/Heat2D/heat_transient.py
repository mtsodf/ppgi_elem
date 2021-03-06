#!/usr/bin/env python
# -*- coding: utf-8 -*-

from heat import ConstructCase
from FiniteElement import *
from PlotMap import *
import argparse
import glob
import os
from numpy.linalg import norm


def SolveSedo(M, K, F, d0, nsteps, alpha, dt, v0):


    dpreditor = d0 + np.dot(dt * (1 - alpha), v0)

    v0[:] = np.linalg.solve(M + np.dot(alpha * dt, K), F - np.dot(K, dpreditor))
    d0 = dpreditor + np.dot(alpha * dt, v0)

    return d0.copy()


def SolveSedoNewtonImplicit(M, K, F, d0, nsteps, dt):

    K = K * dt
    F = F * dt

    dcurrent = d0

    dprevious = dcurrent

    dcurrent = np.linalg.solve(M + K, F + np.dot(M, dprevious))

    return dcurrent


def run_transient_case(entrada, nx, ny, dt, alpha, newton, triangles, nsteps, verbose=False, plot3D=False, plotdelta=1, gif=True):
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

    if verbose:
        print "Tamanho K -> ", K.shape

    M = BuildM(elements, neq)

    if verbose:
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
        t = step * dt

        # Atualizando Condicao de contorno com o tempo
        for node in nodes:
            if node.dirichletBoundary:
                node.p = sol(node.coords[0], node.coords[1], t)

        if verbose:
            print "Calcundo time step %d. Tempo: %f." % (step, t)

        F = CalcF(elements, neq, t)

        if newton:
            if verbose:
                print "Utilizando metodo de newton"
            sols.append(SolveSedoNewtonImplicit(M, K, F, d0, nsteps, dt))
        else:
            if step == 1:
                v0 = np.linalg.solve(M, F - np.dot(K, d0))

            sols.append(SolveSedo(M, K, F, d0, nsteps, alpha, dt, v0[:]))

        d0 = sols[-1]

    zmax = np.amax(sols[0])
    zmin = np.amin(sols[0])

    for s in sols:
        aux = np.amax(s)
        zmax = aux if aux > zmax else zmax
        zmin = aux if aux < zmin else zmin

    residues = []

    solarray = np.zeros(neq)

    for step in range(nsteps + 1):

        if step % plotdelta == 0:

            if verbose:
                print "Plotando time step %d" % step

            fig = plt.figure(figsize=(16, 9))

            t = step * dt

            # Atualizando Condicao de contorno com o tempo
            for node in nodes:
                if node.dirichletBoundary:
                    node.p = sol(node.coords[0], node.coords[1], t)

            X, Y, Z = GetGridValues(nodes, sols[step])


            Zsol = np.zeros(len(Z))


            for node in nodes:
                if node.eq is not None:
                    solarray[node.eq] = sol(node.coords[0], node.coords[1], t=t)

            residues.append(norm(sols[step] - solarray) / norm(solarray))

            if verbose:
                print "Norma da diferenca = ", residues[-1]

            if plot3D:

                ax = fig.add_subplot(121, projection='3d')
                ax2 = fig.add_subplot(122, projection='3d')

                ax.set_title("Solucao Calculada")
                ax2.set_title("Solucao Analitica")

                plot_surface(X, Y, Z, ax=ax, fig=fig,
                             zmin=zmin, zmax=zmax + zmax / 1000)

                for i in xrange(len(X)):
                    Zsol[i] = sol(X[i], Y[i], t)

                plot_surface(X, Y, Zsol, ax=ax2, fig=fig,
                             zmin=zmin, zmax=zmax + zmax / 1000)

            else:
                ax = fig.add_subplot(111)
                plot_map(X, Y, Z, ax=ax, fig=fig,
                         zmin=zmin, zmax=zmax + zmax / 1000, nx=nx, ny=ny)
                plot_elements(elements, ax=ax)

            plt.suptitle("t = %6.4f" % (step * dt))

            plt.savefig('step_%s.png' % str(step).zfill(4))
            plt.close()

    if gif:
        import imageio
        images = []

        filenames = glob.glob("step_*.png")
        filenames.sort()
        for filename in filenames:
            images.append(imageio.imread(filename))

        imageio.mimsave('movie.gif', images, fps=2)

    return residues


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
    parser.add_argument('--alpha', type=float, default=0.5,
                        help='alpha para método de solução')
    parser.add_argument('--dt', type=float, default=0.002,
                        help='dt para método de solução')
    parser.add_argument('--nsteps', type=int, default=40,
                        help='numero de passos de tempo')

    parser.add_argument('--triangles', type=float, default=0.0,
                        help='porcentagem aproximada de triangulos da malha')

    parser.add_argument('--plotdelta', type=int, default=1,
                        help='graficos que devem ser plotados')

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
    plotdelta = args.plotdelta
    triangles = args.triangles

    run_transient_case(entrada, nx, ny, dt, alpha, newton, triangles,
                       nsteps, plot3D=plot3D, plotdelta=plotdelta, verbose=True)


if __name__ == '__main__':
    main()
