#!/usr/bin/env python
# -*- coding: utf-8 -*-

from heat import ConstructCase
from FiniteElement import *
from PlotMap import *
import argparse
import glob


def SolveSedo(M, K, F, d0, nsteps, alpha, dt):
    v0 = np.linalg.solve(M, F - np.dot(K, d0))

    print "v0 ", v0.shape

    sols = []

    for step in range(1, nsteps + 1):
        dpreditor = d0 + dt * (1 - alpha) * v0
        vant = v0

        v0 = np.linalg.solve(M + alpha * dt * K, F - np.dot(K, dpreditor))
        d0 = d0 + dt * (1 - alpha) * vant + alpha * v0 * dt

        sols.append(d0.copy())

    return sols


def SolveSedoNewtonImplicit(M, K, F, d0, nsteps, dt):
    sols = []

    K = K * dt
    F = F * dt

    dcurrent = d0
    for step in range(1, nsteps + 1):
        dprevious = dcurrent

        dcurrent = np.linalg.solve(M + K, F + np.dot(M, dprevious))

        sols.append(dcurrent.copy())

    return sols


def main():

    parser = argparse.ArgumentParser(description='Transferencia de Calor 2D')
    parser.add_argument('--entrada', type=int, default=5,
                        help='quantidade de elementos na direção x')
    parser.add_argument('--nx', type=int, default=40,
                        help='quantidade de elementos na direção x')
    parser.add_argument('--ny', type=int, default=40,
                        help='quantidade de elementos na direção y')
    parser.add_argument('--dx', type=float, default=0.1,
                        help='tamanho da célula na direção x')
    parser.add_argument('--dy', type=float, default=0.1,
                        help='tamanho da célula na direção y')
    parser.add_argument('--alpha', type=float, default=0.5,
                        help='alpha para método de solução')
    parser.add_argument('--dt', type=float, default=0.05,
                        help='dt para método de solução')
    parser.add_argument('--nsteps', type=float, default=20,
                        help='numero de passos de tempo')

    parser.add_argument('-n', '--newton', action="store_true")
    parser.add_argument('-p', '--plot', action="store_true")
    args = parser.parse_args()

    nx = args.nx
    ny = args.ny
    entrada = args.entrada
    newton = args.newton
    alpha = args.alpha
    dt = args.dt
    nsteps = args.nsteps

    elements, nodes, neq, sol = ConstructCase(entrada, nx, ny, verbose=False)

    nelem = len(elements)

    def initialCondition(x, y):
        return 0.0

    # ***************************************************************
    #                   Calculo do Lado Direito
    # ***************************************************************

    F = np.zeros(neq)

    for iel in range(nelem):
        elem = elements[iel]
        Flocal = elem.CalcFlocal()

        for i, node in enumerate(elem.nodes):
            if node.eq is not None:
                F[node.eq] += Flocal[i]

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

    if newton:
        sols.extend(SolveSedoNewtonImplicit(M, K, F, d0, nsteps, dt))
    else:
        sols.extend(SolveSedo(M, K, F, d0, nsteps, alpha, dt))

    zmax = np.amax(sols[0])
    zmin = np.amin(sols[0])

    for sol in sols:
        aux = np.amax(sol)
        zmax = aux if aux > zmax else zmax
        zmin = aux if aux < zmin else zmin

    for step in range(nsteps + 1):
        fig = plt.figure(figsize=(16, 9))
        ax = fig.add_subplot(111)

        plot_map(X, Y, sols[step], ax=ax, fig=fig,
                 zmin=zmin, zmax=zmax + zmax / 1000)
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
