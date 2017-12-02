#!/usr/bin/env python
# -*- coding: utf-8 -*-

from heat import ConstructCase
from FiniteElement import *
from PlotMap import *
import argparse


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
    parser.add_argument('-v', '--verbose', action="store_true")
    parser.add_argument('-p', '--plot', action="store_true")
    args = parser.parse_args()

    nx = args.nx
    ny = args.ny
    entrada = args.entrada
    verbose = args.verbose
    plot = args.plot

    elements, nodes, neq, sol = ConstructCase(entrada, nx, ny, verbose=False)

    nelem = len(elements)

    def initialCondition(x, y):
        return 0.0

    alpha = 0.5
    deltaT = 0.1

    nsteps = 10

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

    M = BuildM(elements, neq)

    d0 = np.zeros(neq)

    for node in nodes:
        if node.eq is not None:
            d0[node.eq] = initialCondition(node.coords[0], node.coords[1])

    X = np.array([node.coords[0] for node in nodes if node.eq is not None])
    Y = np.array([node.coords[1] for node in nodes if node.eq is not None])


    plot_surface(X, Y, d0, out='step_0.png')

    v0 = np.linalg.solve(M, F - np.dot(K, d0))

    for step in range(1, nsteps + 1):
        dpreditor = d0 + deltaT * (1 - alpha) * v0

        vant = v0

        v0 = np.linalg.solve(M + alpha * deltaT * K, F - K * dpreditor)

        d0 = d0 + deltaT * (1 - alpha) * vant + alpha * v0 * deltaT

        print d0

        plot_surface(X, Y, d0, out='step_%d.png' % step)


if __name__ == '__main__':
    main()
