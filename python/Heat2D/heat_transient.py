#!/usr/bin/env python
# -*- coding: utf-8 -*-

from heat import ConstructCase
from FiniteElement import *
from PlotMap import *
import argparse


def SolveSedo(M, K, F, d0, nsteps, alpha, deltaT):
    v0 = np.linalg.solve(M, F - np.dot(K, d0))

    print "v0 ", v0.shape

    sols = []

    for step in range(1, nsteps + 1):
        dpreditor = d0 + deltaT * (1 - alpha) * v0
        vant = v0

        v0 = np.linalg.solve(M + alpha * deltaT * K, F - np.dot(K, dpreditor))
        d0 = d0 + deltaT * (1 - alpha) * vant + alpha * v0 * deltaT


        sols.append(d0.copy())

    return sols

def SolveSedoNewtonImplicit(M, K, F, d0, nsteps, deltaT):
    sols = []

    K = K*deltaT
    F = F*deltaT

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
    deltaT = 10.0
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
    print "Tamanho K -> ", K.shape

    M = BuildM(elements, neq)
    print "Tamanho M -> ", M.shape

    d0 = np.zeros(neq)

    for node in nodes:
        if node.eq is not None:
            d0[node.eq] = initialCondition(node.coords[0], node.coords[1])

    X = np.array([node.coords[0] for node in nodes if node.eq is not None])
    Y = np.array([node.coords[1] for node in nodes if node.eq is not None])

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    plot_map(X, Y, d0, ax=ax, fig=fig, zmin=0, zmax=100)
    plt.savefig('step_0.png')
    plt.close()

    #sols = SolveSedo(M, K, F, d0, nsteps, alpha, deltaT)

    sols = SolveSedoNewtonImplicit(M, K, F, d0, nsteps, deltaT)


    for step in range(1, nsteps+1):
        fig = plt.figure(figsize=(16, 9))
        ax = fig.add_subplot(111)


        plot_map(X, Y, sols[step-1], ax=ax, fig=fig, zmin=0, zmax=100)

        plt.savefig('step_%d.png' % step)    
        plt.close()




if __name__ == '__main__':
    main()
