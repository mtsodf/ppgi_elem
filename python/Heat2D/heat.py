#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
from math import *
import pylab as plt
import mpl_toolkits.mplot3d.axes3d as p3
from Element import *


DIRICHLET = 1
NEUMANN = 2


def getBoundary(i, j, nx, ny):
    if  j == 0:
        return 0
    if i == nx:
        return 1
    if j == ny:
        return 2
    if i == 0:
        return 3

    return None


def run_case(entrada, nx, ny, verbose = False, plot=False):

    nelem = nx * ny
    nnodes = (nx + 1) * (ny + 1)

    sol = []

    if entrada == 0:
        def ffunc(x, y):
            return 2 * pi * pi * sin(pi * x) * sin(pi * y)

        def solfunc(x, y):
            return sin(pi * x) * sin(pi * y)

        lx = 1.0
        ly = 1.0

        x = np.linspace(0.0, lx, num=nx, endpoint=True)
        y = np.linspace(0.0, ly, num=ny, endpoint=True)
        dx = lx / nx
        dy = ly / ny


        contorno = [DIRICHLET, DIRICHLET, DIRICHLET, DIRICHLET]

    if entrada == 1:
        def ffunc(x, y):
            return 2 * pi * pi * sin(pi * x) * sin(pi * y)

        def solfunc(x, y):
            return sin(pi * x) * sin(pi * y)

        def solfuncDx(x):
            return pi*cos(pi * x[0]) * sin(pi * x[1])

        def solfuncDy(x):
            return sin(pi * x[0]) * pi * cos(pi * x[1])

        lx = 0.5
        ly = 1.0

        x = np.linspace(0.0, lx, num=nx, endpoint=True)
        y = np.linspace(0.0, ly, num=ny, endpoint=True)
        dx = lx / nx
        dy = ly / ny

        contorno = [DIRICHLET, NEUMANN, DIRICHLET, DIRICHLET]



    if entrada == 2:
        def ffunc(x, y):
            return 2 * pi * pi * sin(pi * x) * sin(pi * y)

        def solfunc(x, y):
            return sin(pi * x) * sin(pi * y)

        lx = 0.5
        ly = 1.0

        x = np.linspace(0.0, lx, num=nx, endpoint=True)
        y = np.linspace(0.0, ly, num=ny, endpoint=True)
        dx = lx / nx
        dy = ly / ny

        contorno = [DIRICHLET, DIRICHLET, DIRICHLET, DIRICHLET]


    if entrada == 3:
        def ffunc(x, y):
            return 2 * pi * pi * sin(pi * x) * sin(pi * y)

        def solfunc(x, y):
            return sin(pi * x) * sin(pi * y)

        def solfuncDx(x):
            return pi*cos(pi * x[0]) * sin(pi * x[1])

        def solfuncDy(x):
            return sin(pi * x[0]) * pi * cos(pi * x[1])

        lx = 1.0
        ly = 1.0

        x = np.linspace(0.0, lx, num=nx, endpoint=True)
        y = np.linspace(0.0, ly, num=ny, endpoint=True)
        dx = lx / nx
        dy = ly / ny

        contorno = [DIRICHLET, NEUMANN, DIRICHLET, DIRICHLET]

# ***************************************************************
#                        Criando os nós
# ***************************************************************

    nodes = []

    eqCurrent = 0
    y = 0.0
    inode = 0
    for j in range(ny + 1):

        x = 0.0

        for i in range(nx + 1):

            nodes.append(Node(inode, x, y))

            if getBoundary(i, j, nx, ny) is None:
                nodes[-1].eq = eqCurrent
                sol.append(solfunc(x, y))
                eqCurrent += 1
            else:
                c = contorno[getBoundary(i, j, nx, ny)]
                if c == NEUMANN:
                    nodes[-1].eq = eqCurrent
                    sol.append(solfunc(x, y))
                    eqCurrent += 1
                else:
                    nodes[inode].p = solfunc(x, y)

            nodes[inode].f = ffunc(x, y)

            x += dx
            inode += 1

        y += dy

    sol = np.array(sol)

# ***************************************************************
#                     Criando os elementos
# ***************************************************************
    iel = 0
    elements = []

    for j in range(ny):
        for i in range(nx):

            firstnode = i + j * (nx + 1)

            elements.append(Quadrilateral(iel))
            elements[-1].AddNode(nodes[firstnode])
            elements[-1].AddNode(nodes[firstnode + 1])
            elements[-1].AddNode(nodes[firstnode + nx + 2])
            elements[-1].AddNode(nodes[firstnode + nx + 1])


            iel += 1

    neq = eqCurrent

    if verbose:
        print "Coordenadas dos nos"
        print nodesCoord

        print "Eq"
        print eq

        print "Nos dos elementos"
        print lg

        print "Numero de equacoes"
        print neq

# ***************************************************************
#                Setando condicao de Neummann
# ***************************************************************

    if contorno[0] == NEUMANN:
        for i in range(nx):
            elem = elements[i + nx*0]
            node1, node2 = elem.getBoundary(0)

            n = np.array([0,-1])

            du = np.array([solfuncDx(node1.coords), solfuncDy(node1.coords)])
            qvetor = - np.dot(elem.Q, du)
            q1 = np.dot(qvetor, n)

            du = np.array([solfuncDx(node2.coords), solfuncDy(node2.coords)])
            qvetor = - np.dot(elem.Q, du)
            q2 = np.dot(qvetor, n)

            elem.setBoundary(0, q1, q2)

    if contorno[1] == NEUMANN:
        for j in range(ny):
            elem = elements[nx - 1+ nx*j]
            node1, node2 = elem.getBoundary(1)

            n = np.array([1, 0])

            du = np.array([solfuncDx(node1.coords), solfuncDy(node1.coords)])
            qvetor = - np.dot(elem.Q, du)
            q1 = np.dot(qvetor, n)

            du = np.array([solfuncDx(node2.coords), solfuncDy(node2.coords)])
            qvetor = - np.dot(elem.Q, du)
            q2 = np.dot(qvetor, n)


            elem.setBoundary(1, q1, q2)

    if contorno[2] == NEUMANN:
        for i in range(nx):
            elem = elements[i + nx*(ny-1)]
            elem.setBoundary(2, 0.0, 1.0)

    if contorno[3] == NEUMANN:
        for j in range(ny):
            elem = elements[nx*j]
            elem.setBoundary(3, 0.0, 1.0)


# ***************************************************************
#                Construindo Matriz de Rigidez
# ***************************************************************

    K = np.zeros((neq, neq))

    for iel in range(nelem):
        elem = elements[iel]
        Klocal = elem.CalcKlocal()

        for j, nodej in enumerate(elem.nodes):
            for i, nodei in enumerate(elem.nodes):
                globalI = nodei.eq
                globalJ = nodej.eq

                if globalI is not None and globalJ is not None:
                    K[globalI, globalJ] += Klocal[i, j]

    if verbose:
        print K

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
#                 Adicao de Condicao de Neumann
# ***************************************************************


# ***************************************************************
#                 Solucao do Sistema Linear
# ***************************************************************
    calcsol = np.linalg.solve(K, F)


# ***************************************************************
#                 Diferenca entre solucoes
# ***************************************************************
    if verbose:
        print calcsol
        print sol
    print "Diferenca entre solucoes: ", np.linalg.norm(calcsol - sol)/np.linalg.norm(sol)

# ***************************************************************
#                Plot das Solucoes
# ***************************************************************
    if plot:
        x = []
        y = []
        z = []
        z2 = []
        for node in nodes:
            if node.eq is not None:
                x.append(node.coords[0])
                y.append(node.coords[1])
                z.append(calcsol[node.eq])
                z2.append(sol[node.eq])

        x = np.array(x)
        y = np.array(y)
        z = np.array(z)

        fig = plt.figure(figsize=(30, 15))
        ax = fig.add_subplot(121, projection='3d')

        ax.scatter3D(np.ravel(x), np.ravel(y), np.ravel(z))
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title("Solucao Calculada")

        ax = fig.add_subplot(122, projection='3d')
        ax.scatter3D(np.ravel(x), np.ravel(y), np.ravel(z2))
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title("Solucao Analitica")

        plt.show()

    return np.linalg.norm(calcsol - sol)/np.linalg.norm(sol)


def main():

    parser = argparse.ArgumentParser(description='Transferencia de Calor 2D')
    parser.add_argument('--entrada', type=int, default=0,
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
    args = parser.parse_args()

    nx = args.nx
    ny = args.ny
    entrada = args.entrada
    verbose = args.verbose



    residue = run_case(entrada, nx, ny, verbose)


    print "Residuo calculado %f" % residue

if __name__ == '__main__':
    main()

