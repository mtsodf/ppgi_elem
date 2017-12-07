#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
from math import *
import pylab as plt
import mpl_toolkits.mplot3d.axes3d as p3
from Element import *
from FiniteElement import *
from PlotMap import *


DIRICHLET = 1
NEUMANN = 2


def getBoundary(i, j, nx, ny):
    if j == 0:
        return 0
    if i == nx:
        return 1
    if j == ny:
        return 2
    if i == 0:
        return 3

    return None


def ConstructCase(entrada, nx, ny, verbose=False):

    rho = None
    c = None

    ffunc = None
    solfunc = None

    sol = []

    if entrada == 0:
        def ffunc(x, y, t=0.0):
            return 2 * pi * pi * sin(pi * x) * sin(pi * y)

        def solfunc(x, y):
            return sin(pi * x) * sin(pi * y)

        lx = 1.0
        ly = 1.0

        x = np.linspace(0.0, lx, num=nx, endpoint=True)
        y = np.linspace(0.0, ly, num=ny, endpoint=True)
        dx = lx / nx
        dy = ly / ny

        rho = 1.0
        c = 1.0

        contorno = [DIRICHLET, DIRICHLET, DIRICHLET, DIRICHLET]

    elif entrada == 1:
        def ffunc(x, y, t=0.0):
            return 2 * pi * pi * sin(pi * x) * sin(pi * y)

        def solfunc(x, y):
            return sin(pi * x) * sin(pi * y)

        def solfuncDx(x):
            return pi * cos(pi * x[0]) * sin(pi * x[1])

        def solfuncDy(x):
            return sin(pi * x[0]) * pi * cos(pi * x[1])

        lx = 0.5
        ly = 1.0

        x = np.linspace(0.0, lx, num=nx, endpoint=True)
        y = np.linspace(0.0, ly, num=ny, endpoint=True)
        dx = lx / nx
        dy = ly / ny

        rho = 1.0
        c = 1.0

        contorno = [DIRICHLET, NEUMANN, DIRICHLET, DIRICHLET]

    elif entrada == 2:
        def ffunc(x, y, t=0.0):
            return 2 * pi * pi * sin(pi * x) * sin(pi * y)

        def solfunc(x, y):
            return sin(pi * x) * sin(pi * y)

        lx = 0.5
        ly = 1.0

        x = np.linspace(0.0, lx, num=nx, endpoint=True)
        y = np.linspace(0.0, ly, num=ny, endpoint=True)
        dx = lx / nx
        dy = ly / ny

        rho = 1.0
        c = 1.0

        contorno = [DIRICHLET, DIRICHLET, DIRICHLET, DIRICHLET]

    elif entrada == 3:
        def ffunc(x, y, t=0.0):
            return 2 * pi * pi * sin(pi * x) * sin(pi * y)

        def solfunc(x, y):
            return sin(pi * x) * sin(pi * y)

        def solfuncDx(x):
            return pi * cos(pi * x[0]) * sin(pi * x[1])

        def solfuncDy(x):
            return sin(pi * x[0]) * pi * cos(pi * x[1])

        lx = 1.0
        ly = 1.0

        x = np.linspace(0.0, lx, num=nx, endpoint=True)
        y = np.linspace(0.0, ly, num=ny, endpoint=True)
        dx = lx / nx
        dy = ly / ny

        rho = 1.0
        c = 1.0

        contorno = [DIRICHLET, NEUMANN, DIRICHLET, DIRICHLET]

    elif entrada == 4:
        def ffunc(x, y, t=0.0):
            return 2 * pi * pi * sin(pi * x) * sin(pi * y)

        def solfunc(x, y):
            return sin(pi * x) * sin(pi * y)

        def solfuncDx(x):
            return pi * cos(pi * x[0]) * sin(pi * x[1])

        def solfuncDy(x):
            return sin(pi * x[0]) * pi * cos(pi * x[1])

        lx = 1.0
        ly = 1.0

        x = np.linspace(0.0, lx, num=nx, endpoint=True)
        y = np.linspace(0.0, ly, num=ny, endpoint=True)
        dx = lx / nx
        dy = ly / ny

        rho = 1.0
        c = 1.0

        contorno = [NEUMANN, DIRICHLET, NEUMANN, NEUMANN]

    elif entrada == 5:
        lx = 1.0
        ly = 1.0

        def ffunc(x, y, t=0.0):
            return 0.0

        def solfunc(x, y):
            return 1.0 if y >= 1.0 or x <= 0.0 else 0.0

        x = np.linspace(0.0, lx, num=nx, endpoint=True)
        y = np.linspace(0.0, ly, num=ny, endpoint=True)

        dx = lx / nx
        dy = ly / ny

        rho = 1.0
        c = 1.0

        contorno = [DIRICHLET, DIRICHLET, DIRICHLET, DIRICHLET]

    if entrada == 6:
        nodes = []

        nodes.append(Node(0, 0.0, 0.0))
        nodes.append(Node(1, 0.5, 0.0))
        nodes.append(Node(2, 1.0, 0.0))
        nodes.append(Node(3, 0.0, 0.5))
        nodes.append(Node(4, 0.5, 0.5))
        nodes.append(Node(5, 1.0, 0.5))

        for node in nodes:
            node.f = 0.0

            nodes[1].eq = 0
            nodes[4].eq = 1

            nodes[2].p = 100.0
            nodes[5].p = 100.0

            nodes[0].p = 0.0
            nodes[3].p = 0.0

            elems = []

            elem0 = Quadrilateral(0)
            elem0.AddNode(nodes[0])
            elem0.AddNode(nodes[1])
            elem0.AddNode(nodes[4])
            elem0.AddNode(nodes[3])

            elems.append(elem0)

            elem1 = Quadrilateral(1)

            elem1.AddNode(nodes[1])
            elem1.AddNode(nodes[2])
            elem1.AddNode(nodes[5])
            elem1.AddNode(nodes[4])

            elems.append(elem1)

            neq = 2

            sol = np.array([50.0, 50.0])

        return elems, nodes, neq, sol

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

            nodes[inode].f = ffunc

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

            elements.append(Quadrilateral(iel, rho=rho, c=c))
            elements[-1].AddNode(nodes[firstnode])
            elements[-1].AddNode(nodes[firstnode + 1])
            elements[-1].AddNode(nodes[firstnode + nx + 2])
            elements[-1].AddNode(nodes[firstnode + nx + 1])

            iel += 1

    neq = eqCurrent

    # ***************************************************************
    #                Setando condicao de Neummann
    # ***************************************************************

    if contorno[0] == NEUMANN:
        for i in range(nx):
            elem = elements[i + nx * 0]
            node1, node2 = elem.getBoundary(0)

            n = np.array([0, -1])

            du = np.array([solfuncDx(node1.coords), solfuncDy(node1.coords)])
            qvetor = - np.dot(elem.Q, du)
            q1 = np.dot(qvetor, n)

            du = np.array([solfuncDx(node2.coords), solfuncDy(node2.coords)])
            qvetor = - np.dot(elem.Q, du)
            q2 = np.dot(qvetor, n)

            elem.setBoundary(0, q1, q2)

    if contorno[1] == NEUMANN:
        for j in range(ny):
            elem = elements[nx - 1 + nx * j]
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
            elem = elements[i + nx * (ny - 1)]

            node1, node2 = elem.getBoundary(2)

            n = np.array([0, 1])

            du = np.array([solfuncDx(node1.coords), solfuncDy(node1.coords)])
            qvetor = - np.dot(elem.Q, du)
            q1 = np.dot(qvetor, n)

            du = np.array([solfuncDx(node2.coords), solfuncDy(node2.coords)])
            qvetor = - np.dot(elem.Q, du)
            q2 = np.dot(qvetor, n)

            elem.setBoundary(2, q1, q2)

    if contorno[3] == NEUMANN:
        for j in range(ny):
            elem = elements[nx * j]

            node1, node2 = elem.getBoundary(3)

            n = np.array([-1, 0])

            du = np.array([solfuncDx(node1.coords), solfuncDy(node1.coords)])
            qvetor = - np.dot(elem.Q, du)
            q1 = np.dot(qvetor, n)

            du = np.array([solfuncDx(node2.coords), solfuncDy(node2.coords)])
            qvetor = - np.dot(elem.Q, du)
            q2 = np.dot(qvetor, n)

            elem.setBoundary(3, q1, q2)

    return elements, nodes, neq, solfunc


def run_case(entrada, nx, ny, verbose=False, plot=False):

    elements, nodes, neq, solfunc = ConstructCase(entrada, nx, ny, verbose)

    sol = np.array([solfunc(node.coords[0], node.coords[1]) for node in nodes if node.eq is not None])

    nelem = len(elements)
    nnodes = len(nodes)

    # ***************************************************************
    #                Construindo Matriz de Rigidez
    # ***************************************************************

    K = BuildStiffness(elements, neq)

    if verbose:
        print "Matriz de Rigidez"
        print K

    # ***************************************************************
    #                   Calculo do Lado Direito
    # ***************************************************************

    F = CalcF(elements, neq)

    if verbose:
        print "Vetor de Carga"
        print F

    # ***************************************************************
    #                 Solucao do Sistema Linear
    # ***************************************************************
    calcsol = np.linalg.solve(K, F)

    # ***************************************************************
    #                 Diferenca entre solucoes
    # ***************************************************************
    if verbose:
        print "Solucao Calculada"
        print calcsol
        print "Solucao Analitica"
        print sol
        print "Diferenca entre solucoes: ", np.linalg.norm(calcsol - sol) / np.linalg.norm(sol)

    # ***************************************************************
    #                Plot das Solucoes
    # ***************************************************************
    if plot:
        z2 = []
        for node in nodes:
            z2.append(solfunc(node.coords[0], node.coords[1]))

        z2 = np.array(z2)

        fig = plt.figure(figsize=(30, 15))
        ax = fig.add_subplot(121, projection='3d')

        X, Y, Z = GetGridValues(nodes, calcsol)

        ax.scatter3D(X, Y, Z)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title("Solucao Calculada")

        ax = fig.add_subplot(122, projection='3d')

        ax.scatter3D(X, Y, np.ravel(z2))
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title("Solucao Analitica")

        plt.show()



        fig = plt.figure(figsize=(16, 9))
        ax = fig.add_subplot(111)
        plot_map(X, Y, Z, ax=ax, fig=fig, zmin=np.amin(z2), zmax=np.amax(z2))
        plt.savefig('permanent.png')
        plt.close()

    return np.linalg.norm(calcsol - sol) / np.linalg.norm(sol)


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
    parser.add_argument('-p', '--plot', action="store_true")
    args = parser.parse_args()

    nx = args.nx
    ny = args.ny
    entrada = args.entrada
    verbose = args.verbose
    plot = args.plot

    residue = run_case(entrada, nx, ny, verbose, plot)

    print "Residuo calculado %f" % residue


if __name__ == '__main__':
    main()
