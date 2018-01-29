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
from math import exp
from random import random


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


def ConstructCase(entrada, nx, ny, triangles=0.0, verbose=False):

    rho = None
    c = None

    ffunc = None
    solfunc = None

    sol = []

    initialCondition = None

    xOrigin = 0.0
    yOrigin = 0.0

    if entrada == 0:
        def ffunc(x, y, t=0.0):
            E = 1
            v = 0.25
            f = np.zeros(2)
            f[0] = pi*pi*(E/(v*v-1))*sin(pi*x)*sin(pi*y) - (pi*pi*E/(2*v+2))*sin(pi*x)*sin(pi*y)
            f[1] = (-pi*pi*E*v/(v*v-1)) *cos(pi*x)*cos(pi*y) + pi*pi*E/(2*v+2) * cos(pi*x) * cos(pi*y)
            return f

        def solfunc(x, y):
            return np.array([sin(pi * x)*sin(pi*y), 0])

        lx = 1.0
        ly = 1.0

        x = np.linspace(0.0, lx, num=nx+1, endpoint=True)
        y = np.linspace(0.0, ly, num=ny+1, endpoint=True)

        dx = lx / nx
        dy = ly / ny


        E = 1.0
        v = 0.25

        contorno = [DIRICHLET, DIRICHLET, DIRICHLET, DIRICHLET]


    # ***************************************************************
    #                        Criando os nós
    # ***************************************************************

    nodes = []

    eqCurrent = 0
    y = yOrigin
    inode = 0
    for j in range(ny + 1):

        x = xOrigin

        for i in range(nx + 1):

            nodes.append(Node(inode, x, y))

            if getBoundary(i, j, nx, ny) is None:
                nodes[-1].eq = eqCurrent
                sol.append(solfunc(x, y))
                eqCurrent += 2
            else:
                c = contorno[getBoundary(i, j, nx, ny)]
                if c == NEUMANN:
                    nodes[-1].eq = eqCurrent
                    sol.append(solfunc(x, y))
                    eqCurrent += 1
                else:
                    nodes[-1].dirichletBoundary = True
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

            elements.append(Quadrilateral(iel, E=E, v=v, ndim=2, elasticidade=True))
            elements[-1].AddNode(nodes[firstnode])
            elements[-1].AddNode(nodes[firstnode + 1])
            elements[-1].AddNode(nodes[firstnode + nx + 2])
            elements[-1].AddNode(nodes[firstnode + nx + 1])

            iel += 1

    elements_old = elements

    elements = []

    for elem in elements_old:

        if random() < triangles:
            tri0, tri1 = elem.ToTriangles()

            elements.append(tri0)
            elements.append(tri1)

        else:
            elements.append(elem)

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


def run_case(entrada, nx, ny, triangles=0.0, verbose=False, plot=False):

    elements, nodes, neq, solfunc = ConstructCase(entrada, nx, ny, triangles, verbose)


    sol = np.zeros(neq)
    for node in nodes:
        if node.eq is not None:
            sol[node.eq], sol[node.eq + 1] = solfunc(node.coords[0], node.coords[1])

    nelem = len(elements)
    nnodes = len(nodes)

    # ***************************************************************
    #                Construindo Matriz de Rigidez
    # ***************************************************************

    K = BuildStiffnessStrain(elements, neq)

    if verbose:
        print "Matriz de Rigidez"
        print K


    # ***************************************************************
    #                   Calculo do Lado Direito
    # ***************************************************************

    F = CalcFStrain(elements, neq)

    if verbose:
        print "Vetor de Carga"
        print F



    # ***************************************************************
    #                 Solucao do Sistema Linear
    # ***************************************************************
    calcsol = np.linalg.solve(K, -F)

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
        plot_elements(elements, ax)
        plt.savefig('permanent.png')
        plt.close()

    return np.linalg.norm(calcsol - sol) / np.linalg.norm(sol)


def main():

    parser = argparse.ArgumentParser(description='Transferencia de Calor 2D')
    parser.add_argument('--entrada', type=int, default=0,
                        help='quantidade de elementos na direção x')
    parser.add_argument('--nx', type=int, default=10,
                        help='quantidade de elementos na direção x')
    parser.add_argument('--ny', type=int, default=10,
                        help='quantidade de elementos na direção y')
    parser.add_argument('--triangles', type=float, default=0.0,
                        help='porcentagem aproximada de triangulos da malha')

    parser.add_argument('-v', '--verbose', action="store_true")
    parser.add_argument('-p', '--plot', action="store_true")
    args = parser.parse_args()

    nx = args.nx
    ny = args.ny
    entrada = args.entrada
    triangles = args.triangles
    verbose = args.verbose
    plot = args.plot

    residue = run_case(entrada, nx, ny, triangles, verbose, plot)

if __name__ == '__main__':
    main()
