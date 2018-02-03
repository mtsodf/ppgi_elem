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
import ipdb as pdb

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


def ConstructCase(entrada, nx, ny, triangles=0.0, verbose=False, distorce=False):

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

        context = FemContext()

        contorno = [DIRICHLET, DIRICHLET, DIRICHLET, DIRICHLET]

    if entrada == 1:
        def ffunc(x, y, t=0.0):
            E = 1
            v = 0.25
            f = np.zeros(2)
            f[0] = (-pi*pi*E*v/(v*v-1)) *cos(pi*x)*cos(pi*y) + pi*pi*E/(2*v+2) * cos(pi*x) * cos(pi*y)
            f[1] = pi*pi*(E/(v*v-1))*sin(pi*x)*sin(pi*y) - (pi*pi*E/(2*v+2))*sin(pi*x)*sin(pi*y)
            return f

        def solfunc(x, y):
            return np.array([0, sin(pi*x)*sin(pi*y)])

        lx = 1.0
        ly = 1.0

        x = np.linspace(0.0, lx, num=nx+1, endpoint=True)
        y = np.linspace(0.0, ly, num=ny+1, endpoint=True)

        dx = lx / nx
        dy = ly / ny


        E = 1.0
        v = 0.25

        context = FemContext()

        contorno = [DIRICHLET, DIRICHLET, DIRICHLET, DIRICHLET]


    if entrada == 2:
        def ffunc(x, y, t=0.0):
            E = 1
            v = 0.25
            f = np.zeros(2)
            f[0] = pi*pi*(E/(v*v-1))*sin(pi*x)
            f[1] = 0.0
            return f

        def solfunc(x, y):
            return np.array([sin(pi*x), 0])

        lx = 1.0
        ly = 1.0

        x = np.linspace(0.0, lx, num=nx+1, endpoint=True)
        y = np.linspace(0.0, ly, num=ny+1, endpoint=True)

        dx = lx / nx
        dy = ly / ny


        E = 1.0
        v = 0.25

        context = FemContext()

        contorno = [DIRICHLET, DIRICHLET, DIRICHLET, DIRICHLET]

    if entrada == 3:
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
        ly = 2.0

        x = np.linspace(0.0, lx, num=nx+1, endpoint=True)
        y = np.linspace(0.0, ly, num=ny+1, endpoint=True)

        dx = lx / nx
        dy = ly / ny


        E = 1.0
        v = 0.25

        context = FemContext()

        contorno = [DIRICHLET, DIRICHLET, DIRICHLET, DIRICHLET]

    # ***************************************************************
    #                        Criando os nós
    # ***************************************************************

    nodes = []

    eqCurrent = 0
    y = yOrigin
    inode = 0

    for j in range(ny + 1):

        if j == ny:
            y = yOrigin + ly

        x = xOrigin

        for i in range(nx + 1):

            if i == nx:
                x = xOrigin + lx

            xNode = x
            yNode = y

            if distorce and i!= 0 and j!= 0 and i!=nx and j!=ny:
                xNode += np.random.uniform(-0.1*dx, 0.1*dx, 1)[0]
                yNode += np.random.uniform(-0.1*dy, 0.1*dy, 1)[0]

            nodes.append(Node(inode, xNode, yNode))

            if getBoundary(i, j, nx, ny) is None:
                context.setEqs(nodes[-1], [eqCurrent, eqCurrent + 1])
                sol.append(solfunc(xNode, yNode))
                context.setP(nodes[-1], [0.0, 0.0])
                eqCurrent += 2
            else:
                c = contorno[getBoundary(i, j, nx, ny)]
                if c == NEUMANN:
                    raise NotImplementedError("Condicao de contorno de Neumann nao implementada ainda.")
                else:
                    context.setP(nodes[-1], solfunc(xNode, yNode))

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

            elements.append(Quadrilateral(iel, E=E, v=v))
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

    context.neq = eqCurrent

    context.allElements = elements
    context.allNodes = nodes


    context.setAllElementsAndNodesToContext()

    return context, solfunc


def run_case(entrada, nx, ny, triangles=0.0, verbose=False, plot=False, distorce=False):

    context, solfunc = ConstructCase(entrada, nx, ny, triangles, verbose, distorce)

    sol = np.zeros(context.neq)
    for node in context.nodesIter():
        eqs = context.getEq(node)
        solCalc = solfunc(node.coords[0], node.coords[1])
        if eqs[0] is not None:
            sol[eqs[0]] = solCalc[0]
        if eqs[1] is not None:
            sol[eqs[1]] = solCalc[1]


    # ***************************************************************
    #                Construindo Matriz de Rigidez
    # ***************************************************************

    K = BuildStiffnessElasticity(context)

    if verbose:
        print "Matriz de Rigidez"
        print K


    # ***************************************************************
    #                   Calculo do Lado Direito
    # ***************************************************************

    F = CalcFElasticity(context)

    if verbose:
        print "Vetor de Carga"
        print F



    # ***************************************************************
    #                 Solucao do Sistema Linear
    # ***************************************************************


    #TODO precisa ser revisado porque esse valor negativo é necessário.
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

        fig = plt.figure(figsize=(30,16))
        ax=fig.add_subplot(1,2,1)

        plot_solution(calcsol, context, uxAxis=ax, fig=fig)
        plot_elements(context.elementsIter(), ax)

        ax=fig.add_subplot(1,2,2)
        plot_solution(calcsol, context, uyAxis=ax, fig=fig)
        plot_elements(context.elementsIter(), ax)

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

    print "Residuo encontrado = ", residue

if __name__ == '__main__':
    main()
