#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
from math import *
import pylab as plt
import mpl_toolkits.mplot3d.axes3d as p3
from Element import *

w = sqrt(3) / 3
intPoints = [[-w, -w, 1.0], [w, -w, 1.0], [w, w, 1.0], [-w, w, 1.0]]


DIRICHLET = 1
NEUMANN = 2

def funcform(e, n, vert):

    if vert == 0:
        return (1.0 - e) * (1.0 - n) / 4.0
    elif vert == 1.0:
        return (1.0 + e) * (1.0 - n) / 4.0
    elif vert == 2:
        return (1.0 + e) * (1.0 + n) / 4.0
    elif vert == 3:
        return (1.0 - e) * (1.0 + n) / 4.0

    return float('nan')


def VecFuncForm(e, n):
    v = np.zeros((4, 1))

    for i in range(4):
        v[i] = funcform(e, n, i)

    return v


def dfuncform(e, n, vert, var):

    if vert == 0:
        if var == 0:
            return -0.25 * (1.0 - n)
        if var == 1:
            return -0.25 * (1.0 - e)

    if vert == 1:
        if var == 0:
            return 0.25 * (1.0 - n)
        if var == 1:
            return -0.25 * (1.0 + e)

    if vert == 2:
        if var == 0:
            return 0.25 * (1.0 + n)
        if var == 1:
            return 0.25 * (1.0 + e)

    if vert == 3:
        if var == 0:
            return -0.25 * (1.0 + n)
        if var == 1:
            return 0.25 * (1.0 - e)

    return float('nan')


def FormsDeriv(e, n):
    r = np.zeros((2, 4))
    for var in range(2):
        for f in range(4):
            r[var, f] = dfuncform(e, n, f, var)

    return r


def CoordM(iel, nodesCoord, lg):
    coord = np.zeros((4, 2))

    for i in range(4):
        coord[i, :] = nodesCoord[lg[iel, i], :]

    return coord


def CalcKlocal(iel, nodesCoord, lg, Q):
    Klocal = np.zeros((4, 4))

    coord = CoordM(iel, nodesCoord, lg)

    for e, n, w in intPoints:

        D = FormsDeriv(e, n)

        J = np.dot(D, coord)

        detJ = np.linalg.det(J)

        invJ = np.linalg.inv(J)

        B = np.dot(invJ, D)

        Bt = np.transpose(B)

        Klocal = Klocal + w*np.dot(np.dot(Bt, Q), B) * detJ

    return Klocal


def CalcFlocal(iel, nodesCoord, lg, f, Q, DirichletValues):
    n1, n2, n3, n4 = lg[iel, :]

    F = np.zeros(4)
    P = np.zeros(4)
    fValues = np.zeros(4)

    fValues[0], fValues[1], fValues[2], fValues[3] = f[n1], f[n2], f[n3], f[n4]

    fValues = np.transpose(fValues)
    coord = CoordM(iel, nodesCoord, lg)


    P[0], P[1], P[2], P[3] = DirichletValues[n1], DirichletValues[n2], DirichletValues[n3], DirichletValues[n4]

    for e, n, w in intPoints:
        phi = VecFuncForm(e, n)

        D = FormsDeriv(e, n)
        J = np.dot(D, coord)
        detJ = np.linalg.det(J)

        F = F + np.dot(phi, np.dot(fValues, phi)) * detJ


    # Parcela da condicao de contorno de Neumann
    F = F - np.dot(CalcKlocal(iel, nodesCoord, lg, Q), P)

    return F


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

    Qlist = []

    nelem = nx * ny
    nnodes = (nx + 1) * (ny + 1)

    nodesCoord = np.zeros((nnodes, 2))

    eq = np.zeros(nnodes, dtype=int) - 1

    lg = np.zeros((nelem, 4), dtype=int)

    f = np.zeros(nnodes)
    sol = []


    NeumannValues = np.zeros(nnodes)
    DirichletValues = np.zeros(nnodes)

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

        def solfuncDx(x, y):
            return pi*cos(pi * x) * sin(pi * y)

        def solfuncDy(x, y):
            return sin(pi * x) * pi * cos(pi * y)

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

            nodesCoord[inode, 0], nodesCoord[inode, 1] = x, y

            nodes.append(Node(inode, x, y))

            if getBoundary(i, j, nx, ny) is None:
                eq[inode] = eqCurrent
                nodes[-1].eq = eqCurrent
                sol.append(solfunc(x, y))
                eqCurrent += 1
            else:
                c = contorno[getBoundary(i, j, nx, ny)]
                if c == NEUMANN:
                    eq[inode] = eqCurrent
                    nodes[-1].eq = eqCurrent
                    sol.append(solfunc(x, y))
                    eqCurrent += 1
                else:
                    nodes[inode].p = solfunc(x, y)

            f[inode] = ffunc(x, y)
            nodes[inode].f = ffunc(x, y)

            x += dx
            inode += 1

        y += dy

    sol = np.array(sol)


# ***************************************************************
#                 Vetor da condicao de contorno
# ***************************************************************


# ***************************************************************
#                     Criando os elementos
# ***************************************************************
    iel = 0
    elements = []

    for j in range(ny):
        for i in range(nx):

            firstnode = i + j * (nx + 1)
            lg[iel, 0] = firstnode
            lg[iel, 1] = firstnode + 1
            lg[iel, 2] = firstnode + nx + 2
            lg[iel, 3] = firstnode + nx + 1

            elements.append(Quadrilateral(iel))
            elements[-1].AddNode(nodes[firstnode])
            elements[-1].AddNode(nodes[firstnode + 1])
            elements[-1].AddNode(nodes[firstnode + nx + 2])
            elements[-1].AddNode(nodes[firstnode + nx + 1])

            # adicionando Q como identidade inicialmente
            Qlist.append(np.identity(2))

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
#                Construindo Matriz de Rigidez
# ***************************************************************

    K = np.zeros((neq, neq))

    for iel in range(nelem):
        elem = elements[iel]
        Klocal = elem.CalcKlocal()

        for j in range(4):
            for i in range(4):
                globalI = elem.nodes[i].eq
                globalJ = elem.nodes[j].eq

                if globalI >= 0 and globalJ >= 0:
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

        for i in range(4):
            pos = eq[lg[iel, i]]
            if(pos >= 0):
                F[pos] += Flocal[i]

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
        for inode in range(nnodes):
            pos = eq[inode]
            if pos >= 0:
                x.append(nodesCoord[inode, 0])
                y.append(nodesCoord[inode, 1])
                z.append(calcsol[pos])
                z2.append(sol[pos])

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


if __name__ == '__main__':
    main()

