import numpy as np


def BuildStiffness(elements, neq):

    K = np.zeros((neq, neq))
    nelem = len(elements)

    for iel in range(nelem):
        elem = elements[iel]
        Klocal = elem.CalcKlocal()

        for j, nodej in enumerate(elem.nodes):
            for i, nodei in enumerate(elem.nodes):
                globalI = nodei.eq
                globalJ = nodej.eq

                if globalI is not None and globalJ is not None:
                    K[globalI, globalJ] += Klocal[i, j]

    return K

def BuildStiffnessElasticity(elements, neq):

    K = np.zeros((neq, neq))
    nelem = len(elements)

    for iel in range(nelem):
        elem = elements[iel]
        Klocal = elem.CalcKlocal()

        for j, nodej in enumerate(elem.nodes):
            for i, nodei in enumerate(elem.nodes):
                globalI = nodei.eq
                globalJ = nodej.eq

                if globalI is not None and globalJ is not None:
                    for di in range(2):
                        for dj in range(2):
                            K[globalI+di, globalJ+dj] += Klocal[2*i+di, 2*j+dj]


    return K

def CalcFElasticity(elements, neq):
    F = np.zeros(neq)

    nelem = len(elements)

    for iel in xrange(nelem):
        elem = elements[iel]
        Flocal = elem.CalcFlocalElasticity()

        for i, node in enumerate(elem.nodes):
            if node.eq is not None:
                F[node.eq] += Flocal[2*i]
                F[node.eq+1] += Flocal[2*i+1]

    return F



def CalcF(elements, neq, t=0.0):
    F = np.zeros(neq)

    nelem = len(elements)

    for iel in range(nelem):
        elem = elements[iel]
        Flocal = elem.CalcFlocal(t=t)

        for i, node in enumerate(elem.nodes):
            if node.eq is not None:
                F[node.eq] += Flocal[i]
    return F


def BuildM(elements, neq):
    M = np.zeros((neq, neq))
    nelem = len(elements)

    for iel in range(nelem):
        elem = elements[iel]
        Mlocal = elem.CalcMLocal()

        for j, nodej in enumerate(elem.nodes):
            for i, nodei in enumerate(elem.nodes):
                globalI = nodei.eq
                globalJ = nodej.eq

                if globalI is not None and globalJ is not None:
                    M[globalI, globalJ] += Mlocal[i, j]

    return M


def GetGridValues(nodes, sol):
    nnodes = len(nodes)
    val = np.zeros(nnodes)
    x = np.zeros(nnodes)
    y = np.zeros(nnodes)

    for inode, node in enumerate(nodes):
        x[inode], y[inode] = node.coords
        if node.eq is None:
            val[inode] = node.p
        else:
            val[inode] = sol[node.eq]

    return x, y, val
