import numpy as np


class FemCase(object):
    """docstring for FemCase"""
    def __init__(self, elements, nodes, ndim=2):
        super(FemCase, self).__init__()
        self.elements = elements
        self.nodes = nodes

        self.P = None

    def SetDirichletBoundary(func=None):
        pass



def BuildStiffnessStrain(elements, neq):

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


def CalcFStrain(elements, neq):
    F = np.zeros(neq)

    nelem = len(elements)

    for iel in xrange(nelem):
        elem = elements[iel]
        Flocal = elem.CalcFlocalStrain()

        for i, node in enumerate(elem.nodes):
            if node.eq is not None:
                F[node.eq] += Flocal[2*i]
                F[node.eq+1] += Flocal[2*i+1]

    return F



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
