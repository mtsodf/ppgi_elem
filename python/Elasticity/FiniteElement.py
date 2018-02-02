import numpy as np


class FemCase(object):
    """docstring for FemCase"""
    def __init__(self, elements, nodes, ndim=2):
        super(FemCase, self).__init__()
        self.elements = elements
        self.nodes = nodes


    def SetDirichletBoundary(func=None):
        pass


class FemContext(object):
    """docstring for FemContext"""
    def __init__(self):
        super(FemContext, self).__init__()

        self.eq = {}
        self.P = {}
        self.neq = None

    def getEq(self, node):
        return self.eq[node.inode] if node.inode in self.eq else [None, None]

    def setEqs(self, node, eqs):
        self.eq[node.inode] = eqs

    def setP(self, node, Pvalue):
        self.P[node.inode] = Pvalue


def BuildStiffnessStrain(context, elements):
    neq = context.neq
    K = np.zeros((neq, neq))
    nelem = len(elements)

    for iel in range(nelem):
        elem = elements[iel]
        Klocal = elem.CalcKlocal()

        for j, nodej in enumerate(elem.nodes):
            for i, nodei in enumerate(elem.nodes):
                globalIx, globalIy = context.getEq(nodei)
                globalJx, globalJy = context.getEq(nodej)

                if globalIx is not None and globalJx is not None:
                    K[globalIx, globalJx] += Klocal[2*i, 2*j]

                if globalIy is not None and globalJx is not None:
                    K[globalIy, globalJx] += Klocal[2*i+1, 2*j]

                if globalIx is not None and globalJy is not None:
                    K[globalIx, globalJy] += Klocal[2*i, 2*j+1]

                if globalIx is not None and globalJx is not None:
                    K[globalIy, globalJy] += Klocal[2*i+1, 2*j+1]

    return K


def CalcFStrain(context, elements):

    neq = context.neq
    F = np.zeros(neq)

    nelem = len(elements)

    for iel in xrange(nelem):
        elem = elements[iel]
        Flocal = elem.CalcFlocalStrain()

        for i, node in enumerate(elem.nodes):
            eqx, eqy = context.getEq(node)

            if eqx is not None:
                F[eqx] += Flocal[2*i]
            if eqy is not None:
                F[eqy] += Flocal[2*i+1]

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
