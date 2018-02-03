import numpy as np


class FemCase(object):
    """docstring for FemCase"""
    def __init__(self, elements, nodes, ndim=2):
        super(FemCase, self).__init__()
        self.elements = elements
        self.nodes = nodes


    def SetDirichletBoundary(func=None):
        pass

    def NodesIterator(self):
        return self.nodes.__iter__()

    def QtdNodes(self):
        return len(self.nodes)


class IteratorByList:
    def __init__(self, listObjects, listIndexes):
        self.current = 0
        self.listObjects = listObjects
        self.listIndexes = listIndexes

    def __iter__(self):
        return self

    def next(self): # Python 3: def __next__(self)
        if self.current >= len(self.listIndexes):
            raise StopIteration
        else:
            self.current += 1
            return self.listObjects[self.current-1]

class FemContext(object):
    """docstring for FemContext"""
    def __init__(self):
        super(FemContext, self).__init__()

        self.eq = {}
        self.P = {}
        self.neq = None

        self.allElements = None
        self.allNodes = None

        self.elements = None
        self.nodes = None



    def getEq(self, node):
        return self.eq[node.inode] if node.inode in self.eq else [None, None]

    def setEqs(self, node, eqs):
        self.eq[node.inode] = eqs

    def getP(self, node):
        if node.inode in self.P:
            return self.P[node.inode]

        return [0.0, 0.0]

    def setP(self, node, Pvalue):
        self.P[node.inode] = Pvalue


    def setAllElementsAndNodesToContext(self):
        self.elements = np.arange(len(self.allElements))
        self.nodes = np.arange(len(self.allNodes))


    def QtdNodes(self):
        return len(self.nodes)

    def elementsIter(self):
        return IteratorByList(self.allElements, self.elements)

    def nodesIter(self):
        return IteratorByList(self.allNodes, self.nodes)


def BuildStiffnessElasticity(context):
    neq = context.neq
    K = np.zeros((neq, neq))

    for elem in context.elementsIter():
        #elem = elements[iel]
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


def CalcFElasticity(context):

    neq = context.neq
    F = np.zeros(neq)


    for elem in context.elementsIter():

        Flocal = elem.CalcFlocalElasticity(context)

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
