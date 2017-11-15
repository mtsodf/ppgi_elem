import numpy as np
from math import sqrt

w = sqrt(3) / 3
intPoints = [[-w, -w, 1.0], [w, -w, 1.0], [w, w, 1.0], [-w, w, 1.0]]

def FuncForm(e, n, vert):

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
        v[i] = FuncForm(e, n, i)

    return v


def dFuncForm(e, n, vert, var):

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
            r[var, f] = dFuncForm(e, n, f, var)

    return r



class Element(object):
    """docstring for Element"""
    def __init__(self, iel, Q=None):
        super(Element, self).__init__()
        self.nodes = []
        self.iel = iel

        if Q is None:
            self.Q = np.identity(2)
        else:
            self.Q = Q


    def CalcKlocal(self):
        pass

    def GetCoords(self):
        coord = np.zeros((4,2))

        for i in range(4):
            coord[i,:] = self.nodes[i].coords

        return coord

    def AddNode(self, node):
        self.nodes.append(node)



class Quadrilateral(Element):

    def CalcKlocal(self):
        Klocal = np.zeros((4, 4))

        coord = self.GetCoords()

        for e, n, w in intPoints:

            D = FormsDeriv(e, n)

            J = np.dot(D, coord)

            detJ = np.linalg.det(J)

            invJ = np.linalg.inv(J)

            B = np.dot(invJ, D)

            Bt = np.transpose(B)

            Klocal = Klocal + w*np.dot(np.dot(Bt, self.Q), B) * detJ

        return Klocal

        pass

    def CalcFlocal(self):

        F = np.zeros(4)
        P = np.zeros(4)
        fValues = np.zeros(4)

        for i, node in enumerate(self.nodes):
            fValues[i] = node.f
            P[i] = node.p


        fValues = np.transpose(fValues)
        coord = self.GetCoords()



        for e, n, w in intPoints:
            phi = VecFuncForm(e, n)

            D = FormsDeriv(e, n)
            J = np.dot(D, coord)
            detJ = np.linalg.det(J)

            F = F + np.dot(phi, np.dot(fValues, phi)) * detJ


        # Parcela da condicao de contorno de Dirichelet
        F = F - np.dot(self.CalcKlocal(), P)

        return F
class Node(object):
    """docstring for Node"""
    def __init__(self, inode, x, y):
        super(Node, self).__init__()
        self.coords = [x, y]
        self.inode = inode
        self.eq = None
        self.f = None
        self.p = 0
        self.q = 0
