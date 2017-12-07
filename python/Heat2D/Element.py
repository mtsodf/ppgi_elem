import numpy as np
from math import sqrt

w = sqrt(3) / 3
intPoints = [[-w, -w, 1.0], [w, -w, 1.0], [w, w, 1.0], [-w, w, 1.0]]


def FuncForm1d(e, vert):
    if vert == 0:
        return (1.0 - e) / 2
    if vert == 1:
        return (1.0 + e) / 2

    return None


def VecFuncForm1d(e):
    return np.array([FuncForm1d(e, 0), FuncForm1d(e, 1)])


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

    def __init__(self, iel, Q=None, rho=None, c=None):
        super(Element, self).__init__()

        self.nodes = []
        self.iel = iel
        self.neumannBoundary = []

        if Q is None:
            self.Q = np.identity(2)
        else:
            self.Q = Q

        self.rho = rho
        self.c = c

        self.qtdNodes = None

    def CalcKlocal(self):
        pass

    def GetCoords(self):
        coord = np.zeros((4, 2))

        for i in range(4):
            coord[i, :] = self.nodes[i].coords

        return coord

    def AddNode(self, node):
        self.nodes.append(node)

    def CalcMLocal(self):
        pass

    def QtdNodes(self):
        return len(self.nodes)

    def GetValue(self, sol):
        val = 0.0

        for node in self.nodes:
            if node.eq is None:
                val += node.p
            else:
                val += sol[node.eq]

        return val / self.QtdNodes()


class Quadrilateral(Element):

    def __init__(self, iel, Q=None, rho=None, c=None):
        Element.__init__(self, iel, Q, rho, c)
        self.qtdNodes = 4

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

            Klocal = Klocal + w * np.dot(np.dot(Bt, self.Q), B) * detJ

        return Klocal

    def CalcMLocal(self):
        Mlocal = np.zeros((self.qtdNodes, self.qtdNodes))

        coord = self.GetCoords()

        for e, n, w in intPoints:

            D = FormsDeriv(e, n)

            J = np.dot(D, coord)

            detJ = np.linalg.det(J)

            v = VecFuncForm(e, n)
            Mlocal += w * (np.outer(v, v)) * detJ

        Mlocal = Mlocal * self.rho * self.c
        return Mlocal

    def CalcFlocal(self, t=0.0):

        F = np.zeros(4)
        P = np.zeros(4)
        fValues = np.zeros(4)

        for i, node in enumerate(self.nodes):
            fValues[i] = node.f(node.coords[0], node.coords[1], t)
            P[i] = node.p

        fValues = np.transpose(fValues)
        coord = self.GetCoords()

        for e, n, w in intPoints:
            phi = VecFuncForm(e, n)

            D = FormsDeriv(e, n)
            J = np.dot(D, coord)
            detJ = np.linalg.det(J)

            F = F + np.dot(phi, np.dot(fValues, phi)) * detJ

        # Parcela da condicao de contorno de Dirichlet
        F = F - np.dot(self.CalcKlocal(), P)

        # Parcela da condicao de contorno de Neumann
        for boundary in self.neumannBoundary:

            node1, node2 = self.getBoundary(boundary)
            inode1 = boundary
            inode2 = (boundary + 1) % self.qtdNodes
            Q = np.array([node1.q, node2.q])

            Qf = np.zeros(2)
            for e in [-w, w]:

                detJ = np.linalg.norm(node1.coords - node2.coords)

                phi = VecFuncForm1d(e)

                Qf += np.dot(phi, np.dot(Q, phi)) * detJ / 2

            F[inode1] -= Qf[0]
            F[inode2] -= Qf[1]

        return F

    def setBoundary(self, boundary, q1, q2):

        self.neumannBoundary.append(boundary)

        inode1 = boundary
        inode2 = (boundary + 1) % self.qtdNodes

        self.nodes[inode1].q = q1
        self.nodes[inode2].q = q2

    def getBoundary(self, boundary):
        inode1 = boundary
        inode2 = (boundary + 1) % self.qtdNodes

        return [self.nodes[inode1], self.nodes[inode2]]


class Node(object):
    """docstring for Node"""

    def __init__(self, inode, x, y):
        super(Node, self).__init__()
        self.coords = np.array([x, y])
        self.inode = inode
        self.eq = None
        self.f = None
        self.p = 0
        self.q = 0
