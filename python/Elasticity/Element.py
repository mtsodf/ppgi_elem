import numpy as np
from math import sqrt

w = sqrt(3) / 3
intPoints = [[-w, -w, 1.0], [w, -w, 1.0], [w, w, 1.0], [-w, w, 1.0]]


class Element(object):
    """docstring for Element"""

    def __init__(self, iel, E=None, v=None, ndim=1):
        super(Element, self).__init__()

        self.nodes = []
        self.iel = iel
        self.neumannBoundary = []

        self.ndim = ndim

        # Propriedades de Elasticidade
        self.E = E
        self.v = v
        self.qtdNodes = None


    def Hooke(self):

        if self.ndim == 2:

            D = np.zeros((3,3))

            D[0,0] = 1
            D[1,1] = 1
            D[2,2] = (1-self.v)/2

            D[0,1] = self.v
            D[1,0] = self.v

            D = D*(self.E/(1-self.v*self.v))

            return D

        raise NotImplementedError("Metodo nao implementado para ndim = %d" % self.ndim)

    def BeeMat(self, dFuncFormCoordOrig):

        B = np.zeros((3, self.ndim*self.QtdNodes()))

        for i in xrange(self.QtdNodes()):

            B[0, 2*i] = dFuncFormCoordOrig[0, i]
            B[1, 2*i+1] = dFuncFormCoordOrig[1, i]

            B[2, 2*i] = dFuncFormCoordOrig[1, i]
            B[2, 2*i+1] = dFuncFormCoordOrig[0, i]

        return B


    def CalcKlocal(self):

        Klocal = np.zeros((self.QtdDegFree(), self.QtdDegFree()))

        coord = self.GetCoords()

        for e, n, w in intPoints:

            D = self.FormsDeriv(e, n)

            J = np.dot(D, coord)

            detJ = np.linalg.det(J)

            invJ = np.linalg.inv(J)

            B = np.dot(invJ, D)

            B = self.BeeMat(B)

            Bt = np.transpose(B)

            H = self.Hooke()

            Klocal = Klocal + w * np.dot(np.dot(Bt, H), B) * abs(detJ)

        return Klocal


    def GetCoords(self):
        coord = np.zeros((self.QtdNodes(), 2))

        for i in range(self.QtdNodes()):
            coord[i, :] = self.nodes[i].coords

        return coord

    def AddNode(self, node):
        self.nodes.append(node)


    def CalcFlocalStrain(self):

        F = np.zeros(self.QtdNodes()*self.ndim)
        P = np.zeros(self.QtdNodes()*self.ndim)
        fValues = np.zeros((self.QtdNodes(), self.ndim))

        for i, node in enumerate(self.nodes):
            fs = node.f(node.coords[0], node.coords[1])
            fValues[i, 0] = fs[0]
            fValues[i, 1] = fs[1]


        coord = self.GetCoords()
        B = np.zeros((self.QtdNodes()*self.ndim, self.ndim))
        for e, n, w in intPoints:
            phi = self.VecFuncForm(e, n)

            for i in xrange(self.QtdNodes()):
                B[2*i, 0] = phi[i]
                B[2*i+1, 1] = phi[i]

            D = self.FormsDeriv(e, n)
            J = np.dot(D, coord)
            detJ = np.linalg.det(J)
            aux = np.transpose(fValues)
            aux = np.dot(aux, phi)

            for i in xrange(self.QtdNodes()*self.ndim):

                F[i] = F[i] + np.asarray(w*abs(detJ)*np.dot(B, aux))[i]

        return F


    def QtdNodes(self):
        return len(self.nodes)

    def QtdDegFree(self):
        return self.QtdNodes()*self.ndim


    def GetValue(self, sol):
        val = 0.0

        for node in self.nodes:
            if node.eq is None:
                val += node.p
            else:
                val += sol[node.eq]

        return val / self.QtdNodes()


    def VecFuncForm(self, e, n):
        v = np.zeros((self.QtdNodes(), 1))

        for i in range(self.QtdNodes()):
            v[i] = self.FuncForm(e, n, i)

        return v

    def FormsDeriv(self, e, n):
        r = np.zeros((2, self.QtdNodes()))
        for var in range(2):
            for f in range(self.QtdNodes()):
                r[var, f] = self.dFuncForm(e, n, f, var)
        return r


class Quadrilateral(Element):

    def __init__(self, iel, E=None, v=None):
        Element.__init__(self, iel, E, v, 2)
        self.qtdNodes = 4

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


    def dFuncForm(self, e, n, vert, var):

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


    def FuncForm1d(self, e, vert):
        if vert == 0:
            return (1.0 - e) / 2
        if vert == 1:
            return (1.0 + e) / 2

        return None


    def VecFuncForm1d(self, e):
        return np.array([self.FuncForm1d(e, 0), self.FuncForm1d(e, 1)])


    def FuncForm(self, e, n, vert):

        if vert == 0:
            return (1.0 - e) * (1.0 - n) / 4.0
        elif vert == 1:
            return (1.0 + e) * (1.0 - n) / 4.0
        elif vert == 2:
            return (1.0 + e) * (1.0 + n) / 4.0
        elif vert == 3:
            return (1.0 - e) * (1.0 + n) / 4.0

        return float('nan')


    def ToTriangles(self):

        tri0 = Triangle(None, self.Q, self.rho, self.c)

        tri0.AddNode(self.nodes[0])
        tri0.AddNode(self.nodes[1])
        tri0.AddNode(self.nodes[2])

        tri1 = Triangle(None, self.Q, self.rho, self.c)

        tri1.AddNode(self.nodes[0])
        tri1.AddNode(self.nodes[2])
        tri1.AddNode(self.nodes[3])

        return tri0, tri1


class Triangle(Element):

    def __init__(self, iel, E=None, v=None):
        Element.__init__(self, iel, Q, rho, c, E, v, elasticidade)
        self.qtdNodes = 3
        self.quadrilateral = Quadrilateral(None, Q, rho, c)


    def AddNode(self, node):

        Element.AddNode(self, node)

        self.quadrilateral.AddNode(node)

        if self.QtdNodes() == 3:
            self.quadrilateral.AddNode(node)


    def FuncForm(self, e, n, vert):

        if vert == 0 or vert == 1:
            return self.quadrilateral.FuncForm(e, n, vert)

        elif vert == 2:
            return self.quadrilateral.FuncForm(e, n, 2) + self.quadrilateral.FuncForm(e, n, 3)

        return float('nan')


    def dFuncForm(self, e, n, vert, var):

        if vert == 0 or vert == 1:

            return self.quadrilateral.dFuncForm(e, n, vert, var)

        if vert == 2:
            return self.quadrilateral.dFuncForm(e, n, 2, var) + self.quadrilateral.dFuncForm(e, n, 3, var)

        return float('nan')


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
        self.dirichletBoundary = False
