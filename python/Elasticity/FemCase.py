

class FemCase(object):
    """docstring for FemCase"""
    def __init__(self):
        super(FemCase, self).__inidvit__()
        self.nodesCoord = None
        self.lg         = None
        self.eq         = None
        self.nodes      = []

    def addNode(self, node):
        self.nodes.append(node)



class Element(object):
    """docstring for Element"""
    def __init__(self, nodes):
        super(Element, self).__init__()
        self.nodes = nodes


class Node(object):
    """docstring for Node"""
    def __init__(self, case):
        super(Node, self).__init__()
        self.num = None
        self.case = case
        self.free = None

    def getCoords(self):
        return self.case.nodesCoord[self.num, :]

    def getEq(self):
        return self.case.eq[self.num]



def CreateCase(nx = 10, ny = 10, dx = 0.1, dy = 0.1):
    case = Case()
    nelem = nx*ny
    nnodes = (nx+1)*(ny+1)
    case.nodesCoord = np.zeros((nnodes,2))

    y = 0.0
    inode = 0
    for j in range(ny+1):
        x = 0.0
        for i in range(nx+1):
            case.nodesCoord = case.nodesCoord[inode, (x, y)]
            self.addNode(case.node)
            x += dx
            inode += 1

        y += dy
