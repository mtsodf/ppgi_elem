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
