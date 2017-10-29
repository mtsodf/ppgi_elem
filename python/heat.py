import numpy as np

def funcform(e, n, vert):

    if vert == 0:
        return (1.0-e)*(1.0-n)/4.0
    elif vert == 1.0:
        return (1.0+e)*(1.0-n)/4.0
    elif vert == 2:
        return (1.0+e)*(1.0+n)/4.0
    elif vert == 3:
        return (1.0-e)*(1.0+n)/4.0

    return float('nan')

def dfuncform(e, n, vert, var):

    if vert == 0:
        if var == 0:
            return -0.25*(1.0-n)
        if var == 1:
            return -0.25*(1.0-e)

    if vert == 1:
        if var == 0:
            return 0.25*(1.0-n)
        if var == 1:
            return -0.25*(1.0+e)

    if vert == 2:
        if var == 0:
            return 0.25*(1.0+n)
        if var == 1:
            return 0.25*(1.0+e)

    if vert == 3:
        if var == 0:
            return -0.25*(1.0+n)
        if var == 1:
            return  0.25(1.0-e)

    return float('nan')



def formsderivs(e, n):
    r = np.zeros((2,4))
    for var in range(2):
        for f in range(4):
            r[var, f] = dfuncform(e, n, vert, var)

