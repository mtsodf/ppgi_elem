import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np

def newline(p1, p2):
    ax = plt.gca()
    xmin, xmax = ax.get_xbound()

    if(p2[0] == p1[0]):
        xmin = xmax = p1[0]
        ymin, ymax = ax.get_ybound()
    else:
        ymax = p1[1]+(p2[1]-p1[1])/(p2[0]-p1[0])*(xmax-p1[0])
        ymin = p1[1]+(p2[1]-p1[1])/(p2[0]-p1[0])*(xmin-p1[0])

    l = mlines.Line2D([xmin,xmax], [ymin,ymax])
    ax.add_line(l)
    return l


xs = []
ys = []
with open("GridGrosso.txt") as f:
    for line in f:
        x, y = line.split()[1:]
        xs.append(float(x))
        ys.append(float(y))

    xs = np.array(xs)
    ys = np.array(ys)

    xOrigin = xs[0]
    yOrigin = ys[0]
    for i in xrange(len(xs)):
        xs[i] = xs[i] - xOrigin
        ys[i] = yOrigin - ys[i]


for i in xrange(len(xs)):
    print xs[i], ys[i]
x = np.linspace(0, np.max(xs))
y = np.linspace(0, np.max(ys))

plt.plot(x, y, color="white")

for i in range(len(xs)-1):
    if i%6 != 5:
        plt.plot((xs[i], xs[i+1]), (ys[i], ys[i+1]), color="blue")

plt.show()