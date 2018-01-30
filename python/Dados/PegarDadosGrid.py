import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import sys


img=mpimg.imread('GridHomogeneo.png')


ax = plt.gca()
fig = plt.gcf()
implot = ax.imshow(img)

coords = []

def onclick(event):
    global coords
    if event.xdata != None and event.ydata != None:
        if event.button == 1:
            print(event.xdata, event.ydata)
            coords.append((event.xdata, event.ydata))
        if event.button == 3:
            if len(coords) > 0:
                print "Deletado ", coords[-1]
                coords = coords[:-1]



        # print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
        #           ('double' if event.dblclick else 'single', event.button,
        #            event.x, event.y, event.xdata, event.ydata))


def press(event):
    global coords
    print('press', event.key)
    sys.stdout.flush()
    if event.key == 'x':
        print "Salvando arquivo"
        with open("saida.txt", "w") as f:
            for i, c in enumerate(coords):
                f.write("%d %f %f\n" % (i, c[0], c[1]))


cid = fig.canvas.mpl_connect('button_press_event', onclick)
cid = fig.canvas.mpl_connect('key_press_event', press)

plt.show()