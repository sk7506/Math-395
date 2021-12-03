#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 15:02:24 2021

@author: skilpatrick
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import animation

x = [0, 1, 2]
y = [0, 1, 2]
yaw = [0.0, 0.5, 1.3]
fig = plt.figure()
plt.axis('equal')
plt.grid()
ax = fig.add_subplot(111)
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)

patch = patches.Rectangle((0, 0), 0, 0, fc='y')

def init():
    ax.add_patch(patch)
    return patch,

def animate(i):
    patch.set_width(1.2)
    patch.set_height(1.0)
    patch.set_xy([x[i], y[i]])
    patch._angle = -np.rad2deg(yaw[i])
    return patch,

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(x), interval=500, blit=True)
plt.show()
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches

fig, ax = plt.subplots()

ax.plot([1,4],[1,4])

ax.add_patch(
     patches.Rectangle(
        (1, 1),
        0.5,
        0.5,
        edgecolor = 'blue',
        facecolor = 'red',
        fill=True
     ) )

plt.show()
"""
"""
import numpy as np
from matplotlib import pyplot as plt
plt.xlim(-30,30)
plt.ylim(-60,5)
plt.vlines(x = 0, ymin = -35, ymax = 0, colors='green', ls='-', lw=2, label='vline_single - full height')
plt.vlines(x = 5, ymin = -55, ymax = -35, colors='green', ls='-', lw=2, label='vline_single - full height')
plt.hlines(y = -35, xmin = 0, xmax = 5, colors='green', linestyles='-', lw=2, label='Single Short Line')
plt.hlines(y = -55, xmin = 2.5, xmax = 5, colors='green', linestyles='-', lw=2, label='Single Short Line')

plt.show()
"""
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D


# create cluster builder
def cluster(center, radius=10, n=50):
    c = np.asarray(center)
    xx, yy, zz = np.random.uniform(c-radius, c+radius, size=(n, 3)).T
    return xx, yy, zz


# create random cluster
xx1, yy1, zz1 = cluster((25, 15, 5))
xx2, yy2, zz2 = cluster((5, 5, 20))

# create a figure
fig = plt.figure()
# initialise 3D Axes
ax = Axes3D(fig)
ax.grid(False)


# create the initialiser with our two clusters
def init():
    ax.scatter(xx1, yy1, zz1, marker='o', s=40, c='#08C6AB', alpha=0.6)
    ax.scatter(xx2, yy2, zz2, marker='o', s=40, c='#37465B', alpha=0.6)
    return fig,


# create animate function, this will adjust the view one step at a time
def animate(i):
    ax.view_init(elev=10.0, azim=i)
    return fig,


# create the animated plot
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=360, interval=20, blit=True)
# save as a GIF
anim.save('clusters.gif', fps=30, writer='imagemagick')
"""

