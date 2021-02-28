#!/usr/bin/env python
import sys

import matplotlib.pyplot as plt
import numpy as np

filename = sys.argv[1]
output = sys.argv[2]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

n = 100

data = np.loadtxt(filename)

for x, y, z in zip(data[:,1], data[:,2], data[:,3]):
    ax.scatter(x, y, z)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
ax.set_zlim(-1,1)
#plt.show()

fig.savefig(output)