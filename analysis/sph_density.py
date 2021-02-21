import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("../build/density.txt")
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(data[:,0],data[:,1])


plt.show()
