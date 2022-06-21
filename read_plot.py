import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

data = np.fromfile('c000000001.out',dtype=np.float32,count=-1)
nx = 100; ny = 100; nz = 100
data = np.reshape(data,(nx,ny,nz))
x = np.linspace(0,0.5,100);y = np.linspace(0,1,100)
fig, ax = plt.subplots()
ax.pcolormesh(x, y, data[:,50,:])
ax.set_aspect('equal', adjustable='datalim')
plt.show()