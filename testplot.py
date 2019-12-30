import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import main
figure = plt.figure()
ax = figure.add_subplot(projection="3d")
x = np.array([1,3,3])
y = [2,3,4]
x = np.add(x,y)
z = [1,1,2]

ax.plot(x,y,z)
# ax.set_axis_off()
plt.show()