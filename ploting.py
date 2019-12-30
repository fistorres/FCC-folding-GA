from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

ax = plt.axes(projection='3d')

zline, xline, yline = [], [], []

test = [[1, 1, 0], [1, -1, 0], [-1, -1, 0], [-1, 1, 0],
                     # up
                     [1, 0, 1], [0, -1, 1], [-1, 0, 1], [0, 1, 1],
                     # down
                     [1, 0, -1], [0, -1, -1], [-1, 0, -1], [0, 1,-1]]
points = [[0,0,0]]

for vect in test:
        points.append(list(np.add(np.array(points[0]),np.array(vect))))

for point in points:
    xline.append(point[0])
    yline.append(point[1])
    zline.append(point[2])


#lines
ax.plot3D(xline, yline, zline)
# points
ax.scatter3D(xline, yline, zline)

# ax.set_axis_off()

plt.show()