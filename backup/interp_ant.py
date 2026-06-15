import numpy as np
import ast
import os
import matplotlib.pyplot as plt
#
# import matplotlib
# matplotlib.use("TkAgg")
# from matplotlib import pyplot as plt
# os.environ['PROJ_LIB'] = '/Users/micheltossaint/Documents/anaconda3/share/proj'

from scipy import interpolate

test = '(0.00, 33.18),(0.72, 32.68),(1.39, 31.05),(2.06, 28.12),(2.53, 24.76),(2.99, 21.05),(3.32, 16.65),(3.62, 10.25),\
       (3.77, 5.84),(3.83, 3.51),(3.87, 3.08),(3.92, 3.51),(4.02, 6.28),(4.16, 9.83),(4.42, 13.7),(4.68, 15.9),\
       (5.00, 17.29),(5.52, 17.90),(6.08, 17.05),(6.74, 14.82),(7.28, 13.01),(7.84, 11.47),(9.01, 8.05 ),(10.44, 4.25),\
       (12.70, 0.22),(13.31, 4.92),(15.56, 3.24),(17.82, 4.91),(21.09, -1.69),(23.14, 1.89),(25.19, -1.14),(26.83, -0.58),\
       (28.26, -1.70),(29.49, 0.99),(31.33, -1.03),(33.17, 1.99),(36.04, 0.20),(38.29, -0.81),(42.39, 1.99),(45.05, -0.14),\
       (48.94, 1.20),(53.24, -0.93),(56.52, -0.48),(59.39, -2.84),(68.81, 2.87),(81.71, -3.86),(88.26, -7.78),(92.15, -7.45),\
       (95.02, -11.14),(103.21, -13.27),(106.48, -11.26),(110.17, -12.83),(113.45, -12.94),(117.75, -15.07),\
       (127.17, -15.86),(142.32, -14.30),(148.87, -18.11),(160.34, -16.33),(169.97, -21.26),(179.59, -25.97)'
points = np.array(list(ast.literal_eval(test)))
x = points[:,0]
y = points[:,1]
print(points)

tck = interpolate.splrep(x, y, s=0)
xnew = np.arange(0,180,180/360)
ynew = interpolate.splev(xnew, tck, der=0)
plt.figure()
plt.plot(xnew, ynew)
plt.plot(x, y)
plt.title('Cubic-spline interpolation')
plt.show()