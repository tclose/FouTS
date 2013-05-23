import numpy as np
import matplotlib.pyplot as plt
intens = np.loadtxt('intensities.txt')
b0s = np.loadtxt('b0s.txt')
plt.scatter(b0s, intens)
plt.show()

