import numpy as np
import matplotlib.pyplot as plt
for sub in ('donald', 'heath', 'lisa'):
    if sub == 'donald':
        ndirs = (150, 60, 20)
    else:
        ndirs = (150, 60)
    for ndir in ndirs:
        a = np.loadtxt('{sub}/intens_stats.{ndir}.txt'.format(sub=sub, ndir=ndir))
        plt.figure()
        plt.hist(a, bins=100, range=(0.0, 150))
        plt.title("{} - {}".format(sub, ndir))
plt.show()

