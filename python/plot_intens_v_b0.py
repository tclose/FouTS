#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os.path

images = {'donald': [150, 60, 20], 'heath':[150, 60], 'lisa': [150, 60]}

points = []
for subject in images.keys():
    for n in images[subject]:
        b0 = np.loadtxt('/home/tclose/fouts/params/image/reference/{}/corpus_callosum.{}.b0.txt'
                .format(subject, n))
        for i in range(1, 6):
            intens = np.loadtxt('/home/tclose/fouts/params/image/reference/{}/intens.{}.{}.txt'.format(subject, n, i))
            points.append((b0, intens))

points = np.array(points)
plt.scatter(points[:, 0], points[:, 1])
plt.show()


