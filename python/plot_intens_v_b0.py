#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os.path

images = {'donald': [150, 60, 20], 'heath':[150, 60], 'lisa': [150, 60]}

tract_points = []
vox_points = []
for subject in images.keys():
    for n in images[subject]:
        b0 = np.loadtxt('/home/tclose/fouts/params/image/reference/{}/corpus_callosum.{}.b0.txt'
                .format(subject, n))
        for i in range(1, 6):
            intens = np.loadtxt('/home/tclose/fouts/params/image/reference/{}/intens.{}.{}.txt'.format(subject, n, i))
            tract_points.append((b0, intens))
        intens = np.loadtxt('/home/tclose/fouts/params/image/reference/{}/intensity.{}.new.txt'.format(subject, n, i))
        vox_points.append((b0, intens))

fig1 = plt.figure()
tract_points = np.array(tract_points)
plt.scatter(tract_points[:, 0], tract_points[:, 1])
fig2 = plt.figure()
vox_points = np.array(vox_points)
plt.scatter(vox_points[:, 0], vox_points[:, 1])
plt.show()


