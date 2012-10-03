#!/usr/bin/env python
"""
 This plots and saves figures for rotation scans generated in PhD for paper

 Author: Tom Close (tclose@oist.jp)
 Created: 6/8/2012
"""
#Required Imports
import numpy as np
import os.path
import math
from matplotlib import pyplot as plt
# Configuration parameters (set by generation script)
WORK_DIR = '/home/tclose/Documents/BaFTrS/data/rotate/'
FIGURE_DIR = '/home/tclose/Documents/BaFTrS/figures/'
SNR_RANGE = [2, 5, 10, 20]
NUM_RUNS = 100
CONF_INTERVAL = 0.99
BAR_INTERVAL = 0.05
# Indices for determining the 95% confidence intervals on generated values
low_index = round(NUM_RUNS * BAR_INTERVAL)
high_index = round(NUM_RUNS * (1 - BAR_INTERVAL))
conf_intervals = []
#Loop through each run for each SNR value
fig = plt.figure()
for snr_i, snr in enumerate(SNR_RANGE):
    axes = fig.add_subplot(2, 2, snr_i + 1)
    for run_i in xrange(1, NUM_RUNS):
        rician = np.loadtxt(os.path.join(WORK_DIR, 'rician_snr{snr}_run{run_i}.txt'.format(snr=snr, run_i=run_i)), comments='%')[:, 0]
        gauss = np.loadtxt(os.path.join(WORK_DIR, 'gaussian_snr{snr}_run{run_i}.txt'.format(snr=snr, run_i=run_i)), comments='%')[:, 0]
        assert(len(rician) == len(gauss))
        # Normalise probabilities
        rician -= np.log(np.sum(np.exp(rician)))
        gauss -= np.log(np.sum(np.exp(gauss)))
        diff = rician - gauss
        angles = np.arange(-math.pi / 2, math.pi / 2, math.pi / (len(rician) - 1))
        axes.plot(angles, rician, 'b')
        axes.plot(angles, gauss, 'r')
        axes.plot(angles, diff, 'g')
#        ninety_five = find(rician_scan > log(1 - conf_interval));
#        min_x = x_axis(min(ninety_five));
#        max_x = x_axis(max(ninety_five));
#        conf_intervals(run_i, snr_i) = max_x - min_x;
        if run_i % 100 == 0:
            print "Finished {} runs".format(run_i)
    axes.set_title('SNR {}'.format(snr))
    print "Finished snr {}".format(snr)
plt.show()
