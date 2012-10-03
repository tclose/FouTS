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
PLOT_COLOURS = {'Rician': 'b', 'Gaussian': 'g', 'Difference': 'r'}
#Loop through each run for each SNR value
# Create proxy artists for legend
plot_index_mins = []
for snr_i, snr in enumerate(SNR_RANGE):
    run_plot_index_mins = [] # The minimums for each index of every plot
    for run_i in xrange(1, NUM_RUNS):
        rician = np.loadtxt(os.path.join(WORK_DIR, 'rician_snr{snr}_run{run_i}.txt'.format(snr=snr, run_i=run_i)), comments='%')[:, 0]
        gauss = np.loadtxt(os.path.join(WORK_DIR, 'gaussian_snr{snr}_run{run_i}.txt'.format(snr=snr, run_i=run_i)), comments='%')[:, 0]
        # Count the number of samples
        num_samples = len(rician)
        assert(len(gauss) == num_samples)
        bin_width = math.pi / float(num_samples - 1) # The width of the bins used in the scan
        bin_width = math.pi / float(num_samples - 1) # The width of the bins used in the scan
        # Normalise probabilities
        rician -= np.log(np.sum(np.exp(rician)) * bin_width)
        gauss -= np.log(np.sum(np.exp(gauss)) * bin_width)
        run_index_mins = np.minimum(rician, gauss)
        run_plot_index_mins = np.minimum(run_plot_index_mins, run_index_mins) if len(run_plot_index_mins) else run_index_mins
    plot_index_mins.append(run_plot_index_mins)
plot_mins = [ np.min(x) for x in plot_index_mins]
for snr_i, snr in enumerate(SNR_RANGE):
    fig = plt.figure()
    axes = fig.add_subplot(111) #fig.add_subplot(2, 2, snr_i + 1)
    # Calculate the minimum value for the current plot
    for run_i in xrange(1, NUM_RUNS):
        rician = np.loadtxt(os.path.join(WORK_DIR, 'rician_snr{snr}_run{run_i}.txt'.format(snr=snr, run_i=run_i)), comments='%')[:, 0]
        gauss = np.loadtxt(os.path.join(WORK_DIR, 'gaussian_snr{snr}_run{run_i}.txt'.format(snr=snr, run_i=run_i)), comments='%')[:, 0]
        # Count the number of samples
        num_samples = len(rician)
        assert(len(gauss) == num_samples)
        bin_width = math.pi / float(num_samples - 1) # The width of the bins used in the scan
        # Normalise probabilities
        rician -= np.log(np.sum(np.exp(rician)) * bin_width)
        gauss -= np.log(np.sum(np.exp(gauss)) * bin_width)
        angles = np.arange(-math.pi / 2.0, math.pi / 2.0, bin_width)
        # Plot the differences
        diff = rician - gauss
        axes.plot(angles, rician, PLOT_COLOURS['Rician'])
        axes.plot(angles, gauss, PLOT_COLOURS['Gaussian'])
        axes.plot(angles, diff, PLOT_COLOURS['Difference'])
        # Calculate the confidence intervals
        conf_intervals = []
        for dist in (gauss, rician):
            peak_index = np.argmax(dist)
            left_index = peak_index
            right_index = peak_index
            prob = math.exp(dist[peak_index]) * bin_width
            step_count = 0
            while prob < CONF_INTERVAL:
                # Get proposed new left and right indices to potentially add to the interval
                prop_left_index = (left_index - 1) % num_samples
                prop_right_index = (right_index + 1) % num_samples
                # Get the values at the proposed indices
                left_val = math.exp(dist[prop_left_index])
                right_val = math.exp(dist[prop_right_index])
                # Add the index that will increase the probability the most to the interval
                if left_val > right_val:
                    prob += left_val * bin_width
                    left_index = prop_left_index
                else:
                    prob += right_val * bin_width
                    right_index = prop_right_index
                step_count += 1
                if step_count > num_samples:
                    left_index = 0
                    right_index = num_samples - 1
                    break
            # Calculate the angles of the the interval start and ends
            conf_intervals.append((left_index, right_index))
        rician_ci, gauss_ci = conf_intervals
        diff_ci = ((rician_ci[0] - gauss_ci[0] + num_samples // 2) % num_samples, (rician_ci[1] - gauss_ci[1] + num_samples // 2) % num_samples)
        for ci, colour in zip((rician_ci, gauss_ci, diff_ci), (PLOT_COLOURS['Rician'], PLOT_COLOURS['Gaussian'],
                                                                                            PLOT_COLOURS['Difference'])):
            for index, style in zip(ci, ('--', ':')):
                angle = (float(index) - float(num_samples) / 2.0) * bin_width
                axes.plot((angle, angle), (plot_mins[snr_i], plot_index_mins[snr_i][index]), style + colour)
        if run_i % 100 == 0:
            print "Finished {} runs".format(run_i)
    axes.set_title('SNR {}'.format(snr))
    axes.set_xlabel('Rotation (Radians)')
    axes.set_ylabel('Probability Density')
    print "Finished snr {}".format(snr)
plt.show()
