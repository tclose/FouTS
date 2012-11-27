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
ALL_RUNS = range(1, 1001)
PLOT_RUNS = range(1, 31)
CONF_INTERVAL = 0.99
NOISE_MODELS = ['rician', 'gaussian', 'diff']
PLOT_COLOURS = {'rician': 'orange', 'gaussian': 'dodgerblue', 'diff': 'indigo'}
SUMMARISE_CONF_INTERVALS = True
for snr_i, snr in enumerate(SNR_RANGE):
    fig = plt.figure()
    axes = fig.add_subplot(111) #fig.add_subplot(2, 2, snr_i + 1)
    # Calculate the minimum value for the current plot
    conf_interval_sizes = {'rician': [], 'gaussian': [], 'diff': []}
    conf_interval_centres = {'rician': [], 'gaussian': [], 'diff': []}
    plot_min = float('inf')
    plot_max = float('-inf')
    for run_i in ALL_RUNS:
        rician = np.loadtxt(os.path.join(WORK_DIR, 'rician_snr{snr}_run{run_i}.txt'.format(snr=snr, run_i=run_i)),
                            comments='%')[:, 0]
        gauss = np.loadtxt(os.path.join(WORK_DIR, 'gaussian_snr{snr}_run{run_i}.txt'.format(snr=snr, run_i=run_i)),
                           comments='%')[:, 0]
        # Count the number of samples
        num_samples = len(rician)
        assert(len(gauss) == num_samples)
        bin_width = math.pi / float(num_samples - 1) # The width of the bins used in the scan
        # Normalise probabilities
        rician -= np.log(np.sum(np.exp(rician)) * bin_width)
        gauss -= np.log(np.sum(np.exp(gauss)) * bin_width)
        angles = np.arange(-math.pi / 2.0, math.pi / 2.0, bin_width)
        diff = rician - gauss
        # Plot the likelihood graphs if it is a run to plot
        if run_i in PLOT_RUNS:
            for prob, nm in zip((rician, gauss, diff), NOISE_MODELS):
                axes.plot(angles, prob, PLOT_COLOURS[nm])
                # Calculate the limits of the plot
                if min(prob) < plot_min:
                    plot_min = min(prob)
                if max(prob) > plot_max:
                    plot_max = max(prob)
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
        # Calculate the median difference and the median centre and plot that for reference only
        diff_ci = ((rician_ci[0] - gauss_ci[0] + num_samples // 2) % num_samples,
                   (rician_ci[1] - gauss_ci[1] + num_samples // 2) % num_samples)
        for noise_model, ci in (('rician', rician_ci), ('gaussian', gauss_ci), ('diff', diff_ci)):
            angles = []
            for i, bin_index in enumerate(ci):
                angles.append((float(bin_index) - float(num_samples) / 2.0) * bin_width)
            conf_interval_sizes[noise_model].append(abs(angles[1] - angles[0]))
            conf_interval_centres[noise_model].append(angles[0] + (angles[1] - angles[0]) / 2.0)
    # Add 5% of the original plot range to both plot min and plot max for buffer around plot
    plot_range = plot_max - plot_min
    plot_max = plot_max + plot_range * 0.05
    plot_min = plot_min - plot_range * 0.05
    print "SNR {}:".format(snr)
    for nm, ci in zip(NOISE_MODELS, (rician_ci, gauss_ci, diff_ci)):
        size = np.average(np.array(conf_interval_sizes[nm]))
        centre = np.average(np.array(conf_interval_centres[nm]))
        print "{nm} - size={size}, centre={centre}".format(nm=nm, size=size, centre=centre)
        for angle in (centre - size / 2.0, centre + size / 2.0):
            axes.plot((angle, angle), (plot_min, plot_max), color=PLOT_COLOURS[nm],
                      linestyle='--')
    axes.set_ylim((plot_min, plot_max))
    axes.set_xlim((-math.pi / 2.0, math.pi / 2.0))
    axes.set_title('SNR {}'.format(snr))
    axes.set_xlabel('Rotation (Radians)')
    axes.set_ylabel('Probability Density')
plt.show()
