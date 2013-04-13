# Copyright 2013 Florey Neuroscience Institutes
#
# Created by Tom Close on 13/03/09.
#
# This file is part of the Fourier Tract Sampling (FouTS) package.
#
# FouTS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FouTS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with FouTS.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def create_along_samples(degree, num_samples, include_endpoints=False):
    # Get time increments minus endpoints
    t = np.linspace(0, 1, num_samples + 2)
    if not include_endpoints:
        t = t[1:-1]
    psi = np.ones((num_samples, degree))
    for d in xrange(1, degree):
        psi[:, d] = np.sqrt(2) * np.cos(d * np.pi * t)
    return psi

def create_tangent_samples(degree, num_samples):
    # Get time increments minus endpoints
    t = np.linspace(0, 1, num_samples + 2)[1:-1] # Remove endpoints which are always [0,0,0]
    dpsi_dt = np.zeros((num_samples, degree))
    for d in xrange(1, degree):
        dpsi_dt[:, d] = -np.sqrt(2) * np.pi * np.sin(d * np.pi * t)
    return dpsi_dt

def create_across_samples(num_samples):
    """
    Creates pairs of fractions that hexagonally pack the unit circle with a maximal number of 
    samples along its radius of 2 * num_samples + 1
    """
    sub_radius = 1.0 / float(2 * num_samples + 1)
    incr1 = sub_radius * 2.0
    incr2 = sub_radius * np.sqrt(3)
    samples = []
    offset_row = False
    for step2 in np.arange(0, 1, incr2):
        offset = sub_radius if offset_row else 0.0
        for k2 in ((-step2, step2) if step2 else (0,)):
            for step1 in (np.arange(0, 1, incr1) + offset):
                for k1 in ((-step1, step1) if step1 else (0,)):
                    if k1 ** 2 + k2 ** 2 <= (1.0 - sub_radius):
                        samples.append((k1, k2))
        offset_row = not offset_row
    return samples

def create_samples(degree, num_along_samples, num_across_samples):

    along_samples = create_along_samples(degree, num_along_samples)
    across_samples = create_across_samples(num_across_samples)
    samples = np.zeros((num_along_samples * len(across_samples), degree * 3))
    for i, (k1, k2) in enumerate(across_samples):
        row_start = i * num_along_samples
        row_end = row_start + num_along_samples
        samples[row_start:row_end, :degree] = along_samples
        samples[row_start:row_end, degree:(degree * 2)] = along_samples * k1
        samples[row_start:row_end, (degree * 2):(degree * 3)] = along_samples * k2
    return samples

def create_tangents(degree, num_tangent_samples, num_across_samples):

    tangent_samples = create_tangent_samples(degree, num_tangent_samples)
    across_samples = create_across_samples(num_across_samples)
    samples = np.zeros((num_tangent_samples * len(across_samples), degree * 3))
    for i, (k1, k2) in enumerate(across_samples):
        row_start = i * num_tangent_samples
        row_end = row_start + num_tangent_samples
        samples[row_start:row_end, :degree] = tangent_samples
        samples[row_start:row_end, degree:(degree * 2)] = tangent_samples * k1
        samples[row_start:row_end, (degree * 2):(degree * 3)] = tangent_samples * k2
    return samples

if __name__ == '__main__':

    DEGREE = 4
    NUM_LENGTH_SAMPLES = 25
    NUM_ACROSS_SAMPLES = 4
    INTERP_LENGTH = (1, 1, 1)

    tract = np.zeros((DEGREE * 3, 3))
    tract[1, 0] = 1
    tract[2, 1] = 0.1
    tract[DEGREE, 1] = 0.25
    tract[DEGREE * 2, 2] = 0.25

    samples_matrix = create_samples(DEGREE, NUM_LENGTH_SAMPLES, NUM_ACROSS_SAMPLES)
    tangents_matrix = create_tangents(DEGREE, NUM_LENGTH_SAMPLES, NUM_ACROSS_SAMPLES)

    positions = samples_matrix.dot(tract)
    tangents = tangents_matrix.dot(tract)
    lengths = np.sqrt(tangents[:, 0] ** 2 + tangents[:, 1] ** 2 + tangents[:, 2] ** 2)

    indices = np.floor(positions)
    remain = positions - indices
    disps = ([], [], [])
    for d in xrange(3):
        for i in xrange((1 - INTERP_LENGTH), INTERP_LENGTH):
            disps[d].append(remain[:, d] - i)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2])
    plt.show()
