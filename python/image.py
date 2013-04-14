# Copyright 2013 Florey Neuroscience Institutes
#
# maked by Tom Close on 13/03/09.
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

import collections
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

NumberOfSamples = collections.namedtuple("NumberOfSamples", "longitudinal radial")

class Image(object):

    def __init__(self, size, resolution, encodings):
        if len(size) != 3 or len(resolution) != 3:
            raise Exception("Image must be three-dimensional (provided: size={}, resolution={})."
                            .format(size, resolution))
        if encodings.shape[1] != 4:
            raise Exception("The second dimension of the encodings matrix must have length 4 (found"
                            " {}) corresponding to [x y z b_value] ".format(encodings.shape[1]))
        self._size = size
        self._resolution = resolution
        self._encodings = encodings
        self._data = np.zeros((np.prod(size), encodings.shape[0]))

    # These are provided as read only versions of the 
    @property
    def size(self):
        return self._size

    @property
    def resolution(self):
        return self._resolution

    @property
    def encodings(self):
        return self._encodings

    def generate(self, tract, N_longitudinal, N_radial, response_coeffs):
        samples_matrix = Tract.position_matrix(len(tract), N_longitudinal, N_radial)
        tangents_matrix = Tract.tangent_matrix(len(tract), N_longitudinal, N_radial)
        positions = samples_matrix.dot(tract)
        tangents = tangents_matrix.dot(tract)
        lengths = np.sqrt(tangents[:, 0] ** 2 + tangents[:, 1] ** 2 + tangents[:, 2] ** 2)

        indices = np.floor(positions)
        remain = positions - indices
        disps = ([], [], [])
        for d in xrange(3):
            for i in xrange((1 - INTERP_LENGTH), INTERP_LENGTH):
                disps[d].append(remain[:, d] - i)


class Tract(object):

    _longitudinal_position_matrices = {}
    _longitudinal_tangent_matrices = {}
    _radial_matrices = {}
    _position_matrices = {}
    _tangent_matrices = {}

    @classmethod
    def position_matrix(cls, degree, N_longitudinal, N_radial):
        N = NumberOfSamples(N_longitudinal, N_radial)
        try:
            matrix = cls._position_matrices[(degree, N)]
        except KeyError:
            try:
                long_matrix = cls._longitudinal_position_matrices[(degree, N.longitudinal)]
            except KeyError:
                long_matrix = cls._make_longitudinal_matrix(degree, N)
                cls._longitudinal_position_matrices[(degree, N.longitudinal)] = long_matrix
            try:
                radial_matrix = cls._radial_matrices[(degree, N.radial)]
            except KeyError:
                radial_matrix = cls._make_radial_matrix(degree, N)
                cls._radial_matrices[(degree, N.raial)] = radial_matrix
            matrix = cls._combine_sample_matrices(long_matrix, radial_matrix)
        return matrix

    @classmethod
    def tangent_matrix(cls, degree, N_longitudinal, N_radial):
        N = NumberOfSamples(N_longitudinal, N_radial)
        try:
            matrix = cls._tangent_matrices[(degree, N)]
        except KeyError:
            try:
                long_matrix = cls._longitudinal_tangent_matrices[(degree, N.longitudinal)]
            except KeyError:
                long_matrix = cls._make_tangent_matrix(degree, N)
                cls._longitudinal_tangent_matrices[(degree, N.longitudinal)] = long_matrix
            try:
                radial_matrix = cls._radial_matrices[(degree, N.radial)]
            except KeyError:
                radial_matrix = cls._make_radial_matrix(degree, N)
                cls._radial_matrices[(degree, N.raial)] = radial_matrix
            matrix = cls._combine_sample_matrices(long_matrix, radial_matrix)
        return matrix

    @classmethod
    def _make_longitudinal_matrix(cls, degree, N, include_endpoints=False):
        # Get time increments minus endpoints
        t = np.linspace(0, 1, N + 2)
        if not include_endpoints:
            t = t[1:-1]
        psi = np.ones((N, degree))
        for d in xrange(1, degree):
            psi[:, d] = np.sqrt(2) * np.cos(d * np.pi * t)
        return psi

    @classmethod
    def _make_tangents_matrix(cls, degree, N):
        # Get time increments minus endpoints
        t = np.linspace(0, 1, N + 2)[1:-1] # Remove endpoints which are always [0,0,0]
        dpsi_dt = np.zeros((N, degree))
        for d in xrange(1, degree):
            dpsi_dt[:, d] = -np.sqrt(2) * np.pi * np.sin(d * np.pi * t)
        return dpsi_dt

    @classmethod
    def _make_radial_matrix(cls, N):
        """
        makes pairs of fractions that hexagonally pack the unit circle with a maximal number of 
        samples along its radius of 2 * num_samples + 1
        """
        sub_radius = 1.0 / float(2 * N + 1)
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

    @classmethod
    def _combine_sample__matrices(cls, longitudinal_matrix, radial_matrix):
        length = longitudinal_matrix.shape[0]
        degree = longitudinal_matrix.shape[1]
        samples = np.zeros((length * len(radial_matrix), degree * 3))
        for i, (k1, k2) in enumerate(radial_matrix):
            row_start = i * length
            row_end = row_start + length
            samples[row_start:row_end, :degree] = longitudinal_matrix
            samples[row_start:row_end, degree:(degree * 2)] = longitudinal_matrix * k1
            samples[row_start:row_end, (degree * 2):(degree * 3)] = longitudinal_matrix * k2
        return samples

if __name__ == '__main__':

    DEGREE = 4
    NUM_LENGTH_SAMPLES = 25
    NUM_RADIAL = 4
    INTERP_LENGTH = (1, 1, 1)

    tract = np.zeros((DEGREE * 3, 3))
    tract[1, 0] = 1
    tract[2, 1] = 0.1
    tract[DEGREE, 1] = 0.25
    tract[DEGREE * 2, 2] = 0.25

    samples_matrix = Tract.position_matrix(DEGREE, NUM_LENGTH_SAMPLES, NUM_RADIAL)
    tangents_matrix = Tract.tangent_matrix(DEGREE, NUM_LENGTH_SAMPLES, NUM_RADIAL)

    positions = samples_matrix.dot(tract)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2])
    plt.show()
