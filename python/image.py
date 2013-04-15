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
from copy import copy
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

NumberOfSamples = collections.namedtuple("NumberOfSamples", "long radial")


class Tract(object):

    def __init__(self, degree, acs=1.0, N_long=25, N_radial=4):
        self._degree = degree
        self._data = np.zeros(degree * 9 + 1)
        self.acs = acs
        self.N = NumberOfSamples(N_long, N_radial)

    @property
    def degree(self):
        return self._degree

    @property
    def matrix(self):
        return self._data[:-1].reshape((-1, 3))

    @property
    def vector(self):
        return self._data[:]

    @property
    def acs(self):
        return self._data[-1] ** 2

    @acs.setter
    def acs(self, acs):
        self._data[-1] = np.sqrt(acs)

    def __getitem__(self, indices):
        if len(indices) == 1:
            item = self._data[indices[0]]
        elif len(indices) == 2:
            offset = (self._degree * indices[0] + indices[1]) * 3
            item = self._data[offset:(offset + 3)]
        elif len(indices) == 3:
            item = self._data[(self._degree * indices[0] + indices[1]) * 3 + indices[2]]
        return item

    def __setitem__(self, indices, value):
        if len(indices) == 1:
            self._data[indices[0]] = value
        elif len(indices) == 2:
            offset = (self._degree * indices[0] + indices[1]) * 3
            self._data[offset:(offset + 3)] = value
        elif len(indices) == 3:
            self._data[(self._degree * indices[0] + indices[1]) * 3 + indices[2]] = value

    def positions(self):
        pos_matrix = Tract.position_matrix(self.degree, self.N)
        return pos_matrix.dot(self.matrix)

    def tangents(self):
        tang_matrix = Tract.tangent_matrix(self.degree, self.N)
        return tang_matrix.dot(self.matrix)

    def segments(self, scale=(1, 1, 1), offset=(0, 0, 0)):
        # Get segment positions, scaling the tract if required
        if not all(scale) or any(offset):
            norm_matrix = copy(self.matrix)
            norm_matrix[0, 0, :] += offset
            norm_matrix *= scale
        else:
            norm_matrix = self.matrix
        positions = Tract.position_matrix(self.degree, self.N).dot(norm_matrix)
        # Get segment orientations
        tangents = self.tangents()
        lengths = np.sqrt(tangents[:, 0] ** 2 + tangents[:, 1] ** 2 + tangents[:, 2] ** 2) # lengths
        # Normalise orientations
        orientations = tangents # rename for readability
        for dim_i in xrange(3):
            orientations[:, dim_i] /= lengths
        # Scale the segment sizes by the ACS and divide by the number of samples over the cross section
        sizes = lengths # rename for readability
        sizes *= (tract.acs / tract.num_strands())
        return (positions, orientations, sizes)

    def num_strands(self):
        return len(self._radial_matrix(self.N.radial))

    @classmethod
    def num_samples_per_cross_section(cls, N_radial):
        try:
            num = len(cls._radial_matrices[N_radial])
        except KeyError:
            cls._radial_matrices[N_radial] = cls._make_radial_matrix(N_radial)
            num = len(cls._radial_matrices[N_radial])
        return num


    # Class Methods ------------------------------------------------------------------------------

    _long_pos_matrices = {}
    _long_tang_matrices = {}
    _radial_matrices = {}
    _pos_matrices = {}
    _tang_matrices = {}

    @classmethod
    def _make_long_matrix(cls, degree, N, include_endpoints=False):
        # Get time increments minus endpoints
        t = np.linspace(0, 1, N + 2)
        if not include_endpoints:
            t = t[1:-1]
        psi = np.ones((N, degree))
        for d in xrange(1, degree):
            psi[:, d] = np.sqrt(2) * np.cos(d * np.pi * t)
        return psi

    @classmethod
    def _make_tang_matrix(cls, degree, N):
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
    def _combine_sample_matrices(cls, long_matrix, radial_matrix):
        length = long_matrix.shape[0]
        degree = long_matrix.shape[1]
        samples = np.zeros((length * len(radial_matrix), degree * 3))
        for i, (k1, k2) in enumerate(radial_matrix):
            row_start = i * length
            row_end = row_start + length
            samples[row_start:row_end, :degree] = long_matrix
            samples[row_start:row_end, degree:(degree * 2)] = long_matrix * k1
            samples[row_start:row_end, (degree * 2):(degree * 3)] = long_matrix * k2
        return samples

    @classmethod
    def _radial_matrix(cls, N_radial):
        try:
            radial_matrix = cls._radial_matrices[N_radial]
        except KeyError:
            radial_matrix = cls._make_radial_matrix(N_radial)
            cls._radial_matrices[N_radial] = radial_matrix
        return radial_matrix

    @classmethod
    def _long_matrix(cls, degree, N_long):
        try:
            long_matrix = cls._long_pos_matrices[(degree, N_long)]
        except KeyError:
            long_matrix = cls._make_long_matrix(degree, N_long)
            cls._long_pos_matrices[(degree, N_long)] = long_matrix
        return long_matrix

    @classmethod
    def _tang_matrix(cls, degree, N_long):
        try:
            tang_matrix = cls._tang_matrices[(degree, N_long)]
        except KeyError:
            tang_matrix = cls._make_tang_matrix(degree, N_long)
            cls._tang_matrices[(degree, N_long)] = tang_matrix
        return tang_matrix

    @classmethod
    def position_matrix(cls, degree, N):
        try:
            matrix = cls._pos_matrices[(degree, N)]
        except KeyError:
            matrix = cls._combine_sample_matrices(cls._long_matrix(degree, N.long),
                                                                   cls._radial_matrix(N.radial))
        return matrix

    @classmethod
    def tangent_matrix(cls, degree, N):
        try:
            matrix = cls._tang_matrices[(degree, N)]
        except KeyError:
            matrix = cls._combine_sample_matrices(cls._tang_matrix(degree, N.long),
                                                                   cls._radial_matrix(N.radial))
        return matrix

class Interpolator(object):

    def __init__(self, extent):
        self.extent = extent

    def range(self, dim_index):
        return np.arange(1 - self.extent[dim_index], self.extent[dim_index] + 1)

    def interpolate(self, X, Y, Z, indices=None):
        raise NotImplementedError("Classes deriving the 'Interpolator' class must implement the "
                                  "'interpolate' method")


class SincInterpolator(Interpolator):

    def interpolate(self, X, Y, Z, indices=None): #@UnusedVariable indices -> uniform interpolation
        self._interp(X, self.extent[0])
        self._interp(Y, self.extent[1])
        self._interp(Z, self.extent[2])

    def _interp(self, D, extent):
        # To stop the interpolation kernel blowing up at 0
        zero_mask = (D == 0.0)
        # Truncate outside of the kernel extent
        trunc_mask = (D > extent) * (D < extent)
        # Get a mask of the remaining displacements
        mask = np.invert(zero_mask + trunc_mask)
        # Change the values of the displacement matrix in-place to the interpolated values
        D[mask] = np.sin(D[mask]) / D[mask]
        D[zero_mask] = 1.0
        D[trunc_mask] = 0.0


class Image(object):

    def __init__(self, size, resolution, offset, directions):
        if len(size) != 3 or len(resolution) != 3:
            raise Exception("Image must be three-dimensional (provided: size={}, resolution={})."
                            .format(size, resolution))
        if directions.shape[1] != 3:
            raise Exception("The second dimension of the encodings matrix must have length 3 (found"
                            " {}) corresponding to [x y z b_value] ".format(directions.shape[1]))
        self._size = np.asarray(size, dtype=int)
        self._resolution = np.asarray(resolution)
        self._offset = np.asarray(offset)
        self._directions = np.asarray(directions)
        self._data = np.zeros(size + directions.shape[0:1])

    # These are provided as read-only versions of the member variables
    @property
    def size(self):
        return self._size

    @property
    def resolution(self):
        return self._resolution

    @property
    def offset(self):
        return self._offset

    @property
    def directions(self):
        return self._directions

    def generate(self, tract, dw_response=(1.0, 0.5, 0.2), interpolator=SincInterpolator((1, 1, 1)),
                 base_intensity=1.0):
        self._data[:] = 0.0
        (positions, orientations, segment_sizes) = tract.segments(scale=(1.0 / self.resolution),
                                                                  offset= -self.offset)
        # Calculate the diffusion weighting from the segment orientations in each DW direction
        Q = orientations.dot(self.directions.T)
        Q2 = Q ** 2
        QY = np.ones(Q.shape)
        R = np.zeros(Q.shape)
        for coeff in dw_response[:-1]:
            R += QY * coeff
            QY *= Q2
        R += QY * coeff # This is performed outside the loop to avoid the last multiplication of QY
        # Transpose the DW responses for broadcasting with interpolation vectors
        R = R.transpose()
        # Scale the diffusion weightings by the segment sizes (the amount of the tract the signal arises from)
        R *= segment_sizes
        # Calculate the interpolation from the segment positions to their neighbouring voxels
        # along each of the dimensions
        floor = np.floor(positions)
        # Get the voxel indices that are affected by each segment
        base_indices = np.array(floor, dtype=np.int)
#        X_ind = base_indices[:, 0:1] + interpolator.range(0)
#        Y_ind = base_indices[:, 1:2] + interpolator.range(1)
#        Z_ind = base_indices[:, 2:3] + interpolator.range(2)
#        # Get image indices corresponding to segment->voxel contributions 
#        indices = (X_ind.T[interpolator.extent[0], None, None, :] +
#                   Y_ind.T[None, interpolator.extent[1], None, :] * self.size[1] +
#                   Z_ind.T[None, None, interpolator.extent[2], :] * np.prod(self.size[0:2])).ravel()
        # Get the displacement to surrounding voxels
        disp = positions - floor # Since the segment positions are already normalised to the image
        X = disp[:, 0:1] - interpolator.range(0)
        Y = disp[:, 1:2] - interpolator.range(1)
        Z = disp[:, 2:3] - interpolator.range(2)
        interpolator.interpolate(X, Y, Z)
#        # Experimental way of doing this that might be quicker but should definitely be paralleizable
#        H = (X.T[interpolator.extent[0], None, None, :] *
#             Y.T[None, interpolator.extent[1], None, :] *
#             Z.T[None, None, interpolator.extent[2], :])
#        # Multiply responses by interpolations to get the contributions to voxels from segments
#        S = np.ravel(R * H)
#        # Assign the contributions to the relevant voxels (this sorting could be done in a separate 
#        # thread while the other calculations are going on)
#        unique_indices = np.unique(indices)
#        index_indices = np.argsort(indices) # indices of indices
#        sorted_indices = indices[index_indices]
#        left = np.searchsorted(sorted_indices, unique_indices, 'left')
#        right = np.searchsorted(sorted_indices, unique_indices, 'right')
#        for ind, l, r in zip(unique_indices, left, right):
#            self._data[ind, :] = np.sum(S[index_indices[l:r], :], axis=0)
        # Combine interpolation and responses to get the generated signal intensities the old way,
        # which might be a little easier to understand
        for z_i, z_offset in enumerate(interpolator.range(2)):
            z_ind = base_indices[:, 2] + z_offset
            # Map indices that are out of range onto the last index (therefore it is 
            # important to ensure this is an empty voxel outside the true image bounds
            z_ind[(z_ind < 0) * (z_ind > self.size[2])] = -1
            for y_i, y_offset in enumerate(interpolator.range(1)):
                YZ = Y[:, y_i] * Z[:, z_i]
                y_ind = base_indices[:, 1] + y_offset
                # Map indices that are out of range onto the last index (therefore it is 
                # important to ensure this is an empty voxel outside the true image bounds
                y_ind[(y_ind < 0) * (y_ind > self.size[2])] = -1
                for x_i, x_offset in enumerate(interpolator.range(0)):
                    XYZ = X[:, x_i] * YZ
                    x_ind = base_indices[:, 0] + x_offset
                    # Map indices that are out of range onto the last index (therefore it is 
                    # important to ensure this is an empty voxel outside the true image bounds 
                    x_ind[(x_ind < 0) * (x_ind > self.size[2])] = -1
                    self._data[x_ind, y_ind, z_ind, :] += (R * XYZ).transpose()
        self._data *= base_intensity



if __name__ == '__main__':

    DEGREE = 4
    NUM_LENGTH_SAMPLES = 25
    NUM_RADIAL = 4
    INTERP_LENGTH = (1, 1, 1)

    tract = Tract(DEGREE)
    tract[0, 1, 0] = 1
    tract[0, 2, 1] = 0.1
    tract[1, 0, 1] = 0.25
    tract[2, 0, 2] = 0.25

    directions = np.random.rand(10, 3)
    length = directions[:, 0] ** 2 + directions[:, 1] ** 2 + directions[:, 2] ** 2
    directions[:] = (directions.T / length).T

    image = Image((10, 10, 10), (0.1, 0.1, 0.1), (0.0, 0.0, 0.0), directions)
    image.generate(tract)

    for i in range(10):
        fig = plt.figure()
        plt.imshow(image._data[:, :, i, 5])
#    ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2])
    plt.show()
