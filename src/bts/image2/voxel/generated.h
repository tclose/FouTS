/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, 22/07/2012.

    This file is part of MRtrix.

    MRtrix is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MRtrix is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MRtrix.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef __bts_image2_voxel_generated_h__
#define __bts_image2_voxel_generated_h__

#include "bts/image2/voxel/reference.h"
#include "bts/image2/interpolator.h"
#include "bts/fibre/segment.h"

namespace BTS {

	namespace Image2 {

    namespace Voxel {

      class Generated : public Reference {

        //Public static variables, nested classes and typedefs
        public:

          static const size_t SEGMENTS_BLOCK_SIZE = 100;

        //Protected member variables
        protected:

          Reference* ref;
          Diffusion::Model* diff_model;
          Interpolator* interp;
          size_t num_segs;
          MR::Math::Matrix<double> segments;
          MR::Math::Matrix<double> work;
          MR::Math::Vector<size_t> fibre_indices;
          MR::Math::Vector<size_t> seg_indices;

        //Public member functions
        public:

          void initialise(size_t num_encodings, const Coord& voxel_centre, Diffusion::Model* diffusion_model,
                                                                 Interpolator* interpolator, Reference* reference = 0) {
            Reference::initialise(num_encodings, voxel_centre);
            ref = reference;
            diff_model = diffusion_model;
            interp = interpolator;
            num_segs = 0;
          }

          void add_segment(const Fibre::Segment& segment, size_t fibre_index, size_t seg_index) {
            assert(segments.columns() == segment.size() || !segments.columns());
            //Mutex to start here
            size_t next_index = num_segs++;
            if (num_segs > segments.rows()) {
              size_t new_num_rows = segments.rows() + SEGMENTS_BLOCK_SIZE;
              segments.resize(new_num_rows, segment.size());
              work.resize(new_num_rows, 5);
              fibre_indices.resize(new_num_rows);
              seg_indices.resize(new_num_rows);
            }
            //Mutex to end here
            segments.row(next_index) = segment;
            fibre_indices[next_index] = fibre_index;
            seg_indices[next_index] = seg_index;

          }

          void  clear_segments() {
            segments.clear();
            fibre_indices.clear();
            seg_indices.clear();
            num_segs = 0;
          }

          void  calc_intensities() {
            MR::Math::Vector<double> interpolations = work.column(3).sub(0,num_segs);
            interp->interpolate(interpolations, seg_positions, centre, work.sub(0,num_segs,0,3));
            MR::Math::mult(segments.sub(6, ))
          }

        //Protected member functions
        protected:

          const MR::Math::Matrix<double>  seg_positions() const
            { return segments.sub(0,num_segs,0,3); }

          const MR::Math::Matrix<double>  seg_tangents() const
            { return segments.sub(0,num_segs,3,6); }

          const MR::Math::Vector<double>  seg_intensities() const
            { return segments.column(6).sub(0,num_segs); }

          /*! Returns the lengths of the segments (from their tangents) in the first column of the work
           *  matrix (so be careful not to reuse the work matrix before you read them).
           *
           * @return The lengths of the stored segments
           */
          MR::Math::Vector<double> seg_lengths() {
            MR::Math::Matrix<double> lengths_work = work.sub(0,num_segs,0,3);
            MR::Math::Vector<double> lengths = lengths_work.column(X);
            lengths_work = seg_tangents();
            lengths_work.column(X) *= lengths_work.column(X);
            lengths_work.column(Y) *= lengths_work.column(Y);
            lengths_work.column(Z) *= lengths_work.column(Z);
            lengths += lengths_work.column(Y); // lengths is column(X) to start with
            lengths += lengths_work.column(Z);
            for (size_t i = 0; i < lengths.size(); ++i)
              lengths[i] = sqrt(lengths[i]);
            return lengths;
          }

      };

    }

	}

}

#endif /* __bts_image2_generated_h__ */
