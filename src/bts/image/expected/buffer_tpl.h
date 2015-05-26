/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, Mar 15, 2011.

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

#ifndef __bts_image_expected_buffertpl_h__
#define __bts_image_expected_buffertpl_h__

namespace FTS {
    
    namespace Image {
        
        namespace Expected {
            
            template<typename T> class Buffer_tpl;
        
        }
    
    }

}

#include "bts/common.h"

#include "bts/image/observed/buffer.h"
#include "bts/image/index.h"
#include "bts/image/container/buffer.h"
#include "bts/image/reference/buffer.h"

#include "bts/fibre/strand.h"
#include "bts/fibre/tractlet.h"

#include "bts/diffusion/model.h"

#include "math/matrix.h"

#include "bts/common.h"

#include "bts/image/index.h"
#include "bts/image/buffer.h"
#include "bts/fibre/strand.h"
#include "bts/fibre/tractlet.h"

#include "bts/diffusion/model.h"

#include "bts/image/observed/buffer.h"
#include "bts/image/container/buffer.h"
#include "bts/image/expected/buffer.h"

#define LOOP(op) \
for (size_t z = 0; z < this->dim(Z); z++) { \
  for (size_t y = 0; y < this->dim(Y); y++) { \
    for (size_t x = 0; x < this->dim(X); x++) { \
      op \
    } \
  } \
}

namespace FTS {
    
    namespace Image {
        
        namespace Expected {
            
            template<typename T> class Buffer_tpl: public Observed::Buffer_tpl<T> {
                    
                protected:
                    
                    Diffusion::Model diffusion_model;

                    size_t num_len_sections;
                    size_t num_wth_sections;

                    double interp_extent;
                    int neigh_extent;

                    std::map<Index, std::set<T*> > neighbourhoods;

                public:
                    //Holds all strand or tractlet sections used to generate the image. Storage here allows data structures to be
                    //reused between image generations. The size_t index indicates the strand index
                    std::map<size_t, std::vector<Fibre::Strand::Section> > strand_sections;
                    std::map<size_t, std::vector<Fibre::Tractlet::Section> > tractlet_sections;

                    //Holds all strand or tractlet sections used to generate the image. Storage here allows data structures to be
                    //reused between image generations.
                    Reference::Buffer<Fibre::Strand::Section>::Set strand_section_references;
                    Reference::Buffer<Fibre::Tractlet::Section>::Set tractlet_section_references;

                public:
                    
                    size_t num_length_sections() const {
                        return num_len_sections;
                    }
                    
                    size_t num_width_sections() const {
                        return num_wth_sections;
                    }
                    
                    template<typename U> void expected_image(const typename U::Set& fibres);

                    template<typename U> void part_image(const U& fibre);

                    template<typename U> typename Image::Reference::Buffer<typename U::Section>::Set& expected_image_with_references(
                            const typename U::Set& fibres);

                    template<typename U> void part_image(
                            const U& fibre, std::vector<typename U::Section>& path,
                            Image::Reference::Buffer<typename U::Section>& section_reference);

                    template<typename U> void expected_image(
                            const typename U::Set& fibres,
                            typename Container::Buffer<U>::Set& gradients);

                    template<typename U> void part_image(const U& fibre,
                                                         Container::Buffer<U>& gradients);

                    template<typename U> void expected_image(
                            const typename U::Set& fibres,
                            typename Container::Buffer<U>::Set& gradients,
                            typename Container::Buffer<typename U::Tensor>::Set& hessians);

                    template<typename U> void part_image(
                            const U& fibre, Container::Buffer<U>& gradients,
                            Container::Buffer<typename U::Tensor>& hessians);

                    template<typename U> void precalculate_section_weighting_gradients();

                    double get_base_intensity(double ref_b0);

                    template<typename U> void
                    precalculate_section_weighting_gradients_and_hessians();

                    const Diffusion::Model& get_diffusion_model() const {
                        return diffusion_model;
                    }
                    
                    double get_extent() {
                        return interp_extent;
                    }
                    
                    void set_extent(double interp_extent) {
                        this->interp_extent = interp_extent;
                        neigh_extent = (size_t) std::ceil(interp_extent);
                        neighbourhoods.clear();
                    }
                    
                    Buffer_tpl& clear() {
                        Image::Buffer_tpl<T>::clear();
                        neighbourhoods.clear();
                        return *this;
                    }
                    
                    size_t num_encodings() const {
                        return diffusion_model.num_encodings();
                    }
                    
                    const Diffusion::Encoding& encoding(size_t index) const {
                        return diffusion_model[index];
                    }
                    
                protected:
                    
                    Buffer_tpl(bool enforce_bounds = true)
                            : Observed::Buffer_tpl<T>(enforce_bounds), num_len_sections(0), num_wth_sections(
                                      0), interp_extent(0.0), neigh_extent(0) {
                    }
                    
                    Buffer_tpl(const Triple<size_t>& dimensions, const Triple<double>& voxel_sizes,
                               const Diffusion::Model& diffusion_model,
                               size_t number_length_sections, size_t number_width_sections,
                               double interp_extent, const Triple<double>& corner_offsets,
                               bool enforce_bounds)
                            : Observed::Buffer_tpl<T>(dimensions, voxel_sizes, corner_offsets,
                                      enforce_bounds), diffusion_model(diffusion_model), num_len_sections(
                                      number_length_sections), num_wth_sections(
                                      number_width_sections)

                    {
                        
                        if (!number_length_sections)
                            throw Exception("Number of length sections in image cannot be zero.");
                        
                        if (!number_width_sections)
                            throw Exception("Number of length sections in image cannot be zero.");
                        
                        set_extent(interp_extent);
                        
                    }
                    
                    Buffer_tpl(const Buffer_tpl& bt)
                            : Observed::Buffer_tpl<T>(bt), diffusion_model(bt.diffusion_model), num_len_sections(
                                      bt.num_len_sections), num_wth_sections(bt.num_wth_sections)

                    {
                        set_extent(bt.interp_extent);
                    }
                    
                    virtual ~Buffer_tpl() {
                    }
                    
                    // Returns the 'neighbourhood' coordinate that the point lies in.  The 'neighbourhood' coordinate is the coordinate
                    // that the point + [0.5,0.5,0.5] lies in. This is used to map the neighbourhood of voxels that lie
                    // withinness the interpolation extent about the given point (see set and get neighbourhoods).
                    Index voxel_centre_coord(const Coord& point) {
                        Coord offset_point = point + Coord::Halves;
                        Index coord((int) floor(offset_point[X]), (int) floor(offset_point[Y]),
                                (int) floor(offset_point[Z]));
                        return coord;
                    }
                    
                    std::set<T*> create_neighbourhood(const Index& coord);

                    std::set<T*>& get_neighbourhood(const Coord& point);

                    template<typename U> std::set<T*> get_neighbourhood(
                            const typename U::Section& section);

                    std::ostream& to_stream(std::ostream& stream) const;

                    //The dummy arguments to the following functions is used in template functions to specify which type of section
                    //to return.
                    
                    std::map<size_t, std::vector<Fibre::Strand::Section> >& get_sections(
                            const Fibre::Strand& dummy) {
                        return strand_sections;
                    }
                    
                    std::map<size_t, std::vector<Fibre::Tractlet::Section> >& get_sections(
                            const Fibre::Tractlet& dummy) {
                        return tractlet_sections;
                    }
                    
                    Reference::Buffer<Fibre::Strand::Section>::Set& get_section_references(
                            const Fibre::Strand& dummy) {
                        return strand_section_references;
                    }
                    
                    Reference::Buffer<Fibre::Tractlet::Section>::Set& get_section_references(
                            const Fibre::Tractlet& dummy) {
                        return tractlet_section_references;
                    }
                    
                    Buffer_tpl& operator+=(const Expected::Buffer& buffer);

                    Buffer_tpl& operator-=(const Expected::Buffer& buffer);

                    Buffer_tpl& operator-=(const Image::Buffer_tpl<Observed::Voxel>& buffer);
                    
            };
        
        }
    
    }

}

#undef LOOP

#endif /*  */
