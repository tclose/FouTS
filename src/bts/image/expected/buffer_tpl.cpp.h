/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, 16/07/2010.

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

#ifndef __bts_image_expected_buffer_cpp_h__
#define __bts_image_expected_buffer_cpp_h__

#include "bts/image/expected/buffer_tpl.h"
#include "bts/image/reference/buffer.h"

#include "bts/fibre/tractlet/set.h"
#include "bts/fibre/strand/set.h"
#include "bts/fibre/strand/section.h"
#include "bts/fibre/tractlet/section.h"
#include "bts/coord/tensor.h"
#include "bts/fibre/strand/tensor.h"
#include "bts/fibre/tractlet/tensor.h"


namespace BTS {

  namespace Image {

    namespace Expected {

#define LOOP(op) \
  for (size_t z = 0; z < this->dim(Z); z++) { \
    for (size_t y = 0; y < this->dim(Y); y++) { \
      for (size_t x = 0; x < this->dim(X); x++) { \
        op \
      } \
    } \
  }


//***************************************//
//********* Non Hessian Versions ********//
//***************************************//

      template <typename T> template <typename U> void         Buffer_tpl<T>::part_image(const U& fibre) {


        std::vector<typename U::Section> path;

        fibre.sections(path, num_len_sections, num_wth_sections, this->voxel_lengths, this->corner_offsets);

        for (typename std::vector<typename U::Section>::iterator section_it = path.begin(); section_it != path.end(); ++section_it) {

          typename U::Section& section = *section_it;

#ifdef OPTIMISED
          diffusion_model.precalculate_weightings(section);
#endif

          std::set<T*> neighbourhood = this->get_neighbourhood(section.position());

          for (typename std::set<T*>::iterator vox_it = neighbourhood.begin(); vox_it != neighbourhood.end(); ++vox_it) {

            T& voxel = **vox_it;

#ifdef OPTIMISED
            voxel.precalculate_interpolation(section);
#endif

            for (size_t encode_i = 0; encode_i < this->num_encodings(); encode_i++) {

              voxel[encode_i] += voxel.direction(encode_i).signal(section);

            }

          }

        }

      }


      template <typename T> template <typename U> void Buffer_tpl<T>::expected_image(const typename U::Set& fibres) {

        this->zero();

        for (size_t fibre_i = 0; fibre_i < fibres.size(); fibre_i++)
          part_image(fibres[fibre_i]);

        for (typename Buffer_tpl<T>::iterator vox_it = this->begin(); vox_it != this->end(); ++vox_it)
          for (size_t encode_i = 0; encode_i < num_encodings(); ++encode_i)
            vox_it->second[encode_i] *= fibres.base_intensity();

      }



      template <typename T> template <typename U> void         Buffer_tpl<T>::part_image(const U& fibre,
                                                                    std::vector<typename U::Section>& path,
                                                                    Reference::Buffer<typename U::Section>& section_reference) {


        fibre.sections(path, num_len_sections, num_wth_sections, this->voxel_lengths, this->corner_offsets, this->num_encodings());

        for (typename std::vector<typename U::Section>::iterator section_it = path.begin(); section_it != path.end(); ++section_it) {

          typename U::Section& section = *section_it;

#ifdef OPTIMISED
          diffusion_model.precalculate_weightings(section);
#endif

          std::set<T*> neighbourhood = this->get_neighbourhood(section.position());

          for (typename std::set<T*>::iterator vox_it = neighbourhood.begin(); vox_it != neighbourhood.end(); ++vox_it) {

            T& voxel = **vox_it;

            section_reference(voxel.coord()).push_back(&section);

#ifdef OPTIMISED
            voxel.precalculate_interpolation(section);
#endif

            for (size_t encode_i = 0; encode_i < this->num_encodings(); encode_i++) {

              voxel[encode_i] += voxel.direction(encode_i).signal(section);

            }

          }

        }

      }


      template <typename T> template <typename U> typename Image::Reference::Buffer<typename U::Section>::Set& Buffer_tpl<T>::expected_image_with_references(const typename U::Set& fibres) {

        this->zero();

        std::map<size_t, std::vector<typename U::Section> >& sections = get_sections(U());
        typename Reference::Buffer<typename U::Section>::Set& section_refs = get_section_references(U());

        for (size_t fibre_i = 0; fibre_i < fibres.size(); fibre_i++) {

          section_refs[fibre_i].clear_references();

          part_image(fibres[fibre_i], sections[fibre_i], section_refs[fibre_i]);

        }

        for (typename Buffer_tpl<T>::iterator vox_it = this->begin(); vox_it != this->end(); ++vox_it)
          for (size_t encode_i = 0; encode_i < num_encodings(); ++encode_i)
            vox_it->second[encode_i] *= fibres.base_intensity();


        return section_refs;


      }


      template <typename T> template <typename U>
              void         Buffer_tpl<T>::part_image( const U& fibre,
                                                                Container::Buffer<U>& gradients) {


        std::vector<typename U::Section> path;

        fibre.sections(path, num_len_sections, num_wth_sections, this->voxel_lengths, this->corner_offsets);

        for (typename std::vector<typename U::Section>::iterator section_it = path.begin(); section_it != path.end(); ++section_it) {

          typename U::Section& section = *section_it;

#ifdef OPTIMISED
          diffusion_model.precalculate_weightings_and_gradients(section);
#endif

          std::set<T*> neighbourhood = this->get_neighbourhood(section.position());

          for (typename std::set<T*>::iterator vox_it = neighbourhood.begin(); vox_it != neighbourhood.end(); ++vox_it) {

            T& voxel = **vox_it;
            const Index& coord = voxel.coord();

            if (gradients.is_empty(coord))
              gradients(coord).reset(fibre);

            Image::Container::Voxel<U>& gradient_voxel = gradients(coord);


#ifdef OPTIMISED
            voxel.precalculate_interpolation_gradient(section);
#endif

            for (size_t encode_i = 0; encode_i < this->num_encodings(); encode_i++) {

              typename U::Section          section_gradient;

              voxel[encode_i] += voxel.direction(encode_i).signal(section, section_gradient);

#ifndef GRADIENT_NOT_REQUIRED
              section_gradient.unnormalize_gradient(this->vox_lengths());
              gradient_voxel[encode_i].add_section_gradient(*section.parent, section, section_gradient);
#endif

            }

          }

        }

      }


      template <typename T> template <typename U> void Buffer_tpl<T>::expected_image(const typename U::Set& fibres,
                                                       typename Container::Buffer<U>::Set& gradients) {


        this->zero();

        for (size_t fibre_i = 0; fibre_i < fibres.size(); fibre_i++) {

          if (gradients.size() <= fibre_i)
            gradients.push_back(Container::Buffer<U> (this->dimensions, this->num_encodings()));
          else if (gradients[fibre_i].dims() != this->dims() || gradients[fibre_i].num_encodings() != this->num_encodings())
            gradients[fibre_i].reset(this->dimensions, this->num_encodings());
          else
            //NB: This does not actually clear the voxels, but instead sets a flag that designates it them as being empty.
            // Therefore the function 'Buffer_tpl<T>::is_empty(const Image::Index&)' should be used to check before using a given voxel.
            gradients[fibre_i].quick_clear();


          part_image(fibres[fibre_i], gradients[fibre_i]);


        }


        for (typename Buffer_tpl<T>::iterator vox_it = this->begin(); vox_it != this->end(); ++vox_it)
          for (size_t encode_i = 0; encode_i < num_encodings(); ++encode_i)
            vox_it->second[encode_i] *= fibres.base_intensity();


        for (size_t fibre_i = 0; fibre_i < fibres.size(); ++fibre_i)
          for (typename Buffer_tpl<T>::iterator vox_it = this->begin(); vox_it != this->end(); ++vox_it)
            for (size_t encode_i = 0; encode_i < num_encodings(); ++encode_i)
                gradients[fibre_i](vox_it->first)[encode_i] *= fibres.base_intensity();

      }


      template <typename T> template <typename U> void     Buffer_tpl<T>::precalculate_section_weighting_gradients() {

        std::map<size_t, std::vector<typename U::Section> >& sections = get_sections(U());

        for (size_t fibre_i = 0; fibre_i < sections.size(); ++fibre_i)
          for (size_t section_i = 0; section_i < sections[fibre_i].size(); ++section_i)
            diffusion_model.precalculate_weightings_and_gradients(sections[fibre_i][section_i]);

      }


      template <typename T> template <typename U> void     Buffer_tpl<T>::precalculate_section_weighting_gradients_and_hessians() {

        std::map<size_t, std::vector<typename U::Section> >& sections = get_sections(U());

        for (size_t fibre_i = 0; fibre_i < sections.size(); ++fibre_i)
          for (size_t section_i = 0; section_i < sections[fibre_i].size(); ++section_i)
            diffusion_model.precalculate_weightings_gradients_and_hessians(sections[fibre_i][section_i]);

      }



//***********************************//
//********* Hessian Versions ********//
//***********************************//



      template <typename T> template <typename U>
        void         Buffer_tpl<T>::part_image( const U& fibre,
                                                          Container::Buffer<U>& gradients,
                                                          Container::Buffer<typename U::Tensor>& hessians) {


        std::vector<typename U::Section> path;

        fibre.sections(path, num_len_sections, num_wth_sections, this->voxel_lengths, this->corner_offsets);

        for (typename std::vector<typename U::Section>::iterator section_it = path.begin(); section_it != path.end(); ++section_it) {

          typename U::Section& section = *section_it;

#ifdef OPTIMISED
          diffusion_model.precalculate_weightings_gradients_and_hessians(section);
#endif

          std::set<T*> neighbourhood = this->get_neighbourhood(section.position());

          for (typename std::set<T*>::iterator vox_it = neighbourhood.begin(); vox_it != neighbourhood.end(); ++vox_it) {

            T& voxel = **vox_it;
            const Index& coord = voxel.coord();

            if (gradients.is_empty(coord))
              gradients(coord).reset(fibre);

            Image::Container::Voxel<U>& gradient_voxel = gradients(coord);

            if (hessians.is_empty(coord))
              hessians(coord).reset(typename U::Tensor(fibre));

            Image::Container::Voxel<typename U::Tensor>& hessian_voxel = hessians(coord);


#ifdef OPTIMISED
            voxel.precalculate_interpolation_gradient_and_hessian(section);
#endif

            for (size_t encode_i = 0; encode_i < this->num_encodings(); encode_i++) {

              typename U::Section          section_gradient;
              typename U::Section::Tensor  section_hessian;

              voxel[encode_i] += voxel.direction(encode_i).signal(section, section_gradient, section_hessian);

#ifndef GRADIENT_NOT_REQUIRED
              section_gradient.unnormalize_gradient(this->vox_lengths());
              gradient_voxel[encode_i].add_section_gradient(*section.parent, section, section_gradient);

#ifndef HESSIAN_NOT_REQUIRED
              section_hessian.unnormalise_hessian(this->vox_lengths());
              hessian_voxel[encode_i].add_section_hessian(*section.parent, section, section_gradient, section_hessian);
#endif

#endif

            }

          }

        }

      }


      template <typename T> template <typename U>
          void Buffer_tpl<T>::expected_image(const typename U::Set& fibres,
                                                       typename Container::Buffer<U>::Set& gradients,
                                                       typename Container::Buffer<typename U::Tensor>::Set& hessians) {


        this->zero();

        for (size_t fibre_i = 0; fibre_i < fibres.size(); fibre_i++) {

          if (gradients.size() <= fibre_i)
            gradients.push_back(Container::Buffer<U> (this->dimensions, this->num_encodings()));
          else if (gradients[fibre_i].dims() != this->dims() || gradients[fibre_i].num_encodings() != this->num_encodings())
            gradients[fibre_i].reset(this->dimensions, this->num_encodings());
          else
            //NB: This does not actually clear the voxels, but instead sets a flag that designates it them as being empty.
            // Therefore the function 'Buffer_tpl<T>::is_empty(const Image::Index&)' should be used to check before using a given voxel.
            gradients[fibre_i].quick_clear();


          if (hessians.size() <= fibre_i)
            hessians.push_back(Container::Buffer<typename U::Tensor> (this->dimensions, this->num_encodings()));
          else if (hessians[fibre_i].dims() != this->dims() || hessians[fibre_i].num_encodings() != this->num_encodings())
            hessians[fibre_i].reset(this->dimensions, this->num_encodings());
          else
            //NB: This does not actually clear the voxels, but instead sets a flag that designates them as being empty.
            // Therefore the function 'Buffer_tpl<T>::is_empty(const Image::Index&)' should be used to check before using a given voxel.
            hessians[fibre_i].quick_clear();


          part_image(fibres[fibre_i], gradients[fibre_i], hessians[fibre_i]);


        }

        for (typename Buffer_tpl<T>::iterator vox_it = this->begin(); vox_it != this->end(); ++vox_it)
          for (size_t encode_i = 0; encode_i < num_encodings(); ++encode_i)
            vox_it->second[encode_i] *= fibres.base_intensity();


        for (size_t fibre_i = 0; fibre_i < fibres.size(); ++fibre_i)
          for (typename Buffer_tpl<T>::iterator vox_it = this->begin(); vox_it != this->end(); ++vox_it)
            for (size_t encode_i = 0; encode_i < num_encodings(); ++encode_i) {
              gradients[fibre_i](vox_it->first)[encode_i] *= fibres.base_intensity();
              hessians[fibre_i](vox_it->first)[encode_i] *= fibres.base_intensity();
            }


      }


      template <typename T> Buffer_tpl<T>&                             Buffer_tpl<T>::operator+=(const Expected::Buffer& buff) {

        if (this->dims() != buff.dims())
          throw Exception ("Buffer dimensions do not match, " + str(buff.dims()) + " and " + str(this->dims()) + ".");

        std::set<Index> non_empty = buff.non_empty();

        for (std::set<Index>::iterator coord_it = non_empty.begin(); coord_it != non_empty.end(); ++coord_it)
          this->operator()(*coord_it) += buff(*coord_it);

        return *this;

      }

      template <typename T> Buffer_tpl<T>&                             Buffer_tpl<T>::operator-=(const Expected::Buffer& buff) {

        if (this->dims() != buff.dims())
          throw Exception ("Buffer dimensions do not match, " + str(buff.dims()) + " and " + str(this->dims()) + ".");

        std::set<Index> non_empty = buff.non_empty();

        for (std::set<Index>::iterator coord_it = non_empty.begin(); coord_it != non_empty.end(); ++coord_it)
          this->operator()(*coord_it) -= buff(*coord_it);

        return *this;

      }

      template <typename T> Buffer_tpl<T>&       Buffer_tpl<T>::operator-=(const Image::Buffer_tpl<Observed::Voxel>& buff) {

        if (this->dims() != buff.dims())
          throw Exception ("Buffer dimensions do not match, " + str(buff.dims()) + " and " + str(this->dims()) + ".");

        std::set<Index> non_empty = buff.non_empty();

        for (std::set<Index>::iterator coord_it = non_empty.begin(); coord_it != non_empty.end(); ++coord_it)
          this->operator()(*coord_it) -= buff(*coord_it);

        return *this;

      }


      template <typename T> std::set<T*>&                              Buffer_tpl<T>::get_neighbourhood(const Coord& point) {

        Index centre_coord = voxel_centre_coord(point);

        std::set<T*>& neighbourhood = neighbourhoods[centre_coord];

        if (!neighbourhood.size())
          neighbourhoods[centre_coord] = create_neighbourhood(centre_coord);

        return neighbourhood;

      }

      template <typename T> template <typename U> std::set<T*>         Buffer_tpl<T>::get_neighbourhood(const typename U::Section& section) {

        std::set<T*> whole_neighbourhood;

        std::vector<Coord> extreme_points = T::extreme_points(section);

        for (std::vector<Coord>::iterator point_it = extreme_points.begin(); point_it != extreme_points.end(); ++point_it) {
          std::set<T*> point_neighbourhood = get_neighbourhood(voxel_centre_coord(*point_it));
          whole_neighbourhood.insert(point_neighbourhood.begin(), point_neighbourhood.end());
        }

        return whole_neighbourhood;

      }


      template <typename T> std::set<T*>                                Buffer_tpl<T>::create_neighbourhood(const Index& coord) {

        std::set<T*> neighbourhood;

        // Loop through all voxels withinness the interpolation extent of the given voxel and add their address to the
        // neighbourhood vector.
        for (int z_neigh = coord[Z] - neigh_extent; z_neigh < coord[Z] + neigh_extent; z_neigh++)
          for (int y_neigh = coord[Y] - neigh_extent; y_neigh < coord[Y] + neigh_extent; y_neigh++)
            for (int x_neigh = coord[X] - neigh_extent; x_neigh < coord[X] + neigh_extent; x_neigh++) {

              Index coord(x_neigh, y_neigh, z_neigh);

              if (!this->enforce_bounds || this->in_bounds(coord))
                neighbourhood.insert(&(this->operator()(coord)));

            }

        return neighbourhood;

      }


      template <typename T> std::ostream&       Buffer_tpl<T>::to_stream (std::ostream& stream) const {

        this->Observed::Buffer_tpl<T>::to_stream(stream);

        stream << "Diffusion model: " << diffusion_model << std::endl;
        stream << "Num sections: " << num_len_sections << std::endl;
        stream << "Num strands: " << num_wth_sections << std::endl;
        stream << "Interpolation extent: " << interp_extent << std::endl;
        stream << "Neighbourhood extent: " << neigh_extent << std::endl;

        stream << "Neighbourhoods: " << std::endl;

        for (typename std::map< Index, std::set<T*> >::const_iterator neigh_it = neighbourhoods.begin(); neigh_it != neighbourhoods.end(); ++neigh_it) {

          stream << neigh_it->first << std::endl;

          for (typename std::set<T*>::const_iterator vox_it = neigh_it->second.begin(); vox_it != neigh_it->second.end(); ++vox_it)
            stream << "  " << (*vox_it)->centre() << std::endl;

        }

        stream << std::endl << std::endl;

        return stream;

      }


      template <typename T> std::ostream&       operator<< (std::ostream& stream, const Buffer_tpl<T>& buffer) {

        return buffer.to_stream(stream);

      }


      template <typename T>  double             Buffer_tpl<T>::get_base_intensity(double ref_b0) {

        double base_intensity = 0.0;
        if (ref_b0) {
          Coord interp_length = this->interp_extent * this->vox_lengths();
          // Create tract spans the interpolation length of the interpolation kernel which is centred on the
          // bottom left voxel.
          Fibre::Tractlet::Set tcts (1,2);
          tcts.zero();
          tcts[0](0,0) = this->corner_offsets + this->vox_lengths() / 2.0;
          tcts[0](0,1) = Coord(interp_length[0] / M_SQRT2, 0.0, 0.0);
          tcts[0](1,0) = Coord(0.0, interp_length[1] * M_SQRT2, 0.0);
          tcts[0](2,0) = Coord(0.0, 0.0, interp_length[2] * M_SQRT2);
          // Normalize the density of the tract and set the base_intensity of the set to 1.0, to calculate the required
          // base intensity value to match that of the reference.
          tcts.normalise_densities();
          tcts.set_base_intensity(1.0);
          expected_image<Fibre::Tractlet>(tcts);
          // Divide the reference b0 by the value in the test voxel in the bottom left corner.
          base_intensity = ref_b0 / this->operator()(0,0,0).b0();
          clear();
        }
        return base_intensity;

      }
    }

  }


}

#undef LOOP


#endif /* __bts_image_expected_buffer_tpl_cpp_h__ */
