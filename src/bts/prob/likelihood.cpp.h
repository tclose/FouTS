/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, 30/06/2010.

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

#ifndef __bts_prob_likelihood_imagediff_cpp_h__
#define __bts_prob_likelihood_imagediff_cpp_h__

#include "bts/math/common.h"

#include "bts/fibre/strand/set/tensor.h"
#include "bts/fibre/tractlet/set/tensor.h"
#include "bts/image/expected/quartic/buffer.h"

namespace BTS {
    
    namespace Prob {
        
        //TODO: Collapse the workings of the different log_prob_tpl's into common subfunctions.
        template<typename T> double Likelihood::log_prob_tpl(const typename T::Set& fibres) {
            
            exp_image->expected_image(fibres);
            
            return log_prob(*exp_image);
            
        }
        
        template<typename T> double Likelihood::log_prob_tpl(const typename T::Set& fibres,
                                                             typename T::Set& gradient) {
            
            typename Image::Reference::Buffer<typename T::Section>::Set& section_references =
                    exp_image->expected_image_with_references(fibres);
            
            exp_image->precalculate_section_weighting_gradients(T());
            
            double lprob = 0.0;
            
            gradient = fibres;
            gradient.zero();
            
            std::set<Image::Index> coords = exp_image->non_empty_or_inbounds();
            
            for (std::set<Image::Index>::iterator coord_it = coords.begin();
                    coord_it != coords.end(); ++coord_it) {
                
                Image::Expected::Voxel& voxel = exp_image->operator()(*coord_it);
                
                std::vector<double> d_lprob(exp_image->num_encodings());
                
                for (size_t encode_i = 0; encode_i < exp_image->num_encodings(); encode_i++) {
                    
                    if (b0_include == "full" || exp_image->encoding(encode_i).b_value()) {
                        
                        double observed;
                        
                        if (obs_image.in_bounds(*coord_it))
                            observed = obs_image(*coord_it)[encode_i];
                        else
                            observed = 0;
                        
                        lprob += log_prob(exp_image->operator()(*coord_it)[encode_i], observed,
                                d_lprob[encode_i]);
                        
                    }
                    
                }
                
                for (size_t fibre_i = 0; fibre_i < fibres.size(); fibre_i++) {
                    
                    if (!section_references[fibre_i].is_empty(*coord_it)) {
                        
                        std::vector<typename T::Section*>& voxel_section_references =
                                section_references[fibre_i](*coord_it);
                        
                        for (typename std::vector<typename T::Section*>::iterator section_it =
                                voxel_section_references.begin();
                                section_it != voxel_section_references.end(); ++section_it) {
                            
                            typename T::Section& section = **section_it;
                            typename T::Section section_gradient;
                            
                            voxel.precalculate_interpolation_gradient(section);
                            
                            for (size_t encode_i = 0; encode_i < exp_image->num_encodings();
                                    encode_i++) {
                                
                                if (b0_include == "full" || exp_image->encoding(encode_i).b_value()) {
                                    
                                    voxel.direction(encode_i).signal(section, section_gradient);
                                    
                                    section_gradient.unnormalize_gradient(exp_image->vox_lengths());
                                    
                                    //Pre-multiply the section gradient with the lprob derivative, so it can be added directly to the
                                    //fibre gradient.
                                    section_gradient *= d_lprob[encode_i];
                                    
                                    gradient[fibre_i].add_section_gradient(fibres[fibre_i], section,
                                            section_gradient);
                                    
                                }
                                
                            }
                            
                        }
                        
                    }
                    
                }
                
            }
            
            return lprob;
            
        }
        
        template<typename T> double Likelihood::log_prob_tpl(const typename T::Set& fibres,
                                                             typename T::Set& gradient,
                                                             typename T::Set::Tensor& hessian) {
            
            typename Image::Container::Buffer<T>::Set& buffer_gradients =
                    get_recycled_gradients<T>();
            typename Image::Container::Buffer<typename T::Tensor>::Set& buffer_hessians =
                    get_recycled_hessians<T>();
            
            exp_image->expected_image(fibres, buffer_gradients, buffer_hessians);
            
            double lprob = 0.0;
            
//      size_t fibres_vector_size = fibres.vsize();
//      size_t fibres_vector_size = fibres.vsize() - 1;
            
            //TODO: Add zero(T::Set& template_value) overload
            gradient = fibres;
            hessian = typename T::Set::Tensor(fibres);
            
            gradient.zero();
            hessian.zero();
            
//      bool has_base_intensity = fibres.has_intrinsic_property("base_intensity");
            
            std::vector<size_t> fibre_block_start;
            std::vector<size_t> fibre_block_end;
            
//      size_t block_start = has_base_intensity;
            size_t block_start = 0;
            
            for (size_t fibre_i = 0; fibre_i < fibres.size(); fibre_i++) {
                fibre_block_start.push_back(block_start);
                fibre_block_end.push_back(block_start += fibres[fibre_i].vsize());
            }
            
            std::set<Image::Index> coords = exp_image->non_empty();
            
            for (size_t x = 0; x < exp_image->dim(X); ++x)
                for (size_t y = 0; y < exp_image->dim(Y); ++y)
                    for (size_t z = 0; z < exp_image->dim(Z); ++z)
                        coords.insert(Image::Index(x, y, z));
            
            for (std::set<Image::Index>::iterator coord_it = coords.begin();
                    coord_it != coords.end(); ++coord_it) {
                
                for (size_t encode_i = 0; encode_i < exp_image->num_encodings(); encode_i++) {
                    
                    if (b0_include == "full" || exp_image->encoding(encode_i).b_value()) {
                        
                        double observed;
                        
                        if (obs_image.in_bounds(*coord_it))
                            observed = obs_image(*coord_it)[encode_i];
                        else
                            observed = 0;
                        
                        double d_lprob, d2_lprob2;
                        
                        lprob += log_prob(exp_image->operator()(*coord_it)[encode_i], observed,
                                d_lprob, d2_lprob2);
                        
                        typename T::Set direction_gradient(fibres);
                        direction_gradient.zero();
                        
                        typename T::Set::Tensor hessian_fibre_blocks(fibres);
                        hessian_fibre_blocks = 0.0;
                        
                        for (size_t fibre_i = 0; fibre_i < fibres.size(); fibre_i++) {
                            
                            if (!buffer_gradients[fibre_i].is_empty(*coord_it)) {
                                
                                direction_gradient[fibre_i] =
                                        buffer_gradients[fibre_i](*coord_it)[encode_i];
                                
                                MR::Math::Matrix<double> hessian_fibre_block_view =
                                        hessian_fibre_blocks.sub(fibre_block_start[fibre_i],
                                                fibre_block_end[fibre_i],
                                                fibre_block_start[fibre_i],
                                                fibre_block_end[fibre_i]);
                                
                                MR::Math::Matrix<double> hessian_fibre_block =
                                        buffer_hessians[fibre_i](*coord_it)[encode_i];
                                
                                hessian_fibre_block *= d_lprob;
                                
                                hessian_fibre_block_view = hessian_fibre_block;
                                
                            } else
                                direction_gradient[fibre_i] = T(fibres[fibre_i].size(), 0.0);
                            
                        }
                        
                        MR::Math::Matrix<double> hessian_outer_prod = Math::outer(
                                direction_gradient, direction_gradient);
                        
                        hessian_outer_prod *= d2_lprob2;
                        
                        hessian += hessian_fibre_blocks;
                        hessian += hessian_outer_prod;
                        
                        gradient += direction_gradient * d_lprob;
                        
                    }
                    
                }
                
            }
            
            return lprob;
            
        }
        
        //! Used to get the right used gradients for the templated type.
        template<> inline Image::Container::Buffer<Fibre::Strand>::Set& Likelihood::get_recycled_gradients<
                Fibre::Strand>() {
            
            return recycled_strand_gradients;
            
        }
        
        //! Used to get the right used gradients for the templated type.
        template<> inline Image::Container::Buffer<Fibre::Tractlet>::Set& Likelihood::get_recycled_gradients<
                Fibre::Tractlet>() {
            
            return recycled_tractlet_gradients;
            
        }
        
        //! Used to get the right used hessians for the templated type.
        template<> inline Image::Container::Buffer<Fibre::Strand::Tensor>::Set& Likelihood::get_recycled_hessians<
                Fibre::Strand>() {
            
            return recycled_strand_hessians;
            
        }
        
        //! Used to get the right used hessians for the templated type.
        template<> inline Image::Container::Buffer<Fibre::Tractlet::Tensor>::Set& Likelihood::get_recycled_hessians<
                Fibre::Tractlet>() {
            
            return recycled_tractlet_hessians;
            
        }
    
    }

}

#endif /* __bts_prob_likelihood_imagediff_cpp_h__ */
