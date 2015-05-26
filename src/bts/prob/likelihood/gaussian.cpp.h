/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, 23/09/2010.

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

#ifndef ___bts_prob_likelihood_imagediff_gaussian_cpp_h__
#define ___bts_prob_likelihood_imagediff_gaussian_cpp_h__

namespace FTS {
    
    namespace Prob {
        
        template<typename T> double Likelihood::Gaussian::log_prob_and_fisher_tpl(
                const typename T::Set& fibres, typename T::Set& gradient,
                typename T::Set::Tensor& fisher_info,
                std::vector<typename T::Set::Tensor>& fisher_info_gradients) {
            
            size_t fibres_vector_size = fibres.vsize();
            
            typename Image::Reference::Buffer<typename T::Section>::Set& section_references =
                    exp_image->expected_image_with_references(fibres);
            
            exp_image->precalculate_section_weighting_gradients_and_hessians(T());
            //
            double lprob = 0.0;
            
            gradient = fibres;
            fisher_info = typename T::Set::Tensor(fibres);
            
            gradient.zero();
            fisher_info.zero();
            
            fisher_info_gradients.clear();
            
            for (size_t elem_i = 0; elem_i < fibres_vector_size; elem_i++)
                fisher_info_gradients.push_back(fisher_info);
            
            //      bool has_base_intensity = fibres.has_intrinsic_property("base_intensity");
            //
            //      double base_intensity = 1.0;
            //
            //      if (has_base_intensity)
            //        base_intensity = fibres.intrinsic_property("base_intensity");
            
            std::vector<size_t> fibre_block_start;
            std::vector<size_t> fibre_block_end;
            
            //      size_t block_start = has_base_intensity;
            size_t block_start = 0;
            
            for (size_t fibre_i = 0; fibre_i < fibres.size(); fibre_i++) {
                fibre_block_start.push_back(block_start);
                fibre_block_end.push_back(block_start += fibres[fibre_i].vsize());
            }
            
            //      assert(block_start == fibre_vector_size);
            
            std::set<Image::Index> coords = exp_image->non_empty_or_inbounds();
//      std::set<Image::Index> coords;
//      coords.insert(Image::Index(1,1,1));
            
            //Testing
            
            for (std::set<Image::Index>::iterator index_it = coords.begin();
                    index_it != coords.end(); ++index_it) {
                
                Image::Expected::Voxel& voxel = exp_image->operator()(*index_it);
                
                std::vector<double> d_lprob(exp_image->num_encodings());
                std::vector<double> d2_lprob2(exp_image->num_encodings());
                
                //--------------------------------------------------------------------------------//
                //   Calculate log probability and first and second derivatives for current voxel //
                //--------------------------------------------------------------------------------//
                for (size_t encode_i = 0; encode_i < exp_image->num_encodings(); encode_i++) {
                    
                    if (exp_image->encoding(encode_i).b_value()) {
                        
                        double observed;
                        
                        if (obs_image.in_bounds(*index_it))
                            observed = obs_image(*index_it)[encode_i];
                        else
                            observed = 0;
                        
                        lprob += log_prob(exp_image->operator()(*index_it)[encode_i], observed,
                                d_lprob[encode_i], d2_lprob2[encode_i], *index_it);
                        
                    }
                    
                }
                
                //--------------------------------------------------------------------------//
                //  Precalculate interpolation gradients between sections and current voxel //
                //--------------------------------------------------------------------------//
                
                std::vector<size_t> nonzero_fibres;
                
                for (size_t fibre_i = 0; fibre_i < fibres.size(); fibre_i++)
                    if (!section_references[fibre_i].is_empty(*index_it)) {
                        
                        for (typename std::vector<typename T::Section*>::iterator section_it =
                                section_references[fibre_i](*index_it).begin();
                                section_it != section_references[fibre_i](*index_it).end();
                                ++section_it)
                            voxel.precalculate_interpolation_gradient_and_hessian(**section_it);
                        
                        nonzero_fibres.push_back(fibre_i);
                        
                    }
                
                for (size_t encode_i = 0; encode_i < exp_image->num_encodings(); encode_i++) {
                    if (b0_include == "full" || exp_image->encoding(encode_i).b_value()) {
                        
//            if (has_base_intensity) {
//              gradient.intrinsic_property("base_intensity") -= signal * diff / base_intensity;
//              fisher_info(0,0) -= MR::Math::pow2(signal / base_intensity);
//            }
                        
                        std::vector<MR::Math::Vector<double> > gradient_vectors;
                        std::vector<MR::Math::Matrix<double> > hessian_matrices;
                        
                        //-----------------------------------------------------------------------------------------------------------//
                        //  Combine to calculate gradient and hessian for each fibre with nonzero contributions to the voxels signal //
                        //-----------------------------------------------------------------------------------------------------------//
                        
                        for (std::vector<size_t>::iterator fibre_i_it = nonzero_fibres.begin();
                                fibre_i_it != nonzero_fibres.end(); ++fibre_i_it) {
                            
                            T direction_gradient(fibres[*fibre_i_it]);
                            typename T::Tensor direction_hessian(fibres[*fibre_i_it]);
                            
                            direction_gradient.zero();
                            direction_hessian.zero();
                            
                            std::vector<typename T::Section*>& fibre_section_references =
                                    section_references[*fibre_i_it](*index_it);
                            
                            for (typename std::vector<typename T::Section*>::iterator section_it =
                                    fibre_section_references.begin();
                                    section_it != fibre_section_references.end(); ++section_it) {
                                
                                Fibre::Strand::BasicSection section_gradient;
                                Fibre::Strand::BasicSection::Tensor section_hessian;
                                
                                voxel.direction(encode_i).signal(**section_it, section_gradient,
                                        section_hessian);
                                
                                section_gradient.unnormalize_gradient(exp_image->vox_lengths());
                                section_hessian.unnormalise_hessian(exp_image->vox_lengths());
                                
                                direction_gradient.add_section_gradient(fibres[*fibre_i_it],
                                        **section_it, section_gradient);
                                direction_hessian.add_section_hessian(fibres[*fibre_i_it],
                                        **section_it, section_gradient, section_hessian);
                                
                            }
                            
                            direction_gradient *= fibres.base_intensity();
                            direction_hessian *= fibres.base_intensity();
                            
                            gradient_vectors.push_back(direction_gradient);
                            hessian_matrices.push_back(direction_hessian);
                            
                            gradient[*fibre_i_it] += direction_gradient * d_lprob[encode_i];
                            
                        }
                        
                        for (size_t nonzero_i1 = 0; nonzero_i1 < nonzero_fibres.size();
                                ++nonzero_i1) {
                            
                            size_t start1 = fibre_block_start[nonzero_fibres[nonzero_i1]];
                            MR::Math::Vector<double>& gradient1 = gradient_vectors[nonzero_i1];
                            MR::Math::Matrix<double>& hessian1 = hessian_matrices[nonzero_i1];
                            
                            for (size_t i1 = 0; i1 < gradient1.size(); i1++)
                                for (size_t i2 = 0; i2 < gradient1.size(); i2++) {
                                    fisher_info(start1 + i1, start1 + i2) += d2_lprob2[encode_i]
                                            * gradient1[i1] * gradient1[i2];
//                  std::cout << "fisher_info(" << start1 + i1 << ", " << start1 + i2 << ") = " << fisher_info(start1 + i1, start1 + i2) << " (" << d2_lprob2[encode_i] * gradient1[i1] * gradient1[i2] << ")" << std::endl;
                                }
                            
                            for (size_t i1 = 0; i1 < gradient1.size(); i1++)
                                for (size_t i2 = 0; i2 < gradient1.size(); i2++)
                                    for (size_t i3 = 0; i3 < gradient1.size(); i3++)
                                        fisher_info_gradients[start1 + i1](start1 + i2, start1 + i3) +=
                                                d2_lprob2[encode_i] * (hessian1(i2, i1)
                                                        * gradient1[i3]
                                                                       + gradient1[i2] * hessian1(
                                                                                 i3, i1));
                            
                            for (size_t nonzero_i2 = nonzero_i1 + 1;
                                    nonzero_i2 < nonzero_fibres.size(); ++nonzero_i2) {
                                
                                size_t start2 = fibre_block_start[nonzero_fibres[nonzero_i2]];
                                MR::Math::Vector<double>& gradient2 = gradient_vectors[nonzero_i2];
                                
                                for (size_t i1 = 0; i1 < gradient1.size(); i1++)
                                    for (size_t i2 = 0; i2 < gradient2.size(); i2++) {
                                        
                                        double incr = d2_lprob2[encode_i] * gradient1[i1]
                                                      * gradient2[i2];
                                        
                                        fisher_info(start1 + i1, start2 + i2) += incr;
                                        fisher_info(start2 + i2, start1 + i1) += incr;
                                        
                                    }
                                
                                for (size_t i1 = 0; i1 < gradient1.size(); i1++)
                                    for (size_t i2 = 0; i2 < gradient1.size(); i2++)
                                        for (size_t i3 = 0; i3 < gradient2.size(); i3++) {
                                            
                                            double incr = d2_lprob2[encode_i] * hessian1(i2, i1)
                                                          * gradient2[i3];
                                            
                                            fisher_info_gradients[start1 + i1](start1 + i2,
                                                    start2 + i3) += incr;
                                            fisher_info_gradients[start1 + i1](start2 + i3,
                                                    start1 + i2) += incr;
                                            
                                        }
                                
                            }
                            
                        }
                        
                    }
                }
            }
            
            return lprob;
            
        }
        
        template<typename T> double Likelihood::Gaussian::log_prob_and_fisher_tpl(
                const typename T::Set& fibres, typename T::Set& gradient,
                typename T::Set::Tensor& fisher_info) {
            
            typename Image::Reference::Buffer<typename T::Section>::Set& section_references =
                    exp_image->expected_image_with_references(fibres);
            
            exp_image->precalculate_section_weighting_gradients(T());
            
            double lprob = 0.0;
            
            gradient = fibres;
            fisher_info = typename T::Set::Tensor(fibres);
            
            gradient.zero();
            fisher_info.zero();
            
            std::vector<size_t> fibre_block_start;
            std::vector<size_t> fibre_block_end;
            
            //      size_t block_start = has_base_intensity;
            size_t block_start = 0;
            
            for (size_t fibre_i = 0; fibre_i < fibres.size(); fibre_i++) {
                fibre_block_start.push_back(block_start);
                fibre_block_end.push_back(block_start += fibres[fibre_i].vsize());
            }
            
            //      assert(block_start == fibre_vector_size);
            
            std::set<Image::Index> coords = exp_image->non_empty_or_inbounds();
            
            for (std::set<Image::Index>::iterator index_it = coords.begin();
                    index_it != coords.end(); ++index_it) {
                
                Image::Expected::Voxel& voxel = exp_image->operator()(*index_it);
                
                std::vector<double> d_lprob(exp_image->num_encodings());
                std::vector<double> d2_lprob2(exp_image->num_encodings());
                
                //--------------------------------------------------------------------------------//
                //   Calculate log probability and first and second derivatives for current voxel //
                //--------------------------------------------------------------------------------//
                for (size_t encode_i = 0; encode_i < exp_image->num_encodings(); encode_i++) {
                    
                    if (b0_include == "full" || exp_image->encoding(encode_i).b_value()) {
                        
                        double observed;
                        
                        if (obs_image.in_bounds(*index_it))
                            observed = obs_image(*index_it)[encode_i];
                        else
                            observed = 0;
                        
                        lprob += log_prob(exp_image->operator()(*index_it)[encode_i], observed,
                                d_lprob[encode_i], d2_lprob2[encode_i], *index_it);
                        
                    }
                    
                }
                
                //--------------------------------------------------------------------------//
                //  Precalculate interpolation gradients between sections and current voxel //
                //--------------------------------------------------------------------------//
                for (size_t fibre_i = 0; fibre_i < fibres.size(); fibre_i++)
                    if (!section_references[fibre_i].is_empty(*index_it))
                        for (typename std::vector<typename T::Section*>::iterator section_it =
                                section_references[fibre_i](*index_it).begin();
                                section_it != section_references[fibre_i](*index_it).end();
                                ++section_it)
                            voxel.precalculate_interpolation_gradient(**section_it);
                
                //-------------------------------------------//
                //  Combine to calculate gradient and fisher //
                //-------------------------------------------//
                
                for (size_t encode_i = 0; encode_i < exp_image->num_encodings(); encode_i++) {
                    if (b0_include == "full" || exp_image->encoding(encode_i).b_value()) {
                        
                        //            if (has_base_intensity) {
                        //              gradient.intrinsic_property("base_intensity") -= signal * diff / base_intensity;
                        //              fisher_info(0,0) -= MR::Math::pow2(signal / base_intensity);
                        //            }
                        
                        for (size_t fibre_i1 = 0; fibre_i1 < fibres.size(); fibre_i1++) {
                            if (!section_references[fibre_i1].is_empty(*index_it)) {
                                
                                //--------------------------------------------//
                                //  Calculate the first fibre gradient vector //
                                //--------------------------------------------//
                                T gradient1 = fibres[fibre_i1];
                                gradient1.zero();
                                
                                for (typename std::vector<typename T::Section*>::iterator section_it =
                                        section_references[fibre_i1](*index_it).begin();
                                        section_it != section_references[fibre_i1](*index_it).end();
                                        ++section_it) {
                                    
                                    Fibre::Strand::BasicSection section_gradient;
                                    
                                    voxel.direction(encode_i).signal(**section_it,
                                            section_gradient);
                                    section_gradient.unnormalize_gradient(exp_image->vox_lengths());
                                    gradient1.add_section_gradient(fibres[fibre_i1], **section_it,
                                            section_gradient);
                                    
                                }
                                
                                gradient1 *= fibres.base_intensity();
                                
                                MR::Math::Vector<double>& gradient_vector1 = gradient1;
                                
                                gradient[fibre_i1] += gradient1 * d_lprob[encode_i];
                                
                                //                if (has_base_intensity)
                                //                  for (size_t fibre_elem_i = 0; fibre_elem_i < gradient_vector1.size(); fibre_elem_i++) {
                                //
                                //                    size_t elem_i = fibre_block_start[fibre_i1] + fibre_elem_i;
                                //
                                //                    fisher_info(0, elem_i) -= gradient_vector1[fibre_elem_i] * signal;
                                //                    fisher_info(elem_i, 0) -= gradient_vector1[fibre_elem_i] * signal;
                                //
                                //                  }
                                
                                for (size_t fibre_i2 = 0; fibre_i2 < fibres.size(); fibre_i2++) {
                                    if (!section_references[fibre_i2].is_empty(*index_it)) {
                                        
                                        //---------------------------------------------//
                                        //  Calculate the second fibre gradient vector //
                                        //---------------------------------------------//
                                        T gradient2 = fibres[fibre_i2];
                                        gradient2.zero();
                                        for (typename std::vector<typename T::Section*>::iterator section_it =
                                                section_references[fibre_i2](*index_it).begin();
                                                section_it != section_references[fibre_i2](
                                                        *index_it).end(); ++section_it) {
                                            
                                            Fibre::Strand::BasicSection section_gradient;
                                            
                                            voxel.direction(encode_i).signal(**section_it,
                                                    section_gradient);
                                            section_gradient.unnormalize_gradient(
                                                    exp_image->vox_lengths());
                                            gradient2.add_section_gradient(fibres[fibre_i2],
                                                    **section_it, section_gradient);
                                            
                                        }
                                        
                                        gradient2 *= fibres.base_intensity();
                                        
                                        MR::Math::Vector<double>& gradient_vector2 = gradient2;
                                        
                                        //--------------------------------------------------------------------------------------------//
                                        //  Add the outer product of the two gradient vectors to the appropriate block of the Hessian //
                                        //--------------------------------------------------------------------------------------------//
                                        
                                        for (size_t fibre_elem_i1 = 0;
                                                fibre_elem_i1 < gradient_vector1.size();
                                                fibre_elem_i1++) {
                                            
                                            size_t elem_i1 = fibre_block_start[fibre_i1]
                                                    + fibre_elem_i1;
                                            
                                            for (size_t fibre_elem_i2 = 0;
                                                    fibre_elem_i2 < gradient_vector2.size();
                                                    fibre_elem_i2++) {
                                                
                                                size_t elem_i2 = fibre_block_start[fibre_i2]
                                                        + fibre_elem_i2;
                                                
                                                fisher_info(elem_i1, elem_i2) +=
                                                        gradient_vector1[fibre_elem_i1] * gradient_vector2[fibre_elem_i2]
                                                        * d2_lprob2[encode_i];
                                                
                                            }
                                        }
                                        
                                        size_t dummy = 0;
                                        dummy++;
                                        
                                    }
                                }
                                
                            }
                        }
                        
                    }
                }
            }
            
            return lprob;
            
        }
    
    }

}

#endif /*  */
