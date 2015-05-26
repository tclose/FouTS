/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Created by Tom Close on 13/03/09.

 This file is part of Fourier Tract Sampling (FouTS).

 FouTS is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 FouTS is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with FTS.  If not, see <http://www.gnu.org/licenses/>.

 */

#include "bts/prob/prior_component/end_on_sphere.h"

#include "bts/image/buffer.cpp.h"

namespace FTS {

  namespace Prob {

    namespace PriorComponent {

        const std::string       EndOnSphere::NAME  = "end_on_sphere";
        const double            EndOnSphere::SCALE_DEFAULT  = 0.0; //100.0;
        const double            EndOnSphere::RADIUS_DEFAULT  = 0.225;


        double              EndOnSphere::radius_default(const Image::Observed::Buffer& obs_image) {

          return radius_default(obs_image.dim(X), obs_image.vox_length(X));
        }

        double              EndOnSphere::radius_default(size_t dim, double vox_length) {

          return (double)dim * vox_length * 2.0 / 3.0;

    //      return (double)dim * vox_length / 2.0;
    //      return (double)dim * vox_length * SQRT_3 / 2.0;

        }


        double      EndOnSphere::log_prob(const Fibre::Strand& strand, Fibre::Strand& gradient) {

          double log_prob = 0.0;
          gradient.zero();

          std::pair< Triple<double>, Triple<double> > endpoints = strand.endpoints();

          double norm1 = endpoints.first.norm();
          double norm2 = endpoints.second.norm();

          double dist_from_sphere1  = norm1 - radius;
          double dist_from_sphere2  = norm2 - radius;

#ifndef GRADIENT_NOT_REQUIRED
          std::pair<MR::Math::Vector<double>, MR::Math::Vector<double> > coeffs = strand.endpoint_coeffs();
#endif


          log_prob -= MR::Math::pow2(dist_from_sphere1) * scale;
          log_prob -= MR::Math::pow2(dist_from_sphere2) * scale;

#ifndef GRADIENT_NOT_REQUIRED
          for (size_t degree_i = 0; degree_i < strand.degree(); degree_i++) {
            gradient[degree_i] -= (2.0 * scale * dist_from_sphere1 * coeffs.first[degree_i]/ norm1) * endpoints.first;
            gradient[degree_i] -= (2.0 * scale * dist_from_sphere2 * coeffs.second[degree_i]/ norm2) * endpoints.second;
          }
#endif

          return log_prob;


        }

//
//        double      EndOnSphere::log_prob(const Fibre::Strand& strand, Fibre::Strand& gradient, Fibre::Strand::Tensor& hessian) {
//
//          double log_prob = 0.0;
//
//          gradient.resize(strand.degree());
//          gradient.zero();
//
//          hessian.resize(strand.degree());
//          hessian.zero();
//
//          std::pair< Triple<double>, Triple<double> > endpoints = strand.endpoints();
//
//          double norm1 = endpoints.first.norm();
//          double norm2 = endpoints.second.norm();
//
//          double dist_from_sphere1  = norm1 - radius;
//          double dist_from_sphere2  = norm2 - radius;
//
//#ifndef GRADIENT_NOT_REQUIRED
//          std::pair<MR::Math::Vector<double>, MR::Math::Vector<double> > coeffs = strand.endpoint_coeffs();
//#endif
//
//
//          log_prob -= MR::Math::pow2(dist_from_sphere1) * scale;
//          log_prob -= MR::Math::pow2(dist_from_sphere2) * scale;
//
//#ifndef GRADIENT_NOT_REQUIRED
//          for (size_t degree_i = 0; degree_i < strand.degree(); degree_i++) {
//            gradient[degree_i] -= (2.0 * scale * dist_from_sphere1 * coeffs.first[degree_i]/ norm1) * endpoints.first;
//            gradient[degree_i] -= (2.0 * scale * dist_from_sphere2 * coeffs.second[degree_i]/ norm2) * endpoints.second;
//          }
//
//
//#ifndef HESSIAN_NOT_REQUIRED
//
//
//
//          for (size_t degree_i1 = 0; degree_i1 < strand.degree(); degree_i1++) {
//            for (size_t dim_i1 = 0; dim_i1 < 3; dim_i1++) {
//              for (size_t degree_i2 = 0; degree_i2 < strand.degree(); degree_i2++) {
//                for (size_t dim_i2 = 0; dim_i2 < 3; dim_i2++) {
//
//                  double& hess_elem = hessian[degree_i1][dim_i1][degree_i2][dim_i2];
//
//                  hess_elem -= 2.0 * scale * coeffs.first[degree_i1] * coeffs.first[degree_i2] *
//                                  endpoints.first[dim_i1] * endpoints.first[dim_i2] *
//                                  (1.0/MR::Math::pow2(norm1) - dist_from_sphere1/MR::Math::pow3(norm1));
//
//                  hess_elem -= 2.0 * scale * coeffs.second[degree_i1] * coeffs.second[degree_i2] *
//                                  endpoints.second[dim_i1] * endpoints.second[dim_i2] *
//                                  (1.0/MR::Math::pow2(norm2) - dist_from_sphere2/MR::Math::pow3(norm2));
//
//
//                  if ((dim_i1 == dim_i2)) {
//                    hess_elem -=  2.0 * scale * coeffs.first[degree_i1] * coeffs.first[degree_i2] * dist_from_sphere1 / norm1;
//                    hess_elem -=  2.0 * scale * coeffs.second[degree_i1] * coeffs.second[degree_i2] * dist_from_sphere2 / norm2;
//                  }
//
//                }
//              }
//            }
//          }
//
//
//#endif
//
//#endif
//
//          return log_prob;
//
//
//        }


//
//        double      EndOnSphere::log_prob_and_fisher(const Fibre::Strand& strand, Fibre::Strand& gradient, Fibre::Strand::Tensor& fisher) {
//
//          double log_prob = 0.0;
//
//          gradient.zero();
//
//          fisher.resize(strand.degree());
//          fisher.zero();
//
//          std::pair< Triple<double>, Triple<double> > endpoints = strand.endpoints();
//
//          double norm1 = endpoints.first.norm();
//          double norm2 = endpoints.second.norm();
//
//          double dist_from_sphere1  = norm1 - radius;
//          double dist_from_sphere2  = norm2 - radius;
//
//#ifndef GRADIENT_NOT_REQUIRED
//          std::pair<MR::Math::Vector<double>, MR::Math::Vector<double> > coeffs = strand.endpoint_coeffs();
//#endif
//
//
//          log_prob -= MR::Math::pow2(dist_from_sphere1) * scale;
//          log_prob -= MR::Math::pow2(dist_from_sphere2) * scale;
//
//#ifndef GRADIENT_NOT_REQUIRED
//          for (size_t degree_i = 0; degree_i < strand.degree(); degree_i++) {
//            gradient[degree_i] -= (2.0 * scale * dist_from_sphere1 * coeffs.first[degree_i]/ norm1) * endpoints.first;
//            gradient[degree_i] -= (2.0 * scale * dist_from_sphere2 * coeffs.second[degree_i]/ norm2) * endpoints.second;
//          }
//
//
//          for (size_t degree_i1 = 0; degree_i1 < strand.degree(); degree_i1++) {
//            for (size_t dim_i1 = 0; dim_i1 < 3; dim_i1++) {
//              for (size_t degree_i2 = 0; degree_i2 < strand.degree(); degree_i2++) {
//                for (size_t dim_i2 = 0; dim_i2 < 3; dim_i2++) {
//
//                  double& fish_elem = fisher[degree_i1][dim_i1][degree_i2][dim_i2];
//
//                  fish_elem +=  2.0 * scale * coeffs.first[degree_i1] * coeffs.first[degree_i2] *
//                                  endpoints.first[dim_i1] * endpoints.first[dim_i2] / MR::Math::pow2(norm1);
//
//                  fish_elem +=  2.0 * scale * coeffs.second[degree_i1] * coeffs.second[degree_i2] *
//                                  endpoints.second[dim_i1] * endpoints.second[dim_i2] / MR::Math::pow2(norm2);
//
//                }
//              }
//            }
//          }
//
//#endif
//
//          return log_prob;
//
//
//        }


//
//        double      EndOnSphere::log_prob_and_fisher(const Fibre::Strand& strand, Fibre::Strand& gradient, Fibre::Strand::Tensor& fisher, std::vector<Fibre::Strand::Tensor>& fisher_gradient) {
//
//          size_t vector_size = strand.vector_size();
//          size_t has_density = strand.has_property("density");
//
//          double log_prob = 0.0;
//
//          gradient.zero();
//
//          fisher.resize(strand.degree());
//          fisher.zero();
//
//          fisher_gradient.resize(vector_size);
//          for (size_t elem_i = 0; elem_i < vector_size; ++elem_i) {
//            fisher_gradient[elem_i].resize(strand.degree());
//            fisher_gradient[elem_i].zero();
//          }
//
//          std::pair< Triple<double>, Triple<double> > endpoints = strand.endpoints();
//
//          double norm1 = endpoints.first.norm();
//          double norm2 = endpoints.second.norm();
//
//          double dist_from_sphere1  = norm1 - radius;
//          double dist_from_sphere2  = norm2 - radius;
//
//#ifndef GRADIENT_NOT_REQUIRED
//          std::pair<MR::Math::Vector<double>, MR::Math::Vector<double> > coeffs = strand.endpoint_coeffs();
//#endif
//
//
//          log_prob -= MR::Math::pow2(dist_from_sphere1) * scale;
//          log_prob -= MR::Math::pow2(dist_from_sphere2) * scale;
//
//#ifndef GRADIENT_NOT_REQUIRED
//          for (size_t degree_i = 0; degree_i < strand.degree(); degree_i++) {
//            gradient[degree_i] -= (2.0 * scale * dist_from_sphere1 * coeffs.first[degree_i]/ norm1) * endpoints.first;
//            gradient[degree_i] -= (2.0 * scale * dist_from_sphere2 * coeffs.second[degree_i]/ norm2) * endpoints.second;
//          }
//
//
//          for (size_t degree_i1 = 0; degree_i1 < strand.degree(); degree_i1++) {
//            for (size_t dim_i1 = 0; dim_i1 < 3; dim_i1++) {
//              for (size_t degree_i2 = 0; degree_i2 < strand.degree(); degree_i2++) {
//                for (size_t dim_i2 = 0; dim_i2 < 3; dim_i2++) {
//
//                  double& fish_elem = fisher[degree_i1][dim_i1][degree_i2][dim_i2];
//
//                  fish_elem +=  2.0 * scale * coeffs.first[degree_i1] * coeffs.first[degree_i2] *
//                                  endpoints.first[dim_i1] * endpoints.first[dim_i2] / MR::Math::pow2(norm1);
//
//                  fish_elem +=  2.0 * scale * coeffs.second[degree_i1] * coeffs.second[degree_i2] *
//                                  endpoints.second[dim_i1] * endpoints.second[dim_i2] / MR::Math::pow2(norm2);
//
//                  for (size_t degree_i3 = 0; degree_i3 < strand.degree(); degree_i3++) {
//                    for (size_t dim_i3 = 0; dim_i3 < 3; dim_i3++) {
//                      double& fish_grad_elem = fisher_gradient[degree_i3 * 3 + dim_i3 + has_density][degree_i1][dim_i1][degree_i2][dim_i2];
//
//                      fish_grad_elem -= 4.0 * scale * coeffs.first[degree_i1] * coeffs.first[degree_i2] * coeffs.first[degree_i3]
//                                        * endpoints.first[dim_i1] * endpoints.first[dim_i2] * endpoints.first[dim_i3] / MR::Math::pow4(norm1);
//
//                      fish_grad_elem -= 4.0 * scale * coeffs.second[degree_i1] * coeffs.second[degree_i2] * coeffs.second[degree_i3]
//                                        * endpoints.second[dim_i1] * endpoints.second[dim_i2] * endpoints.second[dim_i3] / MR::Math::pow4(norm2);
//
//                      if (dim_i3 == dim_i1) {
//                        fish_grad_elem +=  2.0 * scale * coeffs.first[degree_i1] * coeffs.first[degree_i2] * coeffs.first[degree_i3] *
//                                                      endpoints.first[dim_i2] / MR::Math::pow2(norm1);
//
//                        fish_grad_elem +=  2.0 * scale * coeffs.second[degree_i1] * coeffs.second[degree_i2] * coeffs.second[degree_i3] *
//                                        endpoints.second[dim_i2] / MR::Math::pow2(norm2);
//                      }
//
//                      if (dim_i3 == dim_i2) {
//                        fish_grad_elem +=  2.0 * scale * coeffs.first[degree_i1] * coeffs.first[degree_i2] * coeffs.first[degree_i3] *
//                                                      endpoints.first[dim_i1] / MR::Math::pow2(norm1);
//
//                        fish_grad_elem +=  2.0 * scale * coeffs.second[degree_i1] * coeffs.second[degree_i2] * coeffs.second[degree_i3] *
//                                        endpoints.second[dim_i1] / MR::Math::pow2(norm2);
//                      }
//
//                    }
//                  }
//
//                }
//              }
//            }
//          }
//
//#endif
//
//          return log_prob;
//
//
//        }

    }

  }

}
