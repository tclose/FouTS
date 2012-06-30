/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, 26/01/2011.

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

#ifndef __bts_prob_likelihood_imagediff_rician_h__
#define __bts_prob_likelihood_imagediff_rician_h__

#include "math/math.h"
#include "math/bessel.h"

#include "bts/prob/likelihood.h"

namespace BTS {

	namespace Prob {


	  class Likelihood::Rician : public Likelihood {


      public:

        Rician( const Image::Observed::Buffer& obs_image,
                Image::Expected::Buffer* const exp_image,
                double assumed_snr = Likelihood::SNR_DEFAULT,
                const std::string& b0_include = Likelihood::B0_INCLUDE_DEFAULT,
                double outside_scale = Likelihood::OUTSIDE_SCALE_DEFAULT,
                const std::string& ref_b0 = Likelihood::REF_B0_DEFAULT,
                double ref_signal = NAN)

          : Likelihood(obs_image, exp_image, assumed_snr, b0_include, outside_scale, ref_b0, ref_signal) {}


        Rician(const Rician& s)
          : Likelihood(s) {}


        ~Rician() {}

        Rician&                         operator= (const Rician& s)
          { this->Likelihood::operator=(s); return *this; }


        using Likelihood::log_prob;

        double                          log_prob_and_fisher( const Fibre::Strand::Set& strands,
                                                  Fibre::Strand::Set& gradient,
                                                  Fibre::Strand::Set::Tensor& fisher_info,
                                                  std::vector< MR::Math::Matrix<double> >& fisher_info_gradients)
          { return log_prob_and_fisher_tpl<Fibre::Strand>(strands, gradient, fisher_info, fisher_info_gradients); }

        double                          log_prob_and_fisher( const Fibre::Tractlet::Set& strands,
                                                  Fibre::Tractlet::Set& gradient,
                                                  Fibre::Tractlet::Set::Tensor& fisher_info,
                                                  std::vector< MR::Math::Matrix<double> >& fisher_info_gradients)
          { return log_prob_and_fisher_tpl<Fibre::Tractlet>(strands, gradient, fisher_info, fisher_info_gradients); }


        template <typename T> double     log_prob_and_fisher_tpl(const typename T::Set& fibres,
                                                      typename T::Set& gradient,
                                                      typename T::Set::Tensor& fisher_info,
                                                      std::vector< MR::Math::Matrix<double> >& fisher_info_gradients)
          {  throw Exception ("Not implemented yet."); }


        double                          log_prob_and_fisher( const Fibre::Strand::Set& strands,
                                                  Fibre::Strand::Set& gradient,
                                                  Fibre::Strand::Set::Tensor& fisher_info)
          { return log_prob_and_fisher_tpl<Fibre::Strand>(strands, gradient, fisher_info); }

        double                          log_prob_and_fisher( const Fibre::Tractlet::Set& strands,
                                                  Fibre::Tractlet::Set& gradient,
                                                  Fibre::Tractlet::Set::Tensor& fisher_info)
          { return log_prob_and_fisher_tpl<Fibre::Tractlet>(strands, gradient, fisher_info); }


        template <typename T> double     log_prob_and_fisher_tpl(const typename T::Set& fibres,
                                                      typename T::Set& gradient,
                                                      typename T::Set::Tensor& fisher_info)
          { throw Exception ("Not implemented yet."); }


        double                           log_prob(double expected, double observed) {

          assert(expected > 0.0 && observed > 0.0);

          double bessel_x = expected * observed / this->sigma2;

          double lprob =  MR::Math::log(observed / this->sigma2 *
                          MR::Math::exp (- (MR::Math::pow2(observed) + MR::Math::pow2(expected))/ (2.0 * this->sigma2)) *
                          MR::Math::Bessel::I0_scaled(bessel_x)) + bessel_x;

          return lprob;

        }


        double                           b0_log_prob(double expected, double observed) {

          assert(expected > 0.0 && observed > 0.0);

          double lprob;

          if (b0_include == "full") {

            double bessel_x = expected * observed / this->sigma2;

            lprob =  MR::Math::log(observed / this->sigma2 *
                            MR::Math::exp (- (MR::Math::pow2(observed) + MR::Math::pow2(expected))/ (2.0 * this->sigma2)) *
                            MR::Math::Bessel::I0_scaled(bessel_x)) + bessel_x;

          } else
            lprob = 0.0;

          return lprob;

        }


        double                           log_prob(double expected, double observed, double& d_lprob) {

          throw Exception ("Not implemented yet.");

        }

        double                           log_prob(double expected, double observed, double& d_lprob, double& d2_lprob2) {

          throw Exception ("Not implemented yet.");
        }



	  };

	}

}

#endif /* __bts_prob_likelihood_imagediff_rician_h__ */
