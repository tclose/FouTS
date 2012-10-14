/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, 17/06/2010.

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

#ifndef __bts_image_expected_direction_cpp_h__
#define __bts_image_expected_direction_cpp_h__

#include "bts/image/expected/voxel.h"

namespace BTS {

  namespace Image {


    template <typename T> double      Expected::Direction::signal_tpl (const T& section) {

      //OPTIMISED makes use of the fact that each section is processed sequentially, allowing the weighting and interpolation to be
      //precalculated and stored withinness the Diffusion::Response and Fabricated::Voxel objects respectively. Then quickly combined
      //for all of the diffusion-weighting/voxel combinations.
#ifdef OPTIMISED

      double diffusion_weighting       = section.precalc_weightings[response->get_index()];
      double interpolation_weighting   = section.precalc_interpolation;

#else //ifndef OPTIMISED

      double diffusion_weighting     = response->weighting (section.tangent());
      double interpolation_weighting = voxel->interpolate(section);

#endif


      return section.intensity() * diffusion_weighting * interpolation_weighting;

    }

    template <typename T> double      Expected::Direction::signal_tpl (const T& section, Fibre::Strand::BasicSection& gradient) {

      //NB: diffusion gradient refers to the gradient of the section tangent
      //w.r.t. the diffusion weighting, as opposed to the direction of the
      //magnetic field "gradient".


      //OPTIMISED makes use of the fact that each section is processed sequentially, allowing the weighting and interpolation to be
      //precalculated and stored withinness the Diffusion::Response and Fabricated::Voxel objects respectively. Then quickly combined
      //for all of the diffusion-weighting/voxel combinations.
#ifdef OPTIMISED

      double diffusion_weighting       = section.precalc_weightings[response->get_index()];
      double interpolation_weighting   = section.precalc_interpolation;

#ifndef GRADIENT_NOT_REQUIRED
      const Triple<double>& diffusion_gradient  = section.precalc_weight_gradients[response->get_index()];
      const Fibre::Strand::BasicSection& interpolation_gradient           = section.precalc_interp_gradient;
#endif

#else //ifndef OPTIMISED

      Coord diffusion_gradient;
      T interpolation_gradient;

      double diffusion_weighting     = response->weighting (section.tangent(), diffusion_gradient);
      double interpolation_weighting = voxel->interpolate(section, interpolation_gradient);

#endif


#ifndef GRADIENT_NOT_REQUIRED
      signal_gradient(section, gradient, diffusion_weighting, interpolation_weighting, diffusion_gradient, interpolation_gradient);
#endif

      return section.intensity() * diffusion_weighting * interpolation_weighting;

    }




    template <typename T> double      Expected::Direction::signal_tpl (const T& section, Fibre::Strand::BasicSection& gradient, typename T::Tensor& hessian) {


      //NB: diffusion gradient refers to the gradient of the section tangent
      //w.r.t. the diffusion weighting, as opposed to the direction of the
      //magnetic field "gradient".


      //OPTIMISED makes use of the fact that each section is processed sequentially, allowing the weighting and interpolation to be
      //precalculated and stored withinness the Diffusion::Response and Fabricated::Voxel objects respectively. Then quickly combined
      //for all of the diffusion-weighting/voxel combinations.
#ifdef OPTIMISED

      double diffusion_weighting       = section.precalc_weightings[response->get_index()];
      double interpolation_weighting   = section.precalc_interpolation;

#ifndef GRADIENT_NOT_REQUIRED
      const Triple<double>& diffusion_gradient = section.precalc_weight_gradients[response->get_index()];
      const Fibre::Strand::BasicSection& interpolation_gradient           = section.precalc_interp_gradient;

#ifndef HESSIAN_NOT_REQUIRED
      const Coord::Tensor& diffusion_hessian   = section.precalc_weight_hessians[response->get_index()];
      const typename T::Tensor& interpolation_hessian   = section.precalc_interp_hessian;
#endif

#endif

#else //ifndef OPTIMISED


#ifndef GRADIENT_NOT_REQUIRED

      Coord diffusion_gradient;
      T interpolation_gradient;

#ifndef HESSIAN_NOT_REQUIRED

      typename T::Tensor interpolation_hessian;
      Coord::Tensor diffusion_hessian;

      double diffusion_weighting     = response->weighting (section.tangent(), diffusion_gradient, diffusion_hessian);
      double interpolation_weighting = voxel->interpolate(section, interpolation_gradient, interpolation_hessian);

#else

      double diffusion_weighting     = response->weighting (section.tangent, diffusion_gradient);
      double interpolation_weighting = voxel->interpolate(section, interpolation_gradient);

#endif

#else

      double diffusion_weighting     = response->weighting(section.tangent);
      double interpolation_weighting = voxel->interpolate(section);

#endif

#endif


#ifndef GRADIENT_NOT_REQUIRED
      signal_gradient(section, gradient, diffusion_weighting, interpolation_weighting, diffusion_gradient, interpolation_gradient);

#ifndef HESSIAN_NOT_REQUIRED
      signal_hessian(section, hessian, diffusion_weighting, interpolation_weighting, diffusion_gradient, interpolation_gradient, diffusion_hessian, interpolation_hessian);
#endif

#endif

      return section.intensity() * diffusion_weighting * interpolation_weighting;


    }



    inline void             Expected::Direction::signal_gradient(const Fibre::Strand::Section& section,
                                                    Fibre::Strand::BasicSection& gradient,
                                                    double diffusion_weighting,
                                                    double interpolation_weighting,
                                                    const Triple<double>& diffusion_gradient,
                                                    const Fibre::Strand::BasicSection& interpolation_gradient) {

      gradient.tangent() = (Coord)diffusion_gradient * (section.intensity() * interpolation_weighting) + interpolation_gradient.tangent() * (section.intensity() * diffusion_weighting);

      gradient.position() = interpolation_gradient.position() * (section.intensity() * diffusion_weighting);

      gradient.intensity() = diffusion_weighting * interpolation_weighting;


    }

    inline void              Expected::Direction::signal_gradient(const Fibre::Tractlet::Section& section,
                                                    Fibre::Strand::BasicSection& gradient,
                                                    double diffusion_weighting,
                                                    double interpolation_weighting,
                                                    const Triple<double>& diffusion_gradient,
                                                    const Fibre::Strand::BasicSection& interpolation_gradient) {


      //Call the Fibre::Strand::Section version of the function first, then process the additional parameters.
      signal_gradient ((const Fibre::Strand::Section&)section, (Fibre::Strand::BasicSection&)gradient, diffusion_weighting, interpolation_weighting, diffusion_gradient, (const Fibre::Strand::BasicSection&)interpolation_gradient);

//      //Check to see if the tractlet section 'width' is used in the interpolation calculation.
//      if (!interpolation_gradient.width1 || !interpolation_gradient.width2) {
//        gradient.width1 = Triple<double>::Zeros;
//        gradient.width2 = Triple<double>::Zeros;
//      } else {
//        gradient.width1 = interpolation_gradient.width1 * (section.intensity() * diffusion_weighting);
//        gradient.width2 = interpolation_gradient.width2 * (section.intensity() * diffusion_weighting);
//      }

    }


    inline void             Expected::Direction::signal_hessian(const Fibre::Strand::Section& section,
                                                    Fibre::Strand::BasicSection::Tensor& hessian,
                                                    double diff_weighting,
                                                    double interp_weighting,
                                                    const Triple<double>& diff_gradient,
                                                    const Fibre::Strand::BasicSection& interp_gradient,
                                                    const Coord::Tensor& diff_hessian,
                                                    const Fibre::Strand::BasicSection::Tensor& interp_hessian) {

      //Set hessian intensity elements.

      hessian.intensity().intensity() = 0.0;

      hessian.intensity().tangent() = (Coord)diff_gradient * interp_weighting + interp_gradient.tangent() * diff_weighting;
      hessian.intensity().position() = interp_gradient.position() * diff_weighting;

      //TODO: reformulate this code in terms of MR::Math::Matrix<double>::sub
      hessian.position(X).intensity() = hessian.intensity().position(X);
      hessian.position(Y).intensity() = hessian.intensity().position(Y);
      hessian.position(Z).intensity() = hessian.intensity().position(Z);

      hessian.tangent(X).intensity() = hessian.intensity().tangent(X);
      hessian.tangent(Y).intensity() = hessian.intensity().tangent(Y);
      hessian.tangent(Z).intensity() = hessian.intensity().tangent(Z);


      //Set hessian position()-position() elements.

      hessian.position(X).position() = interp_hessian.position(X).position() * diff_weighting * section.intensity();
      hessian.position(Y).position() = interp_hessian.position(Y).position() * diff_weighting * section.intensity();
      hessian.position(Z).position() = interp_hessian.position(Z).position() * diff_weighting * section.intensity();


      //Set hessian tangent()-tangent() elements.

      Coord::Tensor tangent_hessian = section.intensity() * (
                                                                   diff_hessian * interp_weighting // FIXME: add a operator*(double) function for Coord::Tensor
                                                                 + interp_gradient.tangent().outer(diff_gradient)
                                                                 + ((Coord)diff_gradient).outer(interp_gradient.tangent())
                                                                 + interp_hessian.tangent() * diff_weighting); // FIXME: Make a tangent() function of BasicSection::Tensor

      hessian.tangent(X).tangent() = tangent_hessian[X];
      hessian.tangent(Y).tangent() = tangent_hessian[Y];
      hessian.tangent(Z).tangent() = tangent_hessian[Z];


      //Set hessian position()-tangent() and vice-versa elements.

      Coord::Tensor cross_hessian = section.intensity() * (
                                                interp_gradient.position().outer(diff_gradient)
                                                + interp_hessian.position_tangent() * diff_weighting); // FIXME: Make a cross() function of BasicSection::Tensor


      hessian.position(X).tangent() = cross_hessian[X];
      hessian.position(Y).tangent() = cross_hessian[Y];
      hessian.position(Z).tangent() = cross_hessian[Z];

      hessian.tangent(X).position(X) = cross_hessian[X][X];
      hessian.tangent(X).position(Y) = cross_hessian[Y][X];
      hessian.tangent(X).position(Z) = cross_hessian[Z][X];


      hessian.tangent(Y).position(X) = cross_hessian[X][Y];
      hessian.tangent(Y).position(Y) = cross_hessian[Y][Y];
      hessian.tangent(Y).position(Z) = cross_hessian[Z][Y];

      hessian.tangent(Z).position(X) = cross_hessian[X][Z];
      hessian.tangent(Z).position(Y) = cross_hessian[Y][Z];
      hessian.tangent(Z).position(Z) = cross_hessian[Z][Z];

    }



    inline void              Expected::Direction::signal_hessian(const Fibre::Tractlet::Section& section,
                                                    Fibre::Strand::BasicSection::Tensor& hessian,
                                                    double diffusion_weighting,
                                                    double interpolation_weighting,
                                                    const Triple<double>& diffusion_gradient,
                                                    const Fibre::Strand::BasicSection& interpolation_gradient,
                                                    const Coord::Tensor& diffusion_hessian,
                                                    const Fibre::Strand::BasicSection::Tensor& interpolation_hessian) {


      //Call the Fibre::Strand::Section version of the function first, then process the additional parameters.
      signal_hessian ( (const Fibre::Strand::Section&)section,
                        (Fibre::Strand::BasicSection::Tensor&)hessian,
                        diffusion_weighting,
                        interpolation_weighting,
                        diffusion_gradient,
                        (const Fibre::Strand::Section&)interpolation_gradient,
                        diffusion_hessian,
                        (const Fibre::Strand::BasicSection::Tensor&)interpolation_hessian);

      //Check to see if the tractlet section 'width' is used in the interpolation calculation.
      //TODO: Note Nonzero Tractlet::Section widths are not currently handled in Hessian calculation.

    }



    inline std::ostream&      operator<<(std::ostream& stream, const Expected::Direction& direction) {

      stream << "Voxel: " << direction.voxel << ", Response: " << direction.response << std::endl;

      return stream;

    }

  }

}



#endif /* __bts_image_expected_direction_cpp_h__ */
