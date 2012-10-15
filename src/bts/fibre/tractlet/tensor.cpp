/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, 24/11/2010.

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

#include "bts/fibre/tractlet.h"
#include "bts/fibre/tractlet/tensor.h"
#include "bts/fibre/tractlet/section.h"
#include "bts/fibre/tractlet/section/tensor.h"


namespace BTS {

	namespace Fibre {

    void      Fibre::Tractlet::Tensor::add_section_hessian(const Tractlet& tractlet,
                                                                      const Tractlet::Section& section,
                                                                      const Strand::BasicSection& gradient,
                                                                      const Strand::BasicSection::Tensor& hessian) {


      Fibre::Tractlet d2_intensity = get_template();

      Strand pos_gradient = Strand::outer_product(section.position_coeffs, gradient.position() / section.intensity());

      if (has_var_acs())
        d2_intensity.set_acs(0.0);

//      d2_intensity.set_base_width(0.0);

      //Add position gradient.
      d2_intensity[0] += pos_gradient;
      d2_intensity[1] += pos_gradient * section.ax1_fraction;
      d2_intensity[2] += pos_gradient * section.ax2_fraction;

      //Add tangent gradient.

      Strand tang_gradient = Strand::outer_product(section.tangent_coeffs,  gradient.tangent() * (section.length_fraction / section.intensity()));

      d2_intensity[0] += tang_gradient;
      d2_intensity[1] += tang_gradient * section.ax1_fraction;
      d2_intensity[2] += tang_gradient * section.ax2_fraction;


      Fibre::Tractlet acs_gradient = d2_intensity * section.intensity();


      if (true) { //Needs some serious rewriting.

        throw Exception ("tied_width hessian needs to be adjusted after accounting for tractlet sheer.");

        double norm1 = tractlet[1][0].norm();
        double norm2 = tractlet[2][0].norm();

        Coord d_intensity_d1 = tractlet.acs() * ((norm2/norm1) * tractlet(1,0) - tractlet(2,0));
        Coord d_intensity_d2 = tractlet.acs() * ((norm1/norm2) * tractlet(2,0) - tractlet(1,0));

        for (size_t dim_i1 = 0; dim_i1 < 3; dim_i1++) {


          operator()(1,0,dim_i1,2,0,dim_i1) -= gradient.intensity() * tractlet.acs();
          operator()(2,0,dim_i1,1,0,dim_i1) -= gradient.intensity() * tractlet.acs();

          operator()(1,0,dim_i1,1,0,dim_i1) += gradient.intensity() * tractlet.acs() * norm2 / norm1;
          operator()(2,0,dim_i1,2,0,dim_i1) += gradient.intensity() * tractlet.acs() * norm1 / norm2;


          for (size_t dim_i2 = 0; dim_i2 < 3; dim_i2++) {


            operator()(1,0,dim_i1,1,0,dim_i2) -= tractlet.acs() * gradient.intensity() * tractlet[1][0][dim_i1] * tractlet[1][0][dim_i2] * (norm2 / MR::Math::pow3(norm1));
            operator()(2,0,dim_i1,2,0,dim_i2) -= tractlet.acs() * gradient.intensity() * tractlet[2][0][dim_i1] * tractlet[2][0][dim_i2] * (norm1 / MR::Math::pow3(norm2));

            operator()(1,0,dim_i1,2,0,dim_i2) += tractlet.acs() * gradient.intensity() * tractlet[2][0][dim_i2] * tractlet[1][0][dim_i1] / (norm1 * norm2);
            operator()(2,0,dim_i1,1,0,dim_i2) += tractlet.acs() * gradient.intensity() * tractlet[1][0][dim_i2] * tractlet[2][0][dim_i1] / (norm2 * norm1);

            for (size_t ax_i = 0; ax_i < 3; ++ax_i)
              for (size_t degree_i = 0; degree_i < tractlet.degree(); ++degree_i) {

                operator()(ax_i,degree_i,dim_i1,1,0,dim_i2) += tractlet.acs() * d2_intensity(ax_i,degree_i)[dim_i1] * d_intensity_d1[dim_i2];
                operator()(ax_i,degree_i,dim_i1,2,0,dim_i2) += tractlet.acs() * d2_intensity(ax_i,degree_i)[dim_i1] * d_intensity_d2[dim_i2];

                operator()(1,0,dim_i1,ax_i,degree_i,dim_i2) += tractlet.acs() * d2_intensity(ax_i,degree_i)[dim_i2] * d_intensity_d1[dim_i1];
                operator()(2,0,dim_i1,ax_i,degree_i,dim_i2) += tractlet.acs() * d2_intensity(ax_i,degree_i)[dim_i2] * d_intensity_d2[dim_i1];

              }

          }

        }


        if (has_var_acs()) {
          acs_gradient(1,0) += gradient.intensity() * d_intensity_d1;
          acs_gradient(2,0) += gradient.intensity() * d_intensity_d2;
        }

      }


      if (has_var_acs())
        acs() += acs_gradient;


      for (size_t degree_i1 = 0; degree_i1 < degree(); degree_i1++)
        for (size_t dim_i1 = 0; dim_i1 < 3; dim_i1++) {

          Fibre::Tractlet ax0 = row(0,degree_i1,dim_i1);
          Fibre::Tractlet ax1 = row(1,degree_i1,dim_i1);
          Fibre::Tractlet ax2 = row(2,degree_i1,dim_i1);

          if (has_var_acs()) {
            ax0.alpha() += acs_gradient.operator()(0,degree_i1)[dim_i1];
            ax1.alpha() += acs_gradient.operator()(1,degree_i1)[dim_i1];
            ax2.alpha() += acs_gradient.operator()(2,degree_i1)[dim_i1];
          }


          for (size_t degree_i2 = 0; degree_i2 < degree(); degree_i2++)
            for (size_t dim_i2 = 0; dim_i2 < 3; dim_i2++) {

              double ax0_ax0 = section.position_coeffs[degree_i1] * section.position_coeffs[degree_i2] * hessian.position(dim_i1).position(dim_i2)
                              + section.position_coeffs[degree_i1] * section.tangent_coeffs[degree_i2]  * hessian.position(dim_i1).tangent(dim_i2) * section.length_fraction
                              + section.tangent_coeffs[degree_i1]  * section.position_coeffs[degree_i2] * hessian.tangent(dim_i1).position(dim_i2) * section.length_fraction
                              + section.tangent_coeffs[degree_i1]  * section.tangent_coeffs[degree_i2]  * hessian.tangent(dim_i1).tangent(dim_i2) * MR::Math::pow2(section.length_fraction);

              ax0[0][degree_i2][dim_i2] += ax0_ax0;
              ax0[1][degree_i2][dim_i2] += ax0_ax0 * section.ax1_fraction;
              ax1[0][degree_i2][dim_i2] += ax0_ax0 * section.ax1_fraction;
              ax0[2][degree_i2][dim_i2] += ax0_ax0 * section.ax2_fraction;
              ax2[0][degree_i2][dim_i2] += ax0_ax0 * section.ax2_fraction;

              ax1[1][degree_i2][dim_i2] += ax0_ax0 * MR::Math::pow2(section.ax1_fraction);
              ax1[2][degree_i2][dim_i2] += ax0_ax0 * section.ax1_fraction * section.ax2_fraction;
              ax2[1][degree_i2][dim_i2] += ax0_ax0 * section.ax1_fraction * section.ax2_fraction;

              ax2[2][degree_i2][dim_i2] += ax0_ax0 * MR::Math::pow2(section.ax2_fraction);

            }


        }

    }

	}

}
