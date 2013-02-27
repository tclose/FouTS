/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Created by Tom Close on 13/03/09.

 This file is part of Bayesian Tractlet Sampling (BTS).

 BTS is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 BTS is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with BTS.  If not, see <http://www.gnu.org/licenses/>.

 */



#include "dwi/gradient.h"
#include "math/matrix.h"

#include "math/SH.h"

#include "bts/common.h"

#include "bts/triple.h"
//

#include "bts/diffusion/encoding.h"

#include "bts/diffusion/model.h"
#include "bts/diffusion/response.h"

#include "bts/diffusion/inline_functions.h"

namespace BTS {

  namespace Diffusion {

    const std::string                   Model::ENCODINGS_LOCATION_DEFAULT   = str(getenv("HOME")) + str("/git/FouTS/params/diffusion/encoding_60.b");
    const bool                          Model::ISOTROPIC_DEFAULT            = 0;
    const double                        Model::ADC_DEFAULT                  = 5e-4;
    const double                        Model::FA_DEFAULT                   = 0.8;

    // Maximum order ('l') spherical harmonic that can be calculated;
    const size_t                        Model::LMAX        = 8;

    // 0-8th order m=0 Associated Legendre polynomial coefficients multiplied by the common factor for each order.
    // (Note, The common factor is omitted in the constant for accuracy and is calculated from the expression
    // sqrt[(2*l+1)/(4*pi)] / 2^l). Higher orders can be added, just need to track down the right Legendre coefficients.
    const double                        Model::M0_HARMONIC_L0[5] =  {  1,     0,      0,        0,     0};
    const double                        Model::M0_HARMONIC_L2[5] =  { -2,     6,      0,        0,     0};
    const double                        Model::M0_HARMONIC_L4[5] =  {  6,   -60,     70,        0,     0};
    const double                        Model::M0_HARMONIC_L6[5] =  {-20,   420,  -1260,      924,     0};
    const double                        Model::M0_HARMONIC_L8[5] =  { 70, -2320,  13860,   -24024, 12870};




    Model                               Model::factory(const MR::Math::Matrix<double>& encodings_matrix,
                                                      const MR::Math::Matrix<double>& response_SHs,
                                                      double adc,
                                                      double fa,
                                                      bool include_isotropic,
                                                      bool warn_on_b_mismatch) {

      Model model;

      //If response_SHs is an empty matrix, use adc and fa values to generate response spherical harmonics from tensor model.
      if (!response_SHs.rows() && !response_SHs.columns())
        model = Model(encodings_matrix, adc, fa, include_isotropic);

      //Else if response matrix only has one column assume that it meant to be a column vector.
      else if (response_SHs.columns() == 1)
        model = Model(encodings_matrix, response_SHs.column(0), include_isotropic, warn_on_b_mismatch);

      //Else if response matrix only has one row assume that it meant to be a row vector.
      else if (response_SHs.rows() == 1)
        model = Model(encodings_matrix, response_SHs.row(0), include_isotropic, warn_on_b_mismatch);

      //Otherwise assume that it has separate response coefficients on each row.
      else
        model = Model(encodings_matrix, response_SHs, include_isotropic);

      return model;

    }


    Model::Model (const MR::Math::Matrix<double>& encodings_matrix, const MR::Math::Matrix<double>& response_SHs,
                                                                                               bool include_isotropic) {

      if (response_SHs.rows() != encodings_matrix.rows())
        throw Exception ("If multiple rows to thre response_SHs spherical harmonics are provided then they must match" +
                                                                   str(" the number of rows in the encodings matrix."));

      init(encodings_matrix, response_SHs, include_isotropic);

    }




    void                                Model::precalculate_weightings(Fibre::Strand::Section& section) const {

      section.precalc_weightings.resize(num_encodings());

      for (size_t response_i = 0; response_i < num_encodings(); ++response_i)
        section.precalc_weightings[response_i] = operator[](response_i).weighting(section.tangent());

    }

    void                                Model::precalculate_weightings_and_gradients(Fibre::Strand::Section& section) const {

      section.precalc_weightings.resize(num_encodings());
      section.precalc_weight_gradients.resize(num_encodings());

      for (size_t response_i = 0; response_i < num_encodings(); ++response_i)
        section.precalc_weightings[response_i] = operator[](response_i).weighting(section.tangent(), section.precalc_weight_gradients[response_i]);


    }

    void                                Model::precalculate_weightings_gradients_and_hessians(Fibre::Strand::Section& section) const {

      section.precalc_weightings.resize(num_encodings());
      section.precalc_weight_gradients.resize(num_encodings());
      section.precalc_weight_hessians.resize(num_encodings());

      for (size_t response_i = 0; response_i < num_encodings(); ++response_i)
        section.precalc_weightings[response_i] = operator[](response_i).weighting(section.tangent(), section.precalc_weight_gradients[response_i], section.precalc_weight_hessians[response_i]);

    }



    MR::Math::Vector<double>            Model::tensor_m0_SH(double adc, double fa, double b_value) {

      MR::Math::Vector<double> sh (LMAX/2 + 1);

      MR::Math::SH::FA2SH (sh, fa, adc, b_value, LMAX);

      return sh;

    }


    MR::Math::Vector<double>            Model::SH_to_coeffs(MR::Math::Vector<double> loaded_SH, bool include_isotropic) {

      size_t num_coeffs = LMAX/2 + 1;

      size_t num_loaded_harmonic_coeffs = loaded_SH.size();

      if (num_loaded_harmonic_coeffs > num_coeffs)
        throw Exception ("Can only calculate polynomial coefficients up to spherical harmonic order " + str((num_coeffs-1)*2) + ".");

      MR::Math::Vector<double> SH (loaded_SH);

      SH.resize(num_coeffs);

      // Pad out the remaining coefficients with zeros.
      for (size_t sh_i = num_loaded_harmonic_coeffs; sh_i < num_coeffs; sh_i++)
        SH[sh_i] = 0.0;

      //Pre-weight the spherical harmonics by the common factor for each polynomial.
      for (size_t sh_i = 0; sh_i < num_coeffs; sh_i++) {

        double l = (double)sh_i * 2.0;

        //Normalising scalar for each basis: sqrt[(2*l+1)/(4*pi)] / 2^l
        SH[sh_i] *= sqrt((2.0 * l + 1.0) / (4.0 * M_PI)) / pow(2.0, l);
      }


      //Set up m=0 associated Legendre polynomials (spherical harmonic basis functions)
      std::vector< std::vector<double> > uncollated;

      if (include_isotropic)
        uncollated.push_back(std::vector<double>(M0_HARMONIC_L0, M0_HARMONIC_L0 + num_coeffs));
      else
        uncollated.push_back(std::vector<double>(num_coeffs, 0.0));


      uncollated.push_back(std::vector<double>(M0_HARMONIC_L2, M0_HARMONIC_L2 + num_coeffs));
      uncollated.push_back(std::vector<double>(M0_HARMONIC_L4, M0_HARMONIC_L4 + num_coeffs));
      uncollated.push_back(std::vector<double>(M0_HARMONIC_L6, M0_HARMONIC_L6 + num_coeffs));
      uncollated.push_back(std::vector<double>(M0_HARMONIC_L8, M0_HARMONIC_L8 + num_coeffs));


      // Multiply spherical harmonic basis functions with weighted coefficients.
      for (size_t sh_i = 0; sh_i < num_coeffs; sh_i++)
        for (size_t poly_i = 0; poly_i < num_coeffs; poly_i++)
          uncollated[sh_i][poly_i] *= SH[sh_i];


      //Collect like terms between the polynomial equations to form one set of polynomial coefficients.
      MR::Math::Vector<double> collated(num_coeffs);
      collated.zero();

      for (size_t poly_i = 0; poly_i < num_coeffs; poly_i++)
        for (size_t sh_i = 0; sh_i < num_coeffs; sh_i++)
          collated[poly_i] += uncollated[sh_i][poly_i];

      return collated;

    }


    Model::Model (const MR::Math::Matrix<double>& encodings_matrix, double adc, double fa, bool include_isotropic) {

      MR::Math::Matrix<double> response_SHs(encodings_matrix.rows(), LMAX/2 + 1);

      for (size_t row_i = 0; row_i < encodings_matrix.rows(); row_i++)
        response_SHs.row(row_i) = tensor_m0_SH(adc, fa, encodings_matrix(row_i,3));

      init(encodings_matrix, response_SHs, include_isotropic);

    }




    Model::Model (const MR::Math::Matrix<double>& encodings_matrix, const MR::Math::Vector<double>& response_SH,
                  bool include_isotropic, bool warn_on_b_mismatch) {

      if (!encodings_matrix.rows())
        throw Exception ("No rows found in encodings matrix.");

      if (encodings_matrix.columns() != 4)
        throw Exception ("Four columns required in encodings matrix, found " + str(encodings_matrix.columns()) + ".");


      double consistent_b_value = 0.0;

      MR::Math::Matrix<double> response_SHs(encodings_matrix.rows(), LMAX/2 + 1);

      for (size_t row_i = 0; row_i < encodings_matrix.rows(); row_i++) {

        //If b value == previously set b value or the previous b value hasn't been set yet.
        if (encodings_matrix(row_i,3) == consistent_b_value || consistent_b_value == 0) {
          response_SHs.row(row_i) = response_SH;
          consistent_b_value = encodings_matrix(row_i,3);
        //If b value == 0. Set response coeffs to NAN, as they should be ignored in 'init' method.
        } else if (encodings_matrix(row_i,3) == 0)
          response_SHs.row(row_i) = NAN;
        else {
          if (warn_on_b_mismatch) {
            std::cout << "Inconsistent b values used in encoding matrix. Previously found b value " <<
                      consistent_b_value <<  ", found " << encodings_matrix(row_i,3) << " at row " << str(row_i) + "." <<
                      std::endl;
            response_SHs.row(row_i) = response_SH;
          } else
            throw Exception("Inconsistent b values used in encoding matrix, please either supply a seperate " +
                            str("response harmonics for each row or use the automatic tensor harmonic calculation. ") +
                            "Previously found b value " + str(consistent_b_value) +  ", found "
                            + str(encodings_matrix(row_i,3)) + " at row " + str(row_i) + ".");
        }

      }


      init(encodings_matrix, response_SHs, include_isotropic);


    }


    MR::Math::Vector<double>            Model::weightings(MR::Math::Vector<double> weightings,
                                                                    MR::Math::Matrix<double> orientations) {

      //TODO: Implement this function
      return weightings;
    }


    void                                Model::init (const MR::Math::Matrix<double>& encodings_matrix, const MR::Math::Matrix<double>& response_SHs, bool include_isotropic) {
            
      if (encodings_matrix.columns() <= 3)
        throw Exception ("Insufficient columns in loaded diffusion-weighting matrix");
      
      if (encodings_matrix.columns() > 4)
       throw Exception ("More than 4 columns (" + str(encodings_matrix.columns()) + ") were found in diffusion-weighting matrix");

      if (encodings_matrix.rows() != response_SHs.rows())
        throw Exception ("Number of rows in encodings_matrix (" + str(encodings_matrix.rows()) + ") and response_SH matrix (" + str(response_SHs.rows()) + ") do not match.");

      for (size_t row_i = 0; row_i < encodings_matrix.rows(); row_i++) {

        MR::Math::Vector<double> response_coeffs(1);

        if (encodings_matrix(row_i, 3))
          response_coeffs = SH_to_coeffs(response_SHs.row(row_i), include_isotropic);

        //If b value == 0, the 0th order coefficient is set to one and the remaining are 0.
        else
          response_coeffs[0] = 1.0;

        Response response(encodings_matrix.row(row_i), response_coeffs, row_i);

        responses.push_back(response);

      }

      this->includes_iso = include_isotropic;


    }


      

    void                                Model::scale_coeffs(double scalar) {
    
      for (size_t encode_i = 0; encode_i < this->num_encodings(); encode_i++)
        this->operator[](encode_i).scale_coeffs(scalar);
        
    }
      



    std::ostream& operator<< (std::ostream& stream, const Diffusion::Model& model) {

      stream << std::endl;

      for (std::vector<Response>::const_iterator response_it = model.responses.begin(); response_it != model.responses.end(); ++response_it)
        stream << *response_it << std::endl;

      return (stream);

    }



  }


}
