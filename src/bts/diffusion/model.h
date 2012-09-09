/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Created by Tom Close on 13/03/09.

 This file is part of Bayesian Tractlet Sampling (BTS).

 BTS is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your diff_option) any later version.

 BTS is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with BTS.  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef __bts_diffusion_model_h__
#define __bts_diffusion_model_h__


//Defines the parameters required to initialise a Proposal::Distribution object.
#define DIFFUSION_PARAMETERS \
\
  Option ("diff_response_SH_location", "The location of the axially symmetric (m=0) spherical harmonic coefficients of the diffusion response function.") \
   + Argument ("diff_response_SH_location", "").type_file (), \
\
  Option ("diff_encodings_location", "The location of the gradient encoding used, supplied as a 4xN text file with each line is in the format [ X Y Z b ], where [ X Y Z ] describe the direction of the applied gradient, and b gives the b-value in units (1000 s/mm^2).") \
   + Argument ("diff_encodings_location", "").type_file (), \
\
  Option ("diff_isotropic", "Include isotropic component of the signal."), \
\
  Option ("diff_adc", "The mean apparent diffusion coefficient (mm^2/s) of the tensor model.") \
   + Argument ("diff_adc", "").type_float (SMALL_FLOAT, Diffusion::Model::ADC_DEFAULT, LARGE_FLOAT), \
\
  Option ("diff_fa", "The fractional anisotropy (FA) of the diffusion tensor model.") \
   + Argument ("diff_fa", "").type_float (0.0, Diffusion::Model::FA_DEFAULT, 1.0) \
\

//  Option ("precision", "azimuthal angle precision", "precision of the samples about the azimuthal angle")
//   + Argument ("precision", "precision of the samples about the azimuthal angle").type_integer (3,1e6,100))


//Loads the 'proposal' parameters into variables
#define SET_DIFFUSION_PARAMETERS \
\
  std::string diff_response_SH_location; \
  std::string diff_encodings_location             = Diffusion::Model::ENCODINGS_LOCATION_DEFAULT; \
  bool diff_isotropic                             = Diffusion::Model::ISOTROPIC_DEFAULT; \
  double diff_adc                                 = Diffusion::Model::ADC_DEFAULT; \
  double diff_fa                                  = Diffusion::Model::FA_DEFAULT; \
\
  Options diff_opt = get_options("diff_response_SH_location"); \
  if (diff_opt.size()) \
    diff_response_SH_location = diff_opt[0][0].c_str(); \
\
  diff_opt = get_options("diff_encodings_location"); \
  if (diff_opt.size()) \
    diff_encodings_location = diff_opt[0][0].c_str(); \
\
  diff_opt = get_options("diff_isotropic"); \
  if (diff_opt.size()) \
    diff_isotropic = !Diffusion::Model::ISOTROPIC_DEFAULT; \
\
  diff_opt = get_options("diff_adc"); \
  if (diff_opt.size()) {\
    diff_adc = Diffusion::Model::ADC_DEFAULT; \
    if (diff_response_SH_location.size()) \
      std::cout << "WARNING!! Supplied value for diff_option '-diff_adc' (" + str(diff_adc) + ") will be ignored as supplied value for diff_option '-diff_response_SH' overrides it." << std::endl; \
  } \
\
  diff_opt = get_options("diff_fa"); \
  if (diff_opt.size()) {\
    diff_fa = Diffusion::Model::FA_DEFAULT; \
    if (diff_response_SH_location.size()) \
      std::cout << "WARNING!! Supplied value for diff_option '-diff_fa' (" + str(diff_fa) + ") will be ignored as supplied value for diff_option '-diff_response_SH' overrides it." << std::endl; \
  } \
\
  MR::Math::Matrix<double> diff_response_SH; \
  if (diff_response_SH_location.size()) { \
    diff_response_SH.load (diff_response_SH_location); \
  } \
\
  MR::Math::Matrix<double> diff_encodings (diff_encodings_location); \
\
  

//           = Diffusion::Model::RESPONSE_SH_LOCATION_DEFAULT;
//  std::string diff_spherical_harmonic_location    = Diffusion::Model::SPHERICAL_HARMONIC_LOCATION_DEFAULT;
//  MR::Math::Vector<double> diff_spherical_harmonic (diff_spherical_harmonic_location);

//Adds the 'signal' parameters to the properties to be saved with the data.
#define ADD_DIFFUSION_PROPERTIES(properties) \
  properties["diff_encodings_location"]     = str(diff_encodings_location); \
  properties["diff_isotropic"]              = str(diff_isotropic); \
  if (diff_response_SH.rows()) { \
    properties["diff_response_SH"]            = Math::matlab_str(diff_response_SH); \
    properties["diff_response_SH_location"]   = diff_response_SH_location; \
  } else { \
    properties["diff_adc"]                    = str(diff_adc); \
    properties["diff_fa"]                     = str(diff_fa); \
  } \




namespace BTS {
  namespace Diffusion {
  
    class Model;
    
  }
}

#include <vector>

#include "bts/fibre/strand/set.h"


#include "bts/diffusion/encoding.h"
#include "bts/diffusion/response.h"
 
#define LOOP(op) \
for (std::vector<Response>::iterator response_it = responses.begin(); response_it != responses.end(); ++response_it) { op }

namespace BTS {

  namespace Diffusion {

    class Model {

      //Public const static members
      public:
        
        const static char*                       ENCODINGS_LOCATION_DEFAULT;

        const static bool                        ISOTROPIC_DEFAULT;

        const static double                      ADC_DEFAULT;
        const static double                      FA_DEFAULT;
//        const static double                      B_VALUE_DEFAULT;

        const static double                      AZ_PRECISION_DEFAULT;
//        const static char*                       SPHERICAL_HARMONIC_LOCATION_DEFAULT;
//        const static char*                       RESPONSE_SH_LOCATION_DEFAULT;

        const static size_t LMAX; //Maximum lmax that can be calculated;

        //0-8th order m=0 Associated Legendre polynomial coefficients multiplied by the common factor for each order.
        static const double M0_HARMONIC_L0[5], M0_HARMONIC_L2[5], M0_HARMONIC_L4[5], M0_HARMONIC_L6[5], M0_HARMONIC_L8[5];


      //Protected member variables
      protected:
        
        std::vector<Response>                   responses;
        bool                                    includes_iso;

      public:


        static MR::Math::Vector<double>         tensor_m0_SH(double adc, double fa, double b_value);

        static MR::Math::Vector<double>         SH_to_coeffs(MR::Math::Vector<double> spherical_harmonics, bool include_isotropic);

        static Model                            factory(const MR::Math::Matrix<double>& encodings_matrix,
                                                        const MR::Math::Matrix<double>& response_SHs,
                                                        double adc,
                                                        double fa,
                                                        bool include_isotropic);

      //Public member functions
      public:

        Model() {}

        //TODO: Use b_value from encoding matrix to generate response function per direction.
        Model (const MR::Math::Matrix<double>& encodings_matrix, double adc, double fa, bool include_isotropic);

        //TODO: Use take a matrix of response_SH, allowing a different response function per direction.
        Model (const MR::Math::Matrix<double>& encodings_matrix, const MR::Math::Vector<double>& response_SH, bool include_isotropic);


        Model (const MR::Math::Matrix<double>& encodings_matrix, const MR::Math::Matrix<double>& response_SHs, bool include_isotropic);

        
        Model (const Model& m)
          : responses(m.responses), includes_iso(m.includes_iso)
          {}
      
      
        Model&                              operator= (const Model& m) {

          this->responses = m.responses;
          this->includes_iso = m.includes_iso;
        
          return *this;
        }
      

        const static double                 default_strand_volume_fraction = 0.2262 * .0005 * 1e5 * 2.0;

        size_t                              num_encodings() const
          { return responses.size(); }
        size_t                              size() const
          { return responses.size(); }
        
        Response&                           operator[](size_t index)
          { return responses[index]; }
        const Response&                     operator[](size_t index) const
          { return responses[index]; }

        const Encoding&                     encoding(size_t index) const
          { return responses[index]; }
        
        void                                precalculate_weightings(Fibre::Strand::Section& section) const;

        void                                precalculate_weightings_and_gradients(Fibre::Strand::Section& section) const;

        void                                precalculate_weightings_gradients_and_hessians(Fibre::Strand::Section& section) const;

        void                                scale_coeffs(double scalar);
        
        bool                                includes_isotropic() const
          { return includes_iso; }

        MR::Math::Vector<double>            weightings(MR::Math::Vector<double> weightings,
                                                                                MR::Math::Matrix<double> orientations);

      protected:

        void                                init(const MR::Math::Matrix<double>& encodings_matrix, const MR::Math::Matrix<double>& response_SHs, bool include_isotropic);

      friend std::ostream& operator<< (std::ostream& stream, const Diffusion::Model& model);

    };

  }
}

#undef LOOP

#endif
