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

#ifndef __bts_prob_testpeaks__
#define __bts_prob_testpeaks__

#define TEST_LANDSCAPE_PARAMETERS \
\
  Option ("lnd_location", "") \
    + Argument ("lnd_location", "").type_text (Prob::Test::Landscape::LOCATION_DEFAULT), \
\
  Option ("lnd_scale", "") \
    + Argument ("lnd_scale", "").type_float (SMALL_FLOAT, LARGE_FLOAT, 1.0), \
\
  Option ("lnd_roi_radius", "") \
    + Argument ("lnd_roi_radius", "").type_float (SMALL_FLOAT, LARGE_FLOAT,Prob::Test::Landscape::ROI_RADIUS_DEFAULT), \
\
  Option ("lnd_barrier_rate", "") \
    + Argument ("lnd_barrier_rate", "").type_float (SMALL_FLOAT, LARGE_FLOAT,Prob::Test::Landscape::BARRIER_RATE_DEFAULT) \

//Loads the 'prior' parameters into variables
#define SET_TEST_LANDSCAPE_PARAMETERS \
  std::string lnd_location    = Prob::Test::Landscape::LOCATION_DEFAULT; \
  double lnd_scale            = 1.0; \
  double lnd_roi_radius       = Prob::Test::Landscape::ROI_RADIUS_DEFAULT; \
  double lnd_barrier_rate     = Prob::Test::Landscape::BARRIER_RATE_DEFAULT; \
\
  opt = get_options("lnd_location"); \
  if (opt.size()) \
    lnd_location = opt[0][0].c_str(); \
\
  opt = get_options("lnd_scale"); \
  if (opt.size()) \
    lnd_scale = opt[0][0]; \
\
  opt = get_options("lnd_roi_radius"); \
  if (opt.size()) \
    lnd_roi_radius = opt[0][0]; \
\
  opt = get_options("lnd_barrier_rate"); \
  if (opt.size()) \
    lnd_barrier_rate = opt[0][0]; \

//Adds the 'prior' parameters to the properties to be saved with the data.
#define ADD_TEST_LANDSCAPE_PROPERTIES(properties) \
  properties["lnd_location"]          = lnd_location; \
  properties["lnd_scale"]             = str(lnd_scale); \
  properties["lnd_roi_radius"]        = str(lnd_roi_radius); \
  properties["lnd_barrier_rate"]        = str(lnd_barrier_rate);

#include "bts/prob/likelihood.h"
#include "bts/prob/prior.h"

//TODO: remove
#include "bts/image/expected/buffer.h"

namespace FTS {
    
    namespace Prob {
        
        namespace Test {
            
            class Landscape {
                    
                public:
                    
                    class Peak;

                public:
                    
                    const static char* LOCATION_DEFAULT;

                    const static size_t NUM_DIMS_DEFAULT = 2;
                    const static size_t NUM_PEAKS_DEFAULT = 40;
                    const static double WIDTH_MU_DEFAULT;
                    const static double WIDTH_SIGMA_DEFAULT;
                    const static double WIDTH_MIN_DEFAULT;
                    const static double HEIGHT_VAR_DEFAULT;
                    const static double ROI_RADIUS_DEFAULT;
                    const static double FRACTION_PYRAMID_DEFAULT;
                    const static double BARRIER_RATE_DEFAULT;

                    //TODO: remove
                    //          FTS::Image::Expected::Buffer* exp_image;
                    
                public:
                    
                    static Landscape randomly_generate(size_t num_dims, size_t num_peaks,
                                                       double peak_width_mu,
                                                       double peak_width_sigma,
                                                       double peak_width_min,
                                                       double peak_height_var, double roi_radius,
                                                       double fraction_pyramid,
                                                       double barrier_rate = BARRIER_RATE_DEFAULT,
                                                       size_t seed = time(NULL));

                protected:
                    
                    size_t ndims;

                    double roi_radius;

                    double barrier_rate;

                    std::vector<Peak> peaks;

                public:
                    
                    Landscape(double roi_radius = ROI_RADIUS_DEFAULT, double barrier_rate =
                            BARRIER_RATE_DEFAULT);

                    Landscape(const std::string& location, double roi_radius = ROI_RADIUS_DEFAULT,
                              double barrier_rate = BARRIER_RATE_DEFAULT);

                    ~Landscape();

                    Landscape(double frequency, size_t num_points, double spacing, size_t seed =
                            time(NULL));

                    double log_prob(const MCMC::State& test_state);

                    double log_prob(const MCMC::State& test_state, MCMC::State& gradient);

                    double log_prob(const MCMC::State& test_state, MCMC::State& gradient,
                                    MCMC::State::Tensor& hessian);

                    double log_prob(const MCMC::State& test_state, MCMC::State& gradient,
                                    MCMC::State::Tensor& hessian,
                                    std::vector<MCMC::State::Tensor>& rank3_hessian);

                    void save(const std::string& location);

                    void load(const std::string& location);

                    size_t num_dims() const {
                        return ndims;
                    }
                    
                    Peak& operator[](size_t index) {
                        return peaks[index];
                    }
                    
                    const Peak& operator[](size_t index) const {
                        return peaks[index];
                    }
                    
                    void set_assumed_snr(double assumed_snr, const std::string& ref_b0,
                                         double ref_signal) {
                    }
                    
                    void set_enforce_bounds(bool flag) {
                    }
                    
                    std::vector<std::string> list_components() {
                        return std::vector<std::string>();
                    }
                    
                    template<typename T> std::map<std::string, double> get_component_values(T& t) {
                        return std::map<std::string, double>();
                    }
                    
                    friend std::ostream& operator<<(std::ostream& stream, const Landscape& test);
                    
            };
            
            std::ostream& operator<<(std::ostream& stream, const Landscape& test);
        
        }
    
    }

}

#include "bts/prob/test/landscape/peak.h"

#endif
