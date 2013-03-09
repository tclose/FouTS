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

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
}

#include "bts/prob/test/landscape.h"
#include "bts/prob/test/landscape/peak.h"
#include "bts/prob/test/landscape/peak/writer.h"
#include "bts/prob/test/landscape/peak/reader.h"

#include "bts/prob/inline_functions.h"

#include "bts/image/buffer.cpp.h"
#include "bts/math/common.h"

namespace BTS {
    
    namespace Prob {
        
        namespace Test {
            
            const char* Landscape::LOCATION_DEFAULT = "/home/tclose/data/prob/test/landscape/a.lnd";
            const double Landscape::WIDTH_MU_DEFAULT = 10;
            const double Landscape::WIDTH_SIGMA_DEFAULT = 5;
            const double Landscape::WIDTH_MIN_DEFAULT = 15;
            const double Landscape::HEIGHT_VAR_DEFAULT = 6;
            const double Landscape::ROI_RADIUS_DEFAULT = 10.0;
            const double Landscape::FRACTION_PYRAMID_DEFAULT = 0.0;
            const double Landscape::BARRIER_RATE_DEFAULT = 0.3;
            
            Landscape::Landscape(double roi_radius, double barrier_rate)
                    : roi_radius(roi_radius), barrier_rate(barrier_rate) {
            }
            
            Landscape::Landscape(const std::string& location, double roi_radius,
                                 double barrier_rate)
                    : roi_radius(roi_radius), barrier_rate(barrier_rate) {
                load(location);
            }
            
            Landscape Landscape::randomly_generate(size_t ndims, size_t num_peaks,
                                                   double peak_width_mu, double peak_width_sigma,
                                                   double peak_width_min, double peak_height_var,
                                                   double roi_radius, double fraction_pyramid,
                                                   double barrier_rate, size_t seed) {
                
                Landscape landscape;
                
                landscape.ndims = ndims;
                landscape.roi_radius = roi_radius;
                landscape.barrier_rate = barrier_rate;
                
                gsl_rng* rand_gen = gsl_rng_alloc(gsl_rng_taus);
                gsl_rng_set(rand_gen, seed);
                
                for (size_t peak_i = 0; peak_i < num_peaks; peak_i++) {
                    
                    // Rejection sampling to get peaks randomly distributed withinness the roi_radius.
                    MCMC::State centre(ndims);
                    
                    for (size_t dim_i = 0; dim_i < ndims; ++dim_i)
                        centre[dim_i] = roi_radius;
                    
                    while (centre.norm() >= roi_radius)
                        for (size_t dim_i = 0; dim_i < ndims; ++dim_i)
                            centre[dim_i] = gsl_ran_flat(rand_gen, -roi_radius, roi_radius);
                    
                    MCMC::State::Tensor widths(ndims);
                    
                    MR::Math::Matrix<double> eig_matrix(ndims, ndims), working(ndims, ndims);
                    
                    eig_matrix = 0;
                    
                    for (size_t dim_i = 0; dim_i < 2; dim_i++) {
                        double width =
                                MR::Math::sqrt(
                                        MR::Math::pow2(
                                                gsl_ran_gaussian(rand_gen, peak_width_sigma) + peak_width_mu
                                                        * (1.0 - peak_width_min))
                                        + MR::Math::pow2(
                                                gsl_ran_gaussian(rand_gen, peak_width_sigma) + peak_width_mu
                                                        * (1.0 - peak_width_min)))
                                / M_SQRT2;
                        
                        width += peak_width_min * peak_width_mu;
                        
                        eig_matrix(dim_i, dim_i) = 3.0 / MR::Math::pow2(width);
                    }
                    
                    MR::Math::Matrix<double> rotate_matrix = Math::random_rotation(rand_gen, ndims);
                    
                    MR::Math::mult(working, eig_matrix, rotate_matrix);
                    MR::Math::mult(widths, 1.0, CblasTrans, rotate_matrix, CblasNoTrans, working);
                    
                    // Rayleigh distribution for peak height.
                    double height = gsl_ran_gaussian(rand_gen, peak_height_var);
                    
                    size_t type = 0;
                    
                    // Select the peak to be a pyramid instead of a gaussian with probability = fraction_pyramid.
                    if (gsl_ran_flat(rand_gen, 0.0, 1.0) < fraction_pyramid)
                        type = 1;
                    
                    // Create peak object.
                    Peak peak(ndims, centre, widths, height, type);
                    
                    landscape.peaks.push_back(peak);
                    
                }
                
                gsl_rng_free(rand_gen);
                
                return landscape;
                
            }
            
            Landscape::~Landscape() {
            }
            
            double Landscape::log_prob(const MCMC::State& state) {
                
                MCMC::State dummy;
                
                return log_prob(state, dummy);
                
            }
            
            double Landscape::log_prob(const MCMC::State& state, MCMC::State& gradient) {
                
                double log_prob = 0.0;
                gradient.zero();
                
                MCMC::State peak_gradient;
                
                for (std::vector<Peak>::iterator peak_it = peaks.begin(); peak_it != peaks.end();
                        ++peak_it) {
                    
                    log_prob += peak_it->log_prob(state, peak_gradient);
                    
                    gradient += peak_gradient;
                    
                }
                
                if (state.norm() > roi_radius) {
                    
                    log_prob -= barrier_rate * MR::Math::pow2((state.norm() - roi_radius));
                    gradient -= 2.0 * barrier_rate * state * (state.norm() - roi_radius)
                            / state.norm();
                    
                }
                
                return log_prob;
                
            }
            
            double Landscape::log_prob(const MCMC::State& state, MCMC::State& gradient,
                                       MCMC::State::Tensor& hessian) {
                
                double log_prob = 0.0;
                gradient.resize(state.size());
                hessian.resize(state.size(), state.size());
                
                gradient.zero();
                hessian.zero();
                
                MCMC::State peak_gradient;
                MCMC::State::Tensor peak_hessian;
                
                for (std::vector<Peak>::iterator peak_it = peaks.begin(); peak_it != peaks.end();
                        ++peak_it) {
                    
                    log_prob += peak_it->log_prob(state, peak_gradient, peak_hessian);
                    
                    gradient += peak_gradient;
                    
                    hessian += peak_hessian;
                    
                }
                
                if (state.norm() > roi_radius) {
                    
                    log_prob -= barrier_rate * MR::Math::pow2((state.norm() - roi_radius));
                    gradient -= 2.0 * barrier_rate * state * (state.norm() - roi_radius)
                            / state.norm();
                    
                    MR::Math::Matrix<double> hess_part = BTS::Math::outer(state, state);
                    hess_part *= 2.0
                            * barrier_rate
                            * (-1.0 * (state.norm() - roi_radius) / MR::Math::pow3(state.norm()) + 1
                                    / state.norm2());
                    hessian -= hess_part;
                    
                    for (size_t elem_i = 0; elem_i < state.size(); elem_i++)
                        hessian(elem_i, elem_i) -= 2 * barrier_rate * (state.norm() - roi_radius)
                                / state.norm();
                    
                }
                
                return log_prob;
                
            }
            
            double Landscape::log_prob(const MCMC::State& state, MCMC::State& gradient,
                                       MCMC::State::Tensor& hessian,
                                       std::vector<MCMC::State::Tensor>& rank3_hessian) {
                
                double log_prob = 0.0;
                gradient.resize(state.size());
                hessian.resize(state.size(), state.size());
                
                gradient.zero();
                hessian.zero();
                
                MCMC::State peak_gradient;
                MCMC::State::Tensor peak_hessian;
                std::vector<MCMC::State::Tensor> peak_rank3_hessian;
                
                for (size_t dim_i = 0; dim_i < state.size(); dim_i++)
                    rank3_hessian.push_back(hessian);
                
                for (std::vector<Peak>::iterator peak_it = peaks.begin(); peak_it != peaks.end();
                        ++peak_it) {
                    
                    log_prob += peak_it->log_prob(state, peak_gradient, peak_hessian,
                            peak_rank3_hessian);
                    
                    gradient += peak_gradient;
                    
                    hessian += peak_hessian;
                    
                    for (size_t dim_i = 0; dim_i < state.size(); dim_i++)
                        rank3_hessian[dim_i] += peak_rank3_hessian[dim_i];
                    
                }
                
                if (state.norm() > roi_radius) {
                    
                    log_prob -= barrier_rate * MR::Math::pow2((state.norm() - roi_radius));
                    gradient -= 2.0 * barrier_rate * state * (state.norm() - roi_radius)
                            / state.norm();
                    
                    MR::Math::Matrix<double> hess_part = BTS::Math::outer(state, state);
                    hess_part *= -2.0
                            * barrier_rate
                            * (-1.0 * (state.norm() - roi_radius) / MR::Math::pow3(state.norm()) + 1
                                    / state.norm2());
                    hessian += hess_part;
                    
                    for (size_t elem_i = 0; elem_i < state.size(); elem_i++)
                        hessian(elem_i, elem_i) -= 2 * barrier_rate * (state.norm() - roi_radius)
                                / state.norm();
                    
                    for (size_t dim_i = 0; dim_i < state.size(); dim_i++) {
                        
                        MCMC::State::Tensor& d_hessian = rank3_hessian[dim_i];
                        
                        MR::Math::Matrix<double> d_scalar_hess = Math::outer(state, state);
                        d_scalar_hess *=
                                -2.0 * barrier_rate
                                * state[dim_i]
                                * (3.0 * (state.norm() - roi_radius) / MR::Math::pow5(state.norm()) - 3.0
                                        / MR::Math::pow4(state.norm()));
                        
                        MR::Math::Matrix<double> d_outer_hess(state.size(), state.size());
                        d_outer_hess = 0.0;
                        d_outer_hess.column(dim_i) = state;
                        d_outer_hess.row(dim_i) += state;
                        d_outer_hess *=
                                -2.0 * barrier_rate
                                * (-1.0 * (state.norm() - roi_radius) / MR::Math::pow3(state.norm()) + 1
                                        / state.norm2());
                        
                        d_hessian += d_scalar_hess;
                        d_hessian += d_outer_hess;
                        
                        for (size_t dim_i2 = 0; dim_i2 < state.size(); dim_i2++)
                            d_hessian(dim_i2, dim_i2) -= 2.0
                                    * barrier_rate
                                    * state[dim_i]
                                    * (1 / state.norm2() - (state.norm() - roi_radius)
                                            / MR::Math::pow3(state.norm()));
                        
                    }
                    
                }
                
                return log_prob;
                
            }
            
            void Landscape::save(const std::string& location) {
                
                Peak::Writer writer(location);
                
                for (std::vector<Peak>::iterator peak_it = peaks.begin(); peak_it != peaks.end();
                        ++peak_it)
                    writer.append(*peak_it);
                
            }
            
            void Landscape::load(const std::string& location) {
                
                this->peaks.clear();
                
                Peak::Reader reader(location);
                
                Peak peak;
                
                ndims = 0;
                
                while (reader.next(peak)) {
                    
                    if (!ndims)
                        ndims = peak.num_dims();
                    else if (ndims != peak.num_dims())
                        throw Exception(
                                "Peak dimension (" + str(peak.num_dims())
                                + ") does not match overall dimension (" + str(ndims) + ").");
                    
                    this->peaks.push_back(peak);
                    
                }
                
                if (!ndims)
                    throw Exception(
                            "Loaded test landscape ('" + location
                            + "') does not have any dimension.");
                
            }
            
            std::ostream& operator<<(std::ostream& stream, const Landscape& test) {
                
                for (std::vector<Landscape::Peak>::const_iterator peak_it = test.peaks.begin();
                        peak_it != test.peaks.end(); ++peak_it)
                    stream << *peak_it << std::endl;
                
                return stream;
                
            }
        
        }
    
    }

}
