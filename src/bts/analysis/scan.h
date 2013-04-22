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

#ifndef __bts_analysis_scan_h__
#define __bts_analysis_scan_h__

#include "bts/common.h"

#include "progressbar.h"
#include "bts/prob/prior.h"
#include "bts/prob/likelihood.h"

#include "image/header.h"
#include "dataset/loop.h"

namespace BTS {
    
    namespace Analysis {
        
        template<typename T> void scan(Prob::Likelihood& likelihood, Prob::Prior& prior,
                                       const T& origin, const T& axis, size_t num_steps,
                                       const std::string& output_location,
                                       std::map<std::string, std::string>& run_properties);
        
        template<typename T> void scan(Prob::Likelihood& likelihood, Prob::Prior& prior,
                                       const Fibre::Base::Set<T>& sequence,
                                       const std::string& output_location,
                                       std::map<std::string, std::string>& run_properties);
        
        template<typename T> void scan(Prob::Likelihood& likelihood, Prob::Prior& prior,
                                       const T& origin, const T& axis1, const T& axis2,
                                       size_t num_steps1, size_t num_steps2,
                                       const std::string& output_location,
                                       std::map<std::string, std::string>& properties);
        
        template<typename T> void scan(Prob::Likelihood& likelihood, Prob::Prior& prior,
                                       const T& origin, const T& axis1, const T& axis2,
                                       const T& axis3, size_t num_steps1, size_t num_steps2,
                                       size_t num_steps3, const std::string& output_location,
                                       std::map<std::string, std::string>& properties);
    
    }

}

namespace BTS {
    
    namespace Analysis {
        
        template<typename T> void scan(Prob::Likelihood& likelihood, Prob::Prior& prior,
                                       const T& origin, const T& axis, size_t num_steps,
                                       const std::string& output_location,
                                       std::map<std::string, std::string>& run_properties,
                                       bool save_gradient) {
            
            if (origin.size() != axis.size())
                throw Exception(
                        "Size of origin state (" + str(origin.size())
                        + ") does not match size of axis state (" + str(axis.size()) + ").");
            
            std::string gradient_location = File::strip_extension(output_location) + ".gradient."
                                            + File::extension(output_location);
            std::vector<std::string> header;
            header.push_back("log_px");
            header.push_back("likelihood_px");
            header.push_back("prior_px");
            
            std::vector<std::string> components_list = prior.list_components();
            
            header.insert(header.end(), components_list.begin(), components_list.end());
            
            typename T::Writer writer(output_location, origin, header, run_properties);
            typename T::Writer gradient_writer(gradient_location, origin);    //FIXME: This file shouldn't be created unless 'save_gradient' flag is set.
                    
            MR::ProgressBar progress_bar("Scanning over 1 dimension...", num_steps);
            
            double inc = 2.0 / (double) (num_steps - 1);
            
            for (double frac = -1.0; frac <= 1.0 + inc / 2.0; frac += inc) {    // The +inc/2.0 safeguards against rounding errors.
                    
                T gradient(origin), prior_gradient(origin), likelihood_gradient(origin);
                
                T state = origin + axis * frac;
                
                double prior_px;
                double likelihood_px;
                
                if (save_gradient) {
                    gradient.zero();
                    
                    prior_px = prior.log_prob(state, prior_gradient);
                    gradient += prior_gradient;
                    
                    likelihood_px = likelihood.log_prob(state, likelihood_gradient);
                    gradient += likelihood_gradient;
                    
                    gradient_writer.append(gradient);
                    
                } else {
                    prior_px = prior.log_prob(state);
                    likelihood_px = likelihood.log_prob(state);
                }
                
                double lprob = prior_px + likelihood_px;
                
                std::map<std::string, std::string> properties;
                
                properties["log_px"] = str(lprob);
                properties["likelihood_px"] = str(likelihood_px);
                properties["prior_px"] = str(prior_px);
                
                std::map<std::string, double> component_values = prior.get_component_values(state);
                
                for (std::map<std::string, double>::iterator comp_it = component_values.begin();
                        comp_it != component_values.end(); ++comp_it)
                    properties[comp_it->first] = str(comp_it->second);
                
                writer.append(state, properties);
                
                ++progress_bar;
                
            }
            
        }
        
        template<typename T> void scan(Prob::Likelihood& likelihood, Prob::Prior& prior,
                                       const Fibre::Base::Set<T>& sequence,
                                       const std::string& output_location,
                                       std::map<std::string, std::string>& run_properties) {
            
            std::string gradient_location = File::strip_extension(output_location) + ".gradient."
                                            + File::extension(output_location);
            std::vector<std::string> header;
            header.push_back("log_px");
            header.push_back("likelihood_px");
            header.push_back("prior_px");
            header.push_back("acs_px");
            
            typename T::Writer writer(output_location, run_properties, header);
            typename T::Writer gradient_writer(gradient_location);
            
            MR::ProgressBar progress_bar("Scanning over 1 dimension...", sequence.size());
            
            for (typename Fibre::Base::Set<T>::const_iterator it = sequence.begin();
                    it != sequence.end(); ++it) {
                
                T gradient(*it), prior_gradient(*it), likelihood_gradient(*it);
                
                T& state = *it;
                
                gradient.zero();
                
                double prior_px = prior.log_prob(state, prior_gradient);
                gradient += prior_gradient;
                
                double likelihood_px = likelihood.log_prob(state, likelihood_gradient);
                gradient += likelihood_gradient;
                
                double log_prob = prior_px + likelihood_px;
                
                std::map<std::string, std::string> properties_row;
                properties_row["log_px"] = str(log_prob);
                properties_row["likelihood_px"] = str(likelihood_px);
                properties_row["prior_px"] = str(prior_px);
                
                writer.append(state, properties_row);
                gradient_writer.append(gradient, std::map<std::string, std::string>());
                
                ++progress_bar;
                
            }
            
        }
        
        template<typename T> void scan(Prob::Likelihood& likelihood, Prob::Prior& prior,
                                       const T& origin, const T& axis1, const T& axis2,
                                       size_t num_steps1, size_t num_steps2,
                                       const std::string& output_location,
                                       std::map<std::string, std::string>& properties) {
            
            if (origin.size() != axis1.size())
                throw Exception(
                        "Size of origin state (" + str(origin.size())
                        + ") does not match size of axis1 state (" + str(axis1.size()) + ").");
            
            if (origin.size() != axis2.size())
                throw Exception(
                        "Size of origin state (" + str(origin.size())
                        + ") does not match size of axis2 state (" + str(axis2.size()) + ").");
            
            MR::Image::Header header;
            
            header.insert(properties.begin(), properties.end());
            
            header.set_ndim(2);
            
            header.set_dim(0, num_steps1);
            header.set_dim(1, num_steps2);
            
            header.set_description(0, File::strip_extension(properties["axis1_location"]));
            header.set_description(1, File::strip_extension(properties["axis2_location"]));
            
            header.set_vox(0, 2.0 / (double) num_steps1);
            header.set_vox(1, 2.0 / (double) num_steps2);
            
            File::clear_path(output_location);
            
            header.create(output_location);
            
            double inc1 = 2.0 / (double) (num_steps1 - 1);
            double inc2 = 2.0 / (double) (num_steps2 - 1);
            
            MR::Image::Voxel<double> pixel(header);
            
            MR::DataSet::Loop loop("Scanning over 2 dimensions...", 0, 2);
            loop.start(pixel);
            
            for (double frac1 = -1.0; frac1 <= 1.0 + inc1 / 2.0; frac1 += inc1) {    // The +inc/2.0 safeguards against rounding errors.
                    
                for (double frac2 = -1.0; frac2 <= 1.0 + inc2 / 2.0; frac2 += inc2) {
                    
                    T state = origin + axis1 * frac1 + axis2 * frac2;
                    
                    pixel.value() = prior.log_prob(state) + likelihood.log_prob(state);
                    
                    loop.next(pixel);
                    
                }
                
            }
            
        }
        
        template<typename T> void scan(Prob::Likelihood& likelihood, Prob::Prior& prior,
                                       const T& origin, const T& axis1, const T& axis2,
                                       const T& axis3, size_t num_steps1, size_t num_steps2,
                                       size_t num_steps3, const std::string& output_location,
                                       std::map<std::string, std::string>& properties) {
            
            if (origin.size() != axis1.size())
                throw Exception(
                        "Size of origin state (" + str(origin.size())
                        + ") does not match size of axis1 state (" + str(axis1.size()) + ").");
            
            if (origin.size() != axis2.size())
                throw Exception(
                        "Size of origin state (" + str(origin.size())
                        + ") does not match size of axis2 state (" + str(axis2.size()) + ").");
            
            if (origin.size() != axis3.size())
                throw Exception(
                        "Size of origin state (" + str(origin.size())
                        + ") does not match size of axis3 state (" + str(axis3.size()) + ").");
            
            MR::Image::Header header;
            
            header.insert(properties.begin(), properties.end());
            
            header.set_ndim(3);
            
            header.set_dim(0, num_steps1);
            header.set_dim(1, num_steps2);
            header.set_dim(2, num_steps3);
            
            header.set_description(0, File::strip_extension(properties["axis1_location"]));
            header.set_description(1, File::strip_extension(properties["axis2_location"]));
            header.set_description(2, File::strip_extension(properties["axis3_location"]));
            
            header.set_vox(0, 2.0 / (double) num_steps1);
            header.set_vox(1, 2.0 / (double) num_steps2);
            header.set_vox(2, 2.0 / (double) num_steps3);
            
            File::clear_path(output_location);
            
            header.create(output_location);
            
            double inc1 = 2.0 / (double) (num_steps1 - 1);
            double inc2 = 2.0 / (double) (num_steps2 - 1);
            double inc3 = 2.0 / (double) (num_steps3 - 1);
            
            MR::Image::Voxel<double> pixel(header);
            
            MR::DataSet::Loop loop("Scanning over 3 dimensions...", 0, 3);
            loop.start(pixel);
            
            for (double frac1 = -1.0; frac1 <= 1.0 + inc1 / 2.0; frac1 += inc1) {    // The +inc/2.0 safeguards against rounding errors.
                    
                for (double frac2 = -1.0; frac2 <= 1.0 + inc2 / 2.0; frac2 += inc2) {
                    
                    for (double frac3 = -1.0; frac3 <= 1.0 + inc3 / 2.0; frac3 += inc3) {
                        
                        T state = origin + axis1 * frac1 + axis2 * frac2 + axis3 * frac3;
                        
                        pixel.value() = prior.log_prob(state) + likelihood.log_prob(state);
                        
                        loop.next(pixel);
                        
                    }
                    
                }
                
            }
            
        }
    
    }

}

#endif
