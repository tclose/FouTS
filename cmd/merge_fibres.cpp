/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Written by Thomas G. Close, 04/03/2009.

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

#include "bts/cmd.h"
#include "point.h"
#include "progressbar.h"

#include "bts/common.h"
#include "bts/file.h"

#include "bts/fibre/strand/set.h"
#include "bts/fibre/track/set.h"
#include "bts/fibre/tractlet/set.h"

#include "bts/diffusion/model.h"
#include "bts/image/expected/trilinear/buffer.h"
#include "bts/image/expected/gaussian/buffer.h"
#include "bts/image/expected/quartic/buffer.h"

#include "bts/prob/uniform.h"
#include "bts/prob/prior.h"

#include "bts/math/common.h"

#include "bts/prob/likelihood.h"
#include "bts/prob/likelihood.h"
#include "bts/prob/likelihood/one_sided_gaussian.h"
#include "bts/prob/likelihood/gaussian.h"

#include "bts/image/observed/buffer.h"

#include "bts/image/expected/buffer.h"

#include "bts/math/blossom.h"
#include "bts/math/blossom/edge.h"

#include "bts/inline_functions.h"

using namespace FTS;

const double THRESHOLD_DEFAULT = -0.5;
const size_t NUM_LENGTH_SECTIONS_DEFAULT = 5;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Merge strands that can be better represented by a single strand",
    "",
    NULL
};

ARGUMENTS= {
    Argument ("observed image", "The observed image the strands will be fit against.").type_image_in (),
    Argument ("input", "The strands to be selected from.").type_file (),
    Argument ("output", "The selected strands.").type_file (),
    Argument()
};

OPTIONS= {
//
//    Option ("threshold","Threshold above which to strands will be merged together.")
//    + Argument ("","").type_float (-1, Fibre::Tractlet::STRANDS_PER_AREA_DEFAULT, LARGE_FLOAT),

    Option ("strands_per_acs","Basewidth of the imported tractlets.")
    + Argument ("","").type_float (-1, THRESHOLD_DEFAULT, LARGE_FLOAT),

    Option ("num_tractlets", "Number of tractlets to convert the strands into if outputting tractlets (default to number of input tractlets")
    + Argument ("num_tractlets", "").type_integer (0, 0, LARGE_INT),

    Option ("num_length_sections", "Number of samples to take along the strand path")
    + Argument ("num_length_sections", "").type_integer (1, NUM_LENGTH_SECTIONS_DEFAULT, LARGE_INT),

    Option ("new_degree", "The degree the merged strands will have")
    + Argument ("new_degree", "").type_integer (1, Fibre::Strand::DEFAULT_DEGREE, LARGE_INT),

    DIFFUSION_PARAMETERS,

    IMAGE_PARAMETERS,

    EXPECTED_IMAGE_PARAMETERS,

    LIKELIHOOD_PARAMETERS,

    PRIOR_PARAMETERS,

    Option()};

EXECUTE {
    
        std::string obs_image_location = argument[0];
        std::string input_location = argument[1];
        std::string output_location = argument[2];
        
        MR::Image::Header header(argument[0]);
        
        if (header.ndim() != 4)
            throw Exception("dwi image should contain 4 dimensions");
        
//        float threshold = THRESHOLD_DEFAULT;
        float strands_per_acs = Fibre::Tractlet::STRANDS_PER_AREA_DEFAULT;
        size_t num_tractlets = 0;
        size_t num_length_sections = NUM_LENGTH_SECTIONS_DEFAULT;
        size_t new_degree = Fibre::Strand::DEFAULT_DEGREE;
        
//        Options opt = get_options("threshold");
//        if (opt.size())
//            threshold = opt[0][0];
        
        Options opt = get_options("strands_per_acs");
        if (opt.size())
            strands_per_acs = opt[0][0];
        
        opt = get_options("num_tractlets");
        if (opt.size())
            num_tractlets = opt[0][0];
        
        opt = get_options("num_length_sections");
        if (opt.size())
            num_length_sections = opt[0][0];
        
        opt = get_options("new_degree");
        if (opt.size())
            new_degree = opt[0][0];
        
        // Loads parameters to construct Diffusion::Model ('diff_' prefix)
        SET_DIFFUSION_PARAMETERS;
        
        // Loads parameters to construct Image::Expected::*::Buffer that are inherited from Image::Observed::Buffer ('img_' prefix)
        SET_IMAGE_PARAMETERS;
        
        // Loads parameters to construct Image::Expected::*::Buffer ('img_' prefix)
        SET_EXPECTED_IMAGE_PARAMETERS
        ;
        
        // Loads parameters to construct Prob::Likelihood ('like_' prefix)
        SET_LIKELIHOOD_PARAMETERS
        ;
        
        // Loads parameters to construct Prob::Prior ('prior_' prefix)
        SET_PRIOR_PARAMETERS
        ;
        
        Image::Observed::Buffer obs_image(obs_image_location,
                Diffusion::Encoding::Set(diff_encodings));
        
        if (obs_image.properties().count("diff_response_SH")) {
            
            if (Math::matlab_str(diff_response_SH) != obs_image.properties()["diff_response_SH"])
                std::cout << std::endl << "Warning! Diffusion response function harmonics ("
                          << Math::matlab_str(diff_response_SH)
                          << ") do not match reference image ("
                          << obs_image.properties()["diff_response_SH"] + ")!" << std::endl;
            
        }
        
        //If gradient scheme is included in reference image header, use that instead of default (NB: Will override any gradients passed to '-diff_encodings' option).
        if (header.get_DW_scheme().rows()) {
            diff_encodings = header.get_DW_scheme();
            diff_encodings_location = "From observed image";
        }
        
        if ((header.count("diff_response_SH")) && (Math::matlab_str(diff_response_SH)
                != header["diff_response_SH"]))
            std::cout << std::endl << "Warning! Diffusion response function harmonics ("
                      << Math::matlab_str(diff_response_SH) << ") do not match reference image ("
                      << header["diff_response_SH"] + ")!" << std::endl;
        
        //----------------------------//
        //  Initialize Expected Image //
        //----------------------------//
        
        Diffusion::Model diffusion_model = Diffusion::Model::factory(diff_encodings,
                diff_response_SH, diff_adc, diff_fa, diff_isotropic, diff_warn_b_mismatch);
        
        Image::Expected::Buffer& exp_image = *Image::Expected::Buffer::factory(exp_type, obs_image,
                diffusion_model, exp_num_length_sections, exp_num_width_sections, exp_interp_extent,
                exp_enforce_bounds, exp_half_width);
        
        Image::Expected::Buffer& diff_image = *exp_image.clone();
        
        //-----------------------//
        // Initialize Likelihood //
        //-----------------------//
        
        Prob::Likelihood* likelihood = Prob::Likelihood::factory(like_type, obs_image, &exp_image,
                like_snr, like_b0_include, like_outside_scale, like_ref_b0, like_ref_signal, like_noise_map);
        
        //------------------//
        // Initialize Prior //
        //------------------//
        
        Prob::Prior prior(prior_scale, prior_freq_scale, prior_freq_aux_scale, prior_hook_scale,
                prior_hook_num_points, prior_hook_num_width_sections, prior_density_high_scale,
                prior_density_low_scale, prior_density_num_points, prior_acs_scale, prior_acs_mean,
                prior_length_scale, prior_length_mean, prior_in_image_scale, prior_in_image_power,
                Prob::PriorComponent::InImage::get_offset(obs_image, prior_in_image_border),
                Prob::PriorComponent::InImage::get_extent(obs_image, prior_in_image_border),
                prior_in_image_num_length_sections, prior_in_image_num_width_sections);
        
        Fibre::Strand::Set strands(input_location, strands_per_acs);
        
        if (exp_base_intensity)
            strands.set_base_intensity(exp_base_intensity);
        
//        strands.save("/home/tclose/data/input_strands.str");
        
        if ((num_tractlets == 0) && strands.has_extend_prop(
                    Fibre::Strand::Set::ORIGINAL_NUM_TRACTLETS_PROP))
            num_tractlets = to<size_t>(
                    strands.get_extend_prop(Fibre::Strand::Set::ORIGINAL_NUM_TRACTLETS_PROP));
        
        exp_image.expected_image(strands);
        
//        obs_image.save("/home/tclose/data/orig_observed.mif");
//        exp_image.save("/home/tclose/data/orig_expected.mif");
        
        double orig_like = likelihood->log_prob(strands);
        
        std::set<Math::Blossom::Edge> edges;
        
        MR::ProgressBar progress_bar("Calculating merge weights...", strands.size());
        
        for (size_t row_i = 0; row_i < strands.size(); ++row_i) {
            for (size_t col_i = row_i + 1; col_i < strands.size(); ++col_i) {
                
                bool same_bundle = false;
                
                if (strands.has_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP))
                    if (to<size_t>(
                                strands.get_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP,
                                        row_i))
                        == to<size_t>(
                                strands.get_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP,
                                        col_i)))
                        same_bundle = true;
                
                if (!same_bundle) {
                    
                    Fibre::Strand row_strand = strands[row_i];
                    Fibre::Strand col_strand = strands[col_i];
                    
                    Fibre::Strand merged = row_strand.merge(col_strand, num_length_sections,
                            new_degree);
                    
                    diff_image.zero();
                    diff_image.part_image(row_strand);
                    
                    diff_image.part_image(col_strand);
                    
                    diff_image.negate();
                    
                    diff_image.part_image(merged);
                    
                    exp_image += diff_image;
                    
                    //The distance measure. The original prior is averaged by number of strands, and the merged prior and likelihood
                    //are subtracted from it so that unlikely matches are assigned a large distance.
                    
                    double orig_prior = prior.log_prob(row_strand) + prior.log_prob(col_strand);
                    double merged_prior = prior.log_prob(merged);
                    double merged_like = likelihood->log_prob(exp_image);
                    
//          exp_image.save("/home/tclose/data/merged_expected.mif");
                    
                    double diff_lprob = (merged_prior - orig_prior) + (merged_like - orig_like);
                    
                    //If difference in log prob is greater than zero. i.e. if merging improves probability add to list of edges to
                    //find the best fit.
                    
//          if (!(row_i == 3 && col_i == 8))
                    if (diff_lprob > 0)
                        edges.insert(Math::Blossom::Edge(row_i, col_i, -diff_lprob));
                    
                    exp_image -= diff_image;
                    
                }
                
            }
            
            ++progress_bar;
        }
        
        std::vector<Math::Blossom::Edge> merge_list;
        
        for (std::set<Math::Blossom::Edge>::iterator edge_it = edges.begin();
                edge_it != edges.end();) {
            
            merge_list.push_back(*edge_it);
            
            std::set<Math::Blossom::Edge>::iterator edge_it2 = edge_it;
            
            ++edge_it2;
            
            while (edge_it2 != edges.end())
                if (edge_it2->has_node(edge_it->first) || edge_it2->has_node(edge_it->second))
                    edges.erase(edge_it2++);
                else
                    ++edge_it2;
            
            edges.erase(edge_it++);
        }
        
        std::set<size_t> not_merge_list;
        
        for (size_t strand_i = 0; strand_i < strands.size(); ++strand_i)
            not_merge_list.insert(strand_i);
        
        double expected_improvement = 0;
        
        std::cout << std::endl;
        
        for (std::vector<Math::Blossom::Edge>::iterator merge_it = merge_list.begin();
                merge_it != merge_list.end(); ++merge_it) {
            not_merge_list.erase(merge_it->first);
            not_merge_list.erase(merge_it->second);
            expected_improvement -= merge_it->weight;
            std::cout << merge_it->first << " " << merge_it->second << ": " << -merge_it->weight
                      << std::endl;
        }
        
        Fibre::Strand::Set merged;
        
        for (std::vector<Math::Blossom::Edge>::iterator merge_it = merge_list.begin();
                merge_it != merge_list.end(); ++merge_it) {
            merged.push_back(
                    strands[merge_it->first].merge(strands[merge_it->second], num_length_sections,
                            new_degree));
        }
        
        for (std::set<size_t>::iterator not_merge_it = not_merge_list.begin();
                not_merge_it != not_merge_list.end(); ++not_merge_it) {
            Fibre::Strand strand = strands[*not_merge_it];
            strand.resize(new_degree, 0.0);
            merged.push_back(strand);
            
        }
        
        double improvement = likelihood->log_prob(merged) + prior.log_prob(merged)
                - likelihood->log_prob(strands) - prior.log_prob(strands);
        
//        merged.save("/home/tclose/data/merged.str");
        
//    if (improvement > 0) {
        merged.save(output_location, 0, num_tractlets);
        std::cout << "\n\nMerged " << merge_list.size()
                  << " strand pairs for a combined improvement of the posterior of " << improvement
                  << " (" << expected_improvement << " expected) over the original set."
                  << std::endl;
//    } else {
//      strands.save(output_location,0,num_tractlets);
//      std::cout << "\n\nNo improvement could be found from the original set." << std::endl;
//    }
        
    }
    
//}
    
