/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Written by Thomas G. Close, 04/03/2009.

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

using namespace BTS;

const double THRESHOLD_DEFAULT = -0.5;
const size_t NUM_LENGTH_SECTIONS_DEFAULT = 20;
const double TRACT_WIDTH_DEFAULT = 1.0;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Split Fourier tracts into two to allow them to follow branching structure",
    "",
    NULL
};

ARGUMENTS= {
    Argument ("input", "The strands to be selected from.").type_file (),
    Argument ("output", "The selected strands.").type_file (),
    Argument()
};

OPTIONS= {

    Option ("orientation", "The orientation of the plane to split the Fourier tract into.")
    + Argument ("offsets", "").type_text ("auto"),

    Option ("num_length_sections", "The number of strands to split the tracts into")
    + Argument ("num_length_sections", "").type_integer (1, NUM_LENGTH_SECTIONS_DEFAULT, LARGE_INT),

    Option ("tract_width", "The number of strands to split the tracts into")
    + Argument ("tract_width", "").type_integer (SMALL_FLOAT, TRACT_WIDTH_DEFAULT, LARGE_FLOAT),

    Option ("include","The indices of the strands to include")
    + Argument ("include","").type_text(),

    Option()};

EXECUTE {
    
        std::string input_location = argument[0];
        std::string output_location = argument[1];
        
        Triple<double> orient(1.0, 0.0, 0.0);
        size_t num_length_sections = NUM_LENGTH_SECTIONS_DEFAULT;
        double constant_tract_width = -1.0;
        std::vector<size_t> include;
        
        Options opt = get_options("orientation");
        if (opt.size())
            orient = parse_triple<double>(std::string(opt[0][0])).normalise();
        
        opt = get_options("num_length_sections");
        if (opt.size())
            num_length_sections = opt[0][0];

        opt = get_options("tract_width");
        if (opt.size())
            constant_tract_width = opt[0][0];
        
        opt = get_options("include");
        if (opt.size())
            include = parse_sequence<size_t>(opt[0][0]);
        
        Fibre::Tractlet::Set tracts(input_location);
        Fibre::Tractlet::Set output_tracts(0, tracts.degree());
        output_tracts.add_elem_prop(Fibre::Base::Object::ALPHA_PROP);
        
        MR::ProgressBar progress_bar("Splitting tracts in two.", tracts.size());

        for (size_t tract_i = 0; tract_i < tracts.size(); ++tract_i) {
            
            if (!include.size() || std::find(include.begin(), include.end(), tract_i)
                    != include.end()) {
                
                const Fibre::Tractlet& tract = tracts[tract_i];

                // Get the extent of the auxiliary axis at 'num_points' along its path.
                Fibre::Track backbone = tract[0].to_track(num_length_sections);
                Fibre::Track aux_tck1 = tract[1].to_track(num_length_sections);
                Fibre::Track aux_tck2 = tract[2].to_track(num_length_sections);

                Fibre::Track tck1(backbone);
                Fibre::Track tck2(backbone);

                for (size_t i = 0; i < num_length_sections; ++i) {

                    double norm1 = aux_tck1[i].norm();
                    double norm2 = aux_tck2[i].norm();

                    Triple<double> aux1 = aux_tck1[i] / norm1;
                    Triple<double> aux2 = aux_tck2[i] / norm2;

                    double proj1 = aux1.dot(orient);
                    double proj2 = aux2.dot(orient);

                    aux1 *= Math::sign(proj1);
                    aux2 *= Math::sign(proj2);

                    double theta1 = MR::Math::acos(fabs(proj1));
                    double theta2 = MR::Math::acos(fabs(proj2));

                    double angle_period = theta1 + theta2;

                    double frac1;
                    double frac2;

                    if (angle_period == 0.0) {
                        frac1 = 1.0;
                        frac2 = 1.0;
                    } else {
                        frac1 = cos(theta1 * M_PI / (4.0 * angle_period));
                        frac2 = cos(theta2 * M_PI / (4.0 * angle_period));
                    }

                    double combined_proj = MR::Math::sqrt(MR::Math::pow2(frac1 * proj1 * norm1) +
                                                          MR::Math::pow2(frac2 * proj2 * norm2));

                    tck1[i] += orient * combined_proj;
                    tck2[i] -= orient * combined_proj;

                }

                double tract_width;
                if (constant_tract_width > 0.0)
                    tract_width = constant_tract_width;
                else
                    tract_width = MR::Math::sqrt(tract.acs() / 2.0);
                output_tracts.push_back(tck1.to_strand(tract.degree()).to_tractlet(tract_width).set_acs(tract.acs()));
                output_tracts.push_back(tck2.to_strand(tract.degree()).to_tractlet(tract_width).set_acs(tract.acs()));
            }
            ++progress_bar;
        }
        
        output_tracts.save(output_location);
        
    }
    
