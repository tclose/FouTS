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
const size_t NUM_STRANDS_DEFAULT = 4;
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

    Option ("num_strands", "The number of strands to split the tracts into")
    + Argument ("num_strands", "").type_text ("auto"),

    Option ("include","The indices of the strands to include")
    + Argument ("include","").type_text(),

    Option()};

EXECUTE {
    
        std::string input_location = argument[0];
        std::string output_location = argument[1];
        
        Triple<double> orient(1.0, 0.0, 0.0);
        size_t num_strands = NUM_STRANDS_DEFAULT;
        std::vector<size_t> include;
        
        Options opt = get_options("orientation");
        if (opt.size())
            orient = parse_triple<double>(std::string(opt[0][0])).normalise();
        
        opt = get_options("num_strands");
        if (opt.size())
            num_strands = opt[0][0];
        
        opt = get_options("include");
        if (opt.size())
            include = parse_sequence<size_t>(opt[0][0]);
        
        Fibre::Tractlet::Set tracts(input_location);
        Fibre::Tractlet::Set output_tracts(0, tracts.degree());
        output_tracts.add_elem_prop(Fibre::Base::Object::ALPHA_PROP);
        
        for (size_t tract_i = 0; tract_i < tracts.size(); ++tract_i) {
            
            if (!include.size() || std::find(include.begin(), include.end(), tract_i)
                    != include.end()) {
                
                const Fibre::Tractlet& tract = tracts[tract_i];
                
                Fibre::Strand::Set strands = tract.to_strands(num_strands);
                strands.freeze_elem_degree();
                Fibre::Strand::Set pos(0, strands.degree()), neg(0, strands.degree());
                pos.add_elem_prop(Fibre::Base::Object::ALPHA_PROP);
                neg.add_elem_prop(Fibre::Base::Object::ALPHA_PROP);
                
                Triple<double> centre_midpoint = tract.backbone().midpoint();
                
                for (size_t strand_i = 0; strand_i < strands.size(); ++strand_i) {
                    
                    Triple<double> disp = strands[strand_i].midpoint();
                    disp -= centre_midpoint;
                    
                    if (disp.dot(orient) >= 0)
                        pos.push_back(strands[strand_i]);
                    else
                        neg.push_back(strands[strand_i]);
                    
                }
                
                Fibre::Tractlet::Set pos_tract = pos.to_tractlets((size_t) 1);
                Fibre::Tractlet::Set neg_tract = neg.to_tractlets((size_t) 1);
                
                double acs_value = tract.acs() / 2.0;
                pos_tract[0].set_acs(acs_value);
                neg_tract[0].set_acs(acs_value);
                
                output_tracts.push_back(pos_tract[0]);
                output_tracts.push_back(neg_tract[0]);
                
            }
        }
        
        output_tracts.save(output_location);
        
    }
    
