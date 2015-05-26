/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by J-Donald Tournier, 27/06/08.

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

#include "bts/cmd.h"
#include "bts/common.h"

#include "bts/fibre/strand/set.h"
#include "bts/fibre/tractlet/set.h"
#include "bts/fibre/track/set.h"

#include "bts/inline_functions.h"

using namespace FTS;

template<typename T> void reset_indices(T& input_location);
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Resets bundle indices to a consecutive range",
    "",
    NULL
};

ARGUMENTS= {
    Argument ("input file", "The input fibres file location.").type_file (),
    Argument ("output file", "The output fibres file location.").optional().type_file (),
    Argument()
};

OPTIONS= {

    Option()
};

EXECUTE {
    
        std::string input_location = argument[0];
        std::string output_location = input_location;
        
        if (argument.size() > 1)
            output_location = argument[1].c_str();
        
        if (File::has_extension<Fibre::Strand>(input_location)) {
            
            Fibre::Strand::Set strands(input_location);
            
            reset_indices(strands);
            
            strands.save(output_location);
            
        } else if (File::has_extension<Fibre::Tractlet>(input_location)) {
            
            Fibre::Tractlet::Set tractlets(input_location);
            
            reset_indices(tractlets);
            
            tractlets.save(output_location);
            
        } else if (File::has_extension<Fibre::Track>(input_location)) {
            
            Fibre::Track::Set tcks(input_location);
            
            reset_indices(tcks);
            
            tcks.save(output_location);
            
        } else
            throw Exception("Unrecognised extension of input file '" + input_location + "'.");
        
    }
    
    template<typename T> void reset_indices(T& fibres) {
        
        if (!fibres.has_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP))
            fibres.add_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP, "-1");
        
        std::map<std::string, size_t> indices;
        
        size_t bundle_count = 0;
        
        for (size_t strand_i = 0; strand_i < fibres.size(); ++strand_i) {
            
            std::string bundle_index = fibres.get_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP,
                    strand_i);
            
            if (!indices.count(bundle_index))
                indices[bundle_index] = bundle_count++;
            
        }
        
        for (size_t strand_i = 0; strand_i < fibres.size(); ++strand_i) {
            
            size_t reset_bundle_index = indices[fibres.get_extend_elem_prop(
                    Fibre::Track::BUNDLE_INDEX_EPROP, strand_i)];
            
            fibres.set_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP, str(reset_bundle_index),
                    strand_i);
            
        }
        
    }
