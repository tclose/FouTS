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

#include "bts/common.h"

#include "bts/fibre/strand/set.h"

#include "bts/fibre/strand/set.h"
#include "bts/fibre/track/set.h"
#include "bts/fibre/tractlet/set.h"
#include "bts/fibre/base/set_reader.h"
#include "bts/fibre/base/set_writer.h"

#include "bts/inline_functions.h"

using namespace FTS;

const double SCALE_DEFAULT = 1.0;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Scale the size of all fibres",
    "",
    NULL
};

ARGUMENTS= {
    Argument ("input", "The fibres to be trimmed.").type_file (),
    Argument ("output", "The returned trimmed fibres.").optional().type_file(),
    Argument()
};

OPTIONS= {

    Option ("scale", "The scalar the fibres will be multiplied by.")
    + Argument ("scale", "").type_float (SMALL_FLOAT, SCALE_DEFAULT, LARGE_FLOAT),

    Option ("include_acs", "Scale the acs and base intensities as well as the rest of the strand."),

    Option()};

EXECUTE {
    
        std::string input_location = argument[0];
        std::string output_location;
        
        if (argument.size() > 1)
            output_location = argument[1].c_str();
        else
            output_location = input_location;
        
        double scale = SCALE_DEFAULT;
        bool include_acs = false;
        
        Options opt = get_options("scale");
        if (opt.size())
            scale = opt[0][0];
        
        opt = get_options("include_acs");
        if (opt.size())
            include_acs = true;
        
        if (File::has_extension<Fibre::Strand>(input_location)) {
            
            Fibre::Strand::Set strands(input_location);
            
            Fibre::Strand::Set output = strands * scale;
            
            if (!include_acs) {
                
                output.set_base_intensity(strands.base_intensity());
                
                for (size_t strand_i = 0; strand_i < strands.size(); ++strand_i)
                    output[strand_i].set_acs(strands[strand_i].acs());
                
            }
            
            output.save(output_location);
            
        } else if (File::has_extension<Fibre::Track>(input_location)) {
            
            Fibre::Track::Set tcks(input_location);
            
            Fibre::Track::Set output = tcks * scale;
            
            if (!include_acs) {
                
                output.set_base_intensity(tcks.base_intensity());
                
                for (size_t tck_i = 0; tck_i < tcks.size(); ++tck_i)
                    output[tck_i].set_acs(tcks[tck_i].acs());
                
            }
            
            output.save(output_location);
            
        } else if (File::has_extension<Fibre::Tractlet>(input_location)) {
            
            Fibre::Tractlet::Set tracts(input_location);
            
            Fibre::Tractlet::Set output = tracts * scale;
            
            if (!include_acs) {
                
                output.set_base_intensity(tracts.base_intensity());
                
                for (size_t tract_i = 0; tract_i < tracts.size(); ++tract_i)
                    output[tract_i].set_acs(tracts[tract_i].acs());
                
            }
            
            output.save(output_location);
            
        } else
            throw Exception("Unrecognised file extension of input file '" + input_location + "'.");
        
    }
