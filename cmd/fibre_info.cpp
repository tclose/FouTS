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

#include "bts/inline_functions.h"
#include "bts/fibre/strand.h"

using namespace FTS;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Displays fibres file header.",
    "",
    NULL
};

ARGUMENTS= {
    Argument ("fibres_file", "the input track file.").allow_multiple().type_file (),
    Argument()
};

OPTIONS= {

    Option ("property", "Print only the selected properties.")
    + Argument ("property", "").allow_multiple(),

    Option()
};

EXECUTE {
    
        std::vector<std::string> selected_properties;
        
        Options opt = get_options("property");
        if (opt.size())
            for (size_t prop_i = 0; prop_i < opt.size(); ++prop_i)
                selected_properties.push_back(opt[prop_i][0]);
        
        for (std::vector<ParsedArgument>::iterator arg = argument.begin(); arg != argument.end();
                ++arg) {
            Fibre::Strand::Reader reader;    // Doesn't really matter if it is a Strand/Tractlet or Track
            reader.open(arg->c_str());
            
            if (selected_properties.size()) {
                
                for (std::vector<std::string>::iterator select_it = selected_properties.begin();
                        select_it != selected_properties.end(); ++select_it)
                    std::cout << reader.get_extend_props()[*select_it] << " ";
                
                std::cout << std::endl;
                
            } else {
                
                std::cout << "***********************************\n";
                std::cout << "  Fibre file: \"" << arg->c_str() << "\"\n";
                std::map<std::string, std::string> extended_props = reader.get_extend_props();
                for (std::map<std::string, std::string>::map::iterator i = extended_props.begin();
                        i != extended_props.end(); ++i) {
                    std::string str(i->first + ':');
                    str.resize(40, ' ');
                    std::cout << "    " << str << i->second << "\n";
                }
                
            }
            
        }
        
    }
