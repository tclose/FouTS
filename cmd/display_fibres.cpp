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

#include "bts/common.h"

#include "bts/fibre/strand/set.h"

#include "bts/fibre/strand/set.h"
#include "bts/fibre/track/set.h"
#include "bts/fibre/tractlet/set.h"
#include "bts/fibre/tractlet/section.h"
#include "bts/fibre/base/set_reader.h"
#include "bts/fibre/base/set_writer.h"

#include "bts/inline_functions.h"

using namespace BTS;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Print strands to screen",
    "",
    NULL
};

ARGUMENTS= {
    Argument ("input", "The strands to be displayed.").allow_multiple().type_file (),
    Argument()
};

OPTIONS= {

    Option()};

EXECUTE {
    
        for (size_t arg_i = 0; arg_i < argument.size(); ++arg_i) {
            
            std::string input_location = argument[arg_i];
            
            if (input_location.at(input_location.length() - 1) == 'x')
                continue;
            
            if (File::has_extension<Fibre::Strand>(input_location)) {
                
                Fibre::Strand::Set s_set(input_location);
                
                std::cout << "\n\n\n**************************************\n";
                std::cout << "Strands at: '" << input_location << "'.\n";
                std::cout << "**************************************\n" << std::endl;
                
                std::cout << s_set;
                
            } else if (File::has_extension<Fibre::Track>(input_location)) {
                
                Fibre::Track::Set k_set(input_location);
                
                std::cout << "\n\n\n**************************************\n";
                std::cout << "Tracks at: '" << input_location << "'.\n";
                std::cout << "**************************************\n" << std::endl;
                
                std::cout << k_set;
                
            } else if (File::has_extension<Fibre::Tractlet>(input_location)) {
                
                Fibre::Tractlet::Set t_set(input_location);
                
                std::cout << "\n\n\n**************************************\n";
                std::cout << "Tractlets at: '" << input_location << "'.\n";
                std::cout << "**************************************\n" << std::endl;
                
                std::cout << t_set;
                
            } else if (File::has_extension<Fibre::Strand::Set>(input_location)) {
                
                Fibre::Strand::Set::Reader reader(input_location);
                
                Fibre::Strand::Set s_set;
                
                std::cout << "\n\n\n**************************************\n";
                std::cout << "Strand Sets at: '" << input_location << "'.\n";
                std::cout << "**************************************\n" << std::endl;
                
                size_t count = 0;
                
                while (reader.next(s_set)) {
                    
                    std::cout << "**** STRAND SET " << str(count++) << " ***" << std::endl;
                    
                    std::map<std::string, std::string> extend_props = s_set.get_extend_props();
                    
                    for (std::map<std::string, std::string>::iterator prop_it =
                            extend_props.begin(); prop_it != extend_props.end(); ++prop_it)
                        std::cout << prop_it->first << ": " << prop_it->second << ", ";
                    
                    std::cout << std::endl << std::endl;
                    
                    std::cout << std::endl << s_set << std::endl << std::endl;
                    
                }
                
            } else if (File::has_extension<Fibre::Tractlet::Set>(input_location)) {
                
                Fibre::Tractlet::Set::Reader reader(input_location);
                
                Fibre::Tractlet::Set t_set;
                
                std::cout << "\n\n\n**************************************\n";
                std::cout << "Tractlet Sets at: '" << input_location << "'.\n";
                std::cout << "**************************************\n" << std::endl;
                
                size_t count = 0;
                
                std::map<std::string, std::string> prop_row;
                
                while (reader.next(t_set)) {
                    
                    std::cout << "**** TRACTLET SET " << str(count++) << " ***" << std::endl;
                    
                    std::cout << std::endl << std::endl;
                    
                    std::map<std::string, std::string> extend_props = t_set.get_extend_props();
                    
                    for (std::map<std::string, std::string>::iterator prop_it =
                            extend_props.begin(); prop_it != extend_props.end(); ++prop_it)
                        std::cout << prop_it->first << ": " << prop_it->second << ", ";
                    
                    std::cout << std::endl << std::endl;
                    
                    std::cout << t_set << std::endl << std::endl;
                    
                }
                
            } else if (File::has_extension<Fibre::Track::Set>(input_location)) {
                
                Fibre::Track::Set::Reader reader(input_location);
                
                Fibre::Track::Set k_set;
                
                std::cout << "\n\n\n**************************************\n";
                std::cout << "Track Sets at: '" << input_location << "'.\n";
                std::cout << "**************************************\n" << std::endl;
                
                size_t count = 0;
                
                std::map<std::string, std::string> prop_row;
                
                while (reader.next(k_set)) {
                    
                    std::cout << "**** TRACK SET " << str(count++) << " ***" << std::endl;
                    
                    std::cout << std::endl << std::endl;
                    
                    std::map<std::string, std::string> extend_props = k_set.get_extend_props();
                    
                    for (std::map<std::string, std::string>::iterator prop_it =
                            extend_props.begin(); prop_it != extend_props.end(); ++prop_it)
                        std::cout << prop_it->first << ": " << prop_it->second << ", ";
                    
                    std::cout << std::endl << std::endl;
                    
                    std::cout << k_set << std::endl << std::endl;
                    
                }
                
            } else if (File::has_extension<Fibre::Strand::Section>(input_location)) {
                
                Fibre::Strand::Section::Reader reader(input_location);
                
                Fibre::Strand::Section sec;
                
                std::cout << "\n\n\n**************************************\n";
                std::cout << "Strand Sections at: '" << input_location << "'.\n";
                std::cout << "**************************************\n" << std::endl;
                
                while (reader.next(sec)) {
                    
                    std::cout << sec << std::endl;
                    
                }
                
            } else if (File::has_extension<Fibre::Tractlet::Section>(input_location)) {
                
                Fibre::Tractlet::Section::Reader reader(input_location);
                
                Fibre::Tractlet::Section sec;
                
                std::cout << "\n\n\n**************************************\n";
                std::cout << "Tractlet Sections at: '" << input_location << "'.\n";
                std::cout << "**************************************\n" << std::endl;
                
                while (reader.next(sec)) {
                    
                    std::cout << sec << std::endl;
                    
                }
                
            } else {
                if (argument.size() == 1)
                    throw Exception(
                            "Unrecognised file extension, \"" + File::extension(input_location)
                            + "\", for input file " + input_location + ".");
                else
                    std::cout << "Skipped '" << input_location << "'." << std::endl;
            }
        }
        
    }
