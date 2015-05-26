/*
 Copyright 2010 Brain Research Institute, Melbourne, Australia

 Written by Thomas G. Close, Apr 21, 2010.

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

#include "bts/file.h"

#include "progressbar.h"

#include "bts/fibre/strand/set.h"
#include "bts/fibre/tractlet/set.h"
#include "bts/fibre/track/set.h"
#include "bts/fibre/base/set_reader.h"
#include "bts/fibre/base/set_writer.h"

#include "bts/inline_functions.h"

using namespace FTS;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "set_extend_elem_prop.cpp",
    "",
    NULL
};

ARGUMENTS= {

    Argument ("input", "Input file location").type_file (),
    Argument ("output", "Output file location").optional().type_file (),

    Argument()};

OPTIONS= {

    Option ("set", "Set a property to a new value").allow_multiple()
    + Argument ("name", "Name of property")
    + Argument ("value", "Value the property will be set to."),

    Option ("set_elem", "Set a property to a new value").allow_multiple()
    + Argument ("name", "Name of property")
    + Argument ("value", "Value the property will be set to."),

    Option()};

template<typename T> void set_properties(const std::string& input_location,
                                         const std::string& output_location, App::Options& opt) {
    
    typename T::Set fibres(input_location);
    
    for (size_t prop_i = 0; prop_i < opt.size(); ++prop_i) {
        
        std::string property_name = opt[prop_i][0];
        std::string property_value = opt[prop_i][1];
        
        if (property_name == "base_intensity")
            fibres.set_base_intensity(to<double>(property_value));
        else
            fibres.set_extend_prop(property_name, property_value);
        
    }
    
    fibres.save(output_location);
    
}

template<typename T> void set_set_properties(const std::string& input_location,
                                             const std::string& output_location,
                                             App::Options& opt) {
    
    if (input_location == output_location)
        throw Exception("Input and output locations cannot be the same for fibre sets.");
    
    typename T::Reader reader(input_location);
    
    std::map<std::string, std::string> properties = reader.get_extend_props();
    
    for (size_t prop_i = 0; prop_i < opt.size(); ++prop_i) {
        
        std::string property_name = opt[prop_i][0];
        std::string property_value = opt[prop_i][1];
        
        properties[property_name] = property_value;
        
    }
    
    typename T::Writer writer(output_location, reader, properties);
    
    T strands;
    
    while (reader.next(strands))
        writer.append(strands);
    
}

template<typename T> void set_elem_properties(const std::string& input_location,
                                              const std::string& output_location, App::Options& opt) {

    typename T::Set fibres(input_location);

    for (size_t prop_i = 0; prop_i < opt.size(); ++prop_i) {

        std::string property_name = opt[prop_i][0];
        std::string property_value = opt[prop_i][1];

        if (property_name == "acs" && !fibres.has_elem_prop(Fibre::Tractlet::ALPHA_PROP))
            fibres.add_elem_prop(Fibre::Tractlet::ALPHA_PROP);

        for (size_t fibre_i = 0; fibre_i < fibres.size(); ++fibre_i) {
            if (property_name == "acs") {
                fibres[fibre_i].set_acs(to<double>(property_value));
            } else
                fibres.set_extend_elem_prop(property_name, property_value, fibre_i);
        }
    }

    fibres.save(output_location);

}

EXECUTE {
    
        std::string input_location = argument[0];
        std::string output_location;
        
        if (argument.size() > 1)
            output_location = argument[1].c_str();
        else
            output_location = input_location;
        
        MR::ProgressBar progress_bar("Setting properties...");
        
        Options opt = get_options("set");
        
        if (opt.size()) {
            
            if (File::has_or_txt_extension<Fibre::Strand>(input_location)) {

                set_properties<Fibre::Strand>(input_location, output_location, opt);

            } else if (File::has_or_txt_extension<Fibre::Tractlet>(input_location)) {

                set_properties<Fibre::Tractlet>(input_location, output_location, opt);

            } else if (File::has_or_txt_extension<Fibre::Track>(input_location)) {

                set_properties<Fibre::Track>(input_location, output_location, opt);

            } else if (File::has_extension<Fibre::Strand::Set>(input_location)) {

                set_set_properties<Fibre::Strand::Set>(input_location, output_location, opt);

            } else if (File::has_extension<Fibre::Tractlet::Set>(input_location)) {

                set_set_properties<Fibre::Tractlet::Set>(input_location, output_location, opt);

            } else if (File::has_extension<Fibre::Track::Set>(input_location)) {

                set_set_properties<Fibre::Track::Set>(input_location, output_location, opt);

            } else
                throw Exception("Unrecognised extension '" + File::extension(input_location) + "'.");
        }

        opt = get_options("set_elem");

        if (opt.size()) {

            if (File::has_or_txt_extension<Fibre::Strand>(input_location)) {

                set_elem_properties<Fibre::Strand>(input_location, output_location, opt);

            } else if (File::has_or_txt_extension<Fibre::Tractlet>(input_location)) {

                set_elem_properties<Fibre::Tractlet>(input_location, output_location, opt);

            } else if (File::has_or_txt_extension<Fibre::Track>(input_location)) {

                set_elem_properties<Fibre::Track>(input_location, output_location, opt);

            } else
                throw Exception("Unrecognised extension '" + File::extension(input_location) + "'.");

        }

        //MR::ProgressBar::done();
    }
    
