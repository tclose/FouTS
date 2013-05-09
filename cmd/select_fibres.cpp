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

#include "bts/fibre/base/writer.h"

#include "bts/image/expected/buffer.h"

#include "bts/inline_functions.h"

#define ONE_TO_ONE(From,To) \
  (File::has_or_txt_extension<Fibre::To>(output_location)) { \
    one_to_one<Fibre::From, Fibre::To>(input, input_location, output_location, include, degree, num_length_sections, num_tractlets, width_of_tractlet, strands_per_acs, num_width_sections); \
  } \

#define ONE_TO_SET(From,To) \
  (File::has_extension<Fibre::To::Set>(output_location)) { \
    Fibre::To::Set::Writer writer(output_location, input); \
    one_to_set<Fibre::From, Fibre::To::Set, Fibre::To::Set::Writer>(input, input_location, output_location, writer, include, degree, num_length_sections, num_tractlets, width_of_tractlet, strands_per_acs, num_width_sections, separate_sets); \
  }

#define ONE_TO_TXT_SET(From,To) \
  (File::has_txt_extension<Fibre::To::Set>(output_location)) { \
    Fibre::Base::SetTextWriter<Fibre::To::Set> writer(output_location, input); \
    one_to_set<Fibre::From, Fibre::To::Set, Fibre::Base::SetTextWriter<Fibre::To::Set> >(input, input_location, output_location, writer, include, degree, num_length_sections, num_tractlets, width_of_tractlet, strands_per_acs, num_width_sections, separate_sets); \
  }

#define SET_TO_SET(From,To) \
  (File::has_extension<Fibre::To::Set>(output_location)) { \
    add_extra_element_properties<Fibre::From, Fibre::To>(elem_header); \
    Fibre::To::Set::Writer writer(output_location, reader, reader.extend_prop_keys(), elem_header, reader.get_extend_props()); \
    set_to_set<Fibre::From::Set, Fibre::To::Set, Fibre::From::Set::Reader, Fibre::To::Set::Writer>(input_location, output_location, reader, writer, include, degree, num_length_sections, num_tractlets, width_of_tractlet, strands_per_acs, num_width_sections); \
  }

#define SET_TO_ONE(From,To) \
  (File::has_or_txt_extension<Fibre::To>(output_location)) { \
    set_to_one<Fibre::From::Set, Fibre::To, Fibre::From::Set::Reader>(input_location, output_location, reader, include, degree, num_length_sections, num_tractlets, width_of_tractlet, strands_per_acs, num_width_sections); \
  } \

#define SET_TO_TXT_SET(From,To) \
  (File::has_txt_extension<Fibre::To::Set>(output_location)) { \
    add_extra_element_properties<Fibre::From, Fibre::To>(elem_header); \
    Fibre::Base::SetTextWriter<Fibre::To::Set> writer(output_location, reader, reader.extend_prop_keys(), elem_header, reader.get_extend_props()); \
    set_to_set<Fibre::From::Set, Fibre::To::Set, Fibre::From::Set::Reader, Fibre::Base::SetTextWriter<Fibre::To::Set> >(input_location, output_location, reader, writer, include, degree, num_length_sections, num_tractlets, width_of_tractlet, strands_per_acs, num_width_sections); \
  }

#define TXT_SET_TO_ONE(From,To) \
  (File::has_or_txt_extension<Fibre::To>(output_location)) { \
    set_to_one<Fibre::From::Set, Fibre::To, Fibre::Base::SetTextReader<Fibre::From::Set> >(input_location, output_location, reader, include, degree, num_length_sections, num_tractlets, width_of_tractlet, strands_per_acs, num_width_sections); \
  }

#define TXT_SET_TO_SET(From,To) \
  (File::has_extension<Fibre::To::Set>(output_location)) { \
    add_extra_element_properties<Fibre::From, Fibre::To>(elem_header); \
    Fibre::To::Set::Writer writer(output_location, reader, reader.extend_prop_keys(), elem_header, reader.get_extend_props()); \
    set_to_set<Fibre::From::Set, Fibre::To::Set, Fibre::Base::SetTextReader<Fibre::From::Set> , Fibre::To::Set::Writer>(input_location, output_location, reader, writer, include, degree, num_length_sections, num_tractlets, width_of_tractlet, strands_per_acs, num_width_sections); \
  }

#define TXT_SET_TO_TXT_SET(From,To) \
  (File::has_txt_extension<Fibre::To::Set>(output_location)) { \
    add_extra_element_properties<Fibre::From, Fibre::To>(elem_header); \
    Fibre::Base::SetTextWriter<Fibre::To::Set> writer(output_location, reader, reader.extend_prop_keys(), elem_header, reader.get_extend_props()); \
    set_to_set<Fibre::From::Set, Fibre::To::Set, Fibre::Base::SetTextReader<Fibre::From::Set> , Fibre::Base::SetTextWriter<Fibre::To::Set> >(input_location, output_location, reader, writer, include, degree, num_length_sections, num_tractlets, width_of_tractlet, strands_per_acs, num_width_sections); \
  }

#define ONE_CONVERSIONS(From) \
  (File::has_or_txt_extension<Fibre::From>(input_location)) { \
      Fibre::From::Set input(input_location); \
\
    if \
      ONE_TO_ONE(From,Strand) \
    else if \
      ONE_TO_ONE(From,Tractlet) \
    else if \
      ONE_TO_ONE(From,Track) \
    else if \
      ONE_TO_SET(From,Strand) \
    else if \
      ONE_TO_SET(From,Tractlet) \
    else if \
      ONE_TO_SET(From,Track) \
    else if \
      ONE_TO_TXT_SET(From,Strand) \
    else if \
      ONE_TO_TXT_SET(From,Tractlet) \
    else if \
      ONE_TO_TXT_SET(From,Track) \
    else \
      throw Exception ("Unsupported conversion from '" + File::extension(input_location) + "' to '" + File::extension(output_location) + "'."); \
  }

#define SET_CONVERSIONS(From) \
  (File::has_extension<Fibre::From::Set>(input_location)) { \
    Fibre::From::Set::Reader reader(input_location); \
    std::vector<std::string> elem_header = reader.extend_elem_prop_keys(); \
\
    if \
      SET_TO_ONE(From,Strand) \
    else if \
      SET_TO_ONE(From,Tractlet) \
    else if \
      SET_TO_ONE(From,Track) \
    else if \
      SET_TO_SET(From,Strand) \
    else if \
      SET_TO_SET(From,Tractlet) \
    else if \
      SET_TO_SET(From,Track) \
    else if \
      SET_TO_TXT_SET(From,Strand) \
    else if \
      SET_TO_TXT_SET(From,Tractlet) \
    else if \
      SET_TO_TXT_SET(From,Track) \
    else \
      throw Exception ("Unsupported conversion from '" + File::extension(input_location) + "' to '" + File::extension(output_location) + "'."); \
  }

#define TXT_SET_CONVERSIONS(From) \
  (File::has_txt_extension<Fibre::From::Set>(input_location)) { \
    Fibre::Base::SetTextReader<Fibre::From::Set>  reader(input_location); \
    std::vector<std::string> elem_header = reader.extend_elem_prop_keys(); \
\
    if \
      TXT_SET_TO_ONE(From,Strand) \
    else if \
      TXT_SET_TO_ONE(From,Tractlet) \
    else if \
      TXT_SET_TO_ONE(From,Track) \
    else if \
      TXT_SET_TO_SET(From,Strand) \
    else if \
      TXT_SET_TO_SET(From,Tractlet) \
    else if \
      TXT_SET_TO_SET(From,Track) \
    else if \
      TXT_SET_TO_TXT_SET(From,Strand) \
    else if \
      TXT_SET_TO_TXT_SET(From,Tractlet) \
    else if \
      TXT_SET_TO_TXT_SET(From,Track) \
    else \
      throw Exception ("Unsupported conversion from '" + File::extension(input_location) + "' to '" + File::extension(output_location) + "'."); \
  }

using namespace BTS;
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Select a subset of strands",
    "",
    NULL
};

ARGUMENTS= {
    Argument ("input", "The strands to be selected from.").type_file (),
    Argument ("output", "The selected strands.").type_file (),
    Argument()
};

OPTIONS= {

    Option ("include","The indices of the strands to include")
    + Argument ("include","").type_text(),

    Option ("bundle_include","The indices of the strands to bundle_include")
    + Argument ("bundle_include","").type_text(),

    Option ("num_length_sections", "Number of samples to take along the strand path")
    + Argument ("num_length_sections", "").type_integer (1, 100, 2000),

    Option ("degree", "The degree of the Strand coefficients used to describe the strands")
    + Argument ("degree", "").type_integer (1, Fibre::Strand::DEFAULT_DEGREE, LARGE_INT),

    Option ("num_tractlets", "Used when converting between strands and tractlets when the strands are to first clustered into bundles then converted tractlets.")
    + Argument ("num_tractlets", "").type_integer (0, -1.0, LARGE_INT),

    Option ("width_of_tractlet","Used when converting between strands and tractlets when the strands are to represent the backbone of the converted tractlets.")
    + Argument ("","").type_float (-1, -1, LARGE_FLOAT),

    Option ("strands_per_acs","Basewidth of the imported tractlets.")
    + Argument ("","").type_float (-1, Fibre::Tractlet::STRANDS_PER_AREA_DEFAULT, LARGE_FLOAT),

    Option ("num_width_sections","Number of samples across the tractlet radius.")
    + Argument ("num_width_sections","").type_integer (0, Image::Expected::Buffer::NUM_WIDTH_SECTIONS_DEFAULT, LARGE_INT),

    Option ("separate_sets","Separate selected strands into separate sets."),

    Option()};

//FIXME: Move this into file.h once all the dependencies have been sorted out.

std::string class_name(const std::string& location) {
    
    std::string name;
    
    if (File::has_or_txt_extension<Fibre::Tractlet>(location))
        name = "tractlet";
    else if (File::has_or_txt_extension<Fibre::Strand>(location))
        name = "strand";
    else if (File::has_or_txt_extension<Fibre::Track>(location))
        name = "track";
    else if (File::has_or_txt_extension<Fibre::Tractlet::Set>(location))
        name = "tractlet set";
    else if (File::has_or_txt_extension<Fibre::Strand::Set>(location))
        name = "strand set";
    else if (File::has_or_txt_extension<Fibre::Track::Set>(location))
        name = "track set";
    else
        throw Exception("Unknown extension of file '" + location + "'.");
    
    return name;
    
}

template<typename InputType, typename OutputType> OutputType& convert(const InputType& input,
                                                                      OutputType& output,
                                                                      size_t degree,
                                                                      size_t num_length_sections,
                                                                      size_t num_tractlets,
                                                                      double width_of_tractlet,
                                                                      double strands_per_acs,
                                                                      size_t num_width_sections);

template<typename InputType, typename OutputType> void add_extra_element_properties(
        std::vector<std::string>& elem_header);

template<> Fibre::Strand::Set& convert(const Fibre::Track::Set& input, Fibre::Strand::Set& output,
                                       size_t degree, size_t num_length_sections,
                                       size_t num_tractlets, double width_of_tractlet,
                                       double strands_per_acs, size_t num_width_sections) {
    
    return output = input.to_strands(degree);
    
}

template<> Fibre::Track::Set& convert(const Fibre::Strand::Set& input, Fibre::Track::Set& output,
                                      size_t degree, size_t num_length_sections,
                                      size_t num_tractlets, double width_of_tractlet,
                                      double strands_per_acs, size_t num_width_sections) {
    
    return output = input.to_tracks(num_length_sections);
    
}

template<> Fibre::Tractlet::Set& convert(const Fibre::Strand::Set& input,
                                         Fibre::Tractlet::Set& output, size_t degree,
                                         size_t num_length_sections, size_t num_tractlets,
                                         double width_of_tractlet, double strands_per_acs,
                                         size_t num_width_sections) {
    
    if (num_tractlets)
        output = input.to_tractlets(num_tractlets);
    //Requires bundle index property to be present in strand set.
    else if (width_of_tractlet > 0)
        output = input.to_tractlets(width_of_tractlet);
    //Each strand is treated as the primary strand of a tractlet with base width 'width_of_tractlet'.
    else
        output = input.to_tractlets();
    
    return output;
    
}

template<> Fibre::Strand::Set& convert(const Fibre::Tractlet::Set& input,
                                       Fibre::Strand::Set& output, size_t degree,
                                       size_t num_length_sections, size_t num_tractlets,
                                       double width_of_tractlet, double strands_per_acs,
                                       size_t num_width_sections) {
    
    if (strands_per_acs > 0)
        output = input.to_strands(strands_per_acs);
    else if (num_width_sections > 0)
        output = input.to_strands(num_width_sections);
    else
        throw Exception("either strands_per_acs or num_width_sections needs to be set.");
    
    return output;
    
}

template<> Fibre::Track::Set& convert(const Fibre::Tractlet::Set& input, Fibre::Track::Set& output,
                                      size_t degree, size_t num_length_sections,
                                      size_t num_tractlets, double width_of_tractlet,
                                      double strands_per_acs, size_t num_width_sections) {
    
    Fibre::Strand::Set strands;
    
    if (strands_per_acs > 0)
        strands = input.to_strands(strands_per_acs);
    else if (num_width_sections > 0)
        strands = input.to_strands(num_width_sections);
    else
        throw Exception("either strands_per_acs or num_width_sections needs to be set.");
    
    return output = strands.to_tracks(num_length_sections);
    
}

template<> Fibre::Tractlet::Set& convert(const Fibre::Track::Set& input,
                                         Fibre::Tractlet::Set& output, size_t degree,
                                         size_t num_length_sections, size_t num_tractlets,
                                         double width_of_tractlet, double strands_per_acs,
                                         size_t num_width_sections) {
    
    if (num_tractlets)
        output = input.to_tractlets(degree, num_tractlets);
    //Requires bundle index property to be present in strand set.
    else if (width_of_tractlet > 0)
        output = input.to_tractlets(degree, width_of_tractlet);
    //Each strand is treated as the primary strand of a tractlet with base width 'width_of_tractlet'.
    else
        output = input.to_tractlets(degree);
    
    return output;
    
}

template<typename InputType, typename OutputType> OutputType& convert(const InputType& input,
                                                                      OutputType& output,
                                                                      size_t degree,
                                                                      size_t num_length_sections,
                                                                      size_t num_tractlets,
                                                                      double width_of_tractlet,
                                                                      double strands_per_acs,
                                                                      size_t num_width_sections) {
    
    return output = OutputType(input);
    
}

template<> void add_extra_element_properties<Fibre::Tractlet, Fibre::Strand>(
        std::vector<std::string>& header) {
    
    if (std::find(header.begin(), header.end(), Fibre::Track::BUNDLE_INDEX_EPROP) == header.end())
        header.push_back(Fibre::Track::BUNDLE_INDEX_EPROP);
    
    std::vector<std::string>::iterator prop_it = std::find(header.begin(), header.end(),
            str("tract_volume"));
    
    if (prop_it != header.end())
        header.erase(prop_it);
    
}

template<> void add_extra_element_properties<Fibre::Tractlet, Fibre::Track>(
        std::vector<std::string>& header) {
    
    if (std::find(header.begin(), header.end(), Fibre::Track::BUNDLE_INDEX_EPROP) == header.end())
        header.push_back(Fibre::Track::BUNDLE_INDEX_EPROP);
    
    std::vector<std::string>::iterator prop_it = std::find(header.begin(), header.end(),
            str("tract_volume"));
    
    if (prop_it != header.end())
        header.erase(prop_it);
    
}

template<typename InputType, typename OutputType> void add_extra_element_properties(
        std::vector<std::string>& elem_header) { /* do nothing */
}

template<typename InputType, typename OutputType>
void one_to_one(const typename InputType::Set& input, const std::string& input_location,
                const std::string& output_location, std::vector<size_t>& include, size_t degree,
                size_t num_length_sections, size_t num_tractlets, double width_of_tractlet,
                double strands_per_acs, size_t num_width_sections) {
    
    typename InputType::Set selected;
    typename OutputType::Set output(input.get_extend_props());
    
    if (include.size())
        input.select(selected, include);
    else
        selected = input;
    
    MR::ProgressBar progress_bar(
            "Selecting from " + class_name(input_location) + " set and outputting to "
            + class_name(output_location) + " set...");
    
    convert(selected, output, degree, num_length_sections, num_tractlets, width_of_tractlet,
            strands_per_acs, num_width_sections);
    
    output.save(output_location);
    
}

template<typename InputType, typename OutputType, typename ReaderType, typename WriterType>
void set_to_set(const std::string& input_location, const std::string& output_location,
                ReaderType& reader, WriterType& writer, std::vector<size_t>& include, size_t degree,
                size_t num_length_sections, size_t num_tractlets, double width_of_tractlet,
                double strands_per_acs, size_t num_width_sections) {
    
    InputType input_set;
    
    size_t index = 0;
    
    size_t num_included_sets;
    
    if (include.size())
        num_included_sets = include.size();
    else
        num_included_sets = to<size_t>(reader.get_extend_props()["count"]);
    
    MR::ProgressBar progress_bar(
            "Selecting from set of " + class_name(input_location) + "s and outputting to set of "
            + class_name(output_location) + "s...", num_included_sets);
    
    while (reader.next(input_set)) {
        
        if (!include.size() || std::find(include.begin(), include.end(), index) != include.end()) {
            
            OutputType output_set;
            
            convert(input_set, output_set, degree, num_length_sections, num_tractlets,
                    width_of_tractlet, strands_per_acs, num_width_sections);
            
            writer.append(output_set);
            
            progress_bar++;
        }
        
        ++index;
    }
    
}

template<typename InputType, typename OutputType, typename ReaderType>
void set_to_one(const std::string& input_location, const std::string& output_location,
                ReaderType& reader, std::vector<size_t>& include, size_t degree,
                size_t num_length_sections, size_t num_tractlets, double width_of_tractlet,
                double strands_per_acs, size_t num_width_sections) {
    
    InputType input_set;
    
    if (!reader.next(input_set))
        throw Exception("No sets found in file '" + input_location + "'.");
    
    reader.rewind();
    
    std::map<std::string, std::string> props; //= input_set.get_extend_props();
    input_set.insert_props(props);

    std::vector<std::string> header;
    
    //Add extra element properties that might be generated in the conversion process.
    add_extra_element_properties<typename InputType::Element, OutputType>(header);
    
    Fibre::Base::Writer<OutputType> writer(output_location,
            Fibre::Base::Object::select_props<OutputType>(reader.elem_prop_keys()), header, props);
    
    size_t index = 0;
    
    size_t num_included_sets;
    
    if (include.size())
        num_included_sets = include.size();
    else
        num_included_sets = to<size_t>(reader.get_extend_props()["count"]);
    
    MR::ProgressBar progress_bar(
            "Selecting from set of " + class_name(input_location) + "s and outputting to "
            + class_name(output_location) + " set...", num_included_sets);
    
    while (reader.next(input_set)) {
        
        if (!include.size() || std::find(include.begin(), include.end(), index) != include.end()) {
            
            typename OutputType::Set output_set;
            
            convert(input_set, output_set, degree, num_length_sections, num_tractlets,
                    width_of_tractlet, strands_per_acs, num_width_sections);
            
            for (size_t fibre_i = 0; fibre_i < output_set.size(); ++fibre_i)
                writer.append(output_set[fibre_i], output_set.get_extend_elem_prop_row(fibre_i));
            
            progress_bar++;
        }
        
        ++index;
    }
    
}

template<typename InputType, typename OutputType, typename WriterType>
void one_to_set(const typename InputType::Set& input, const std::string& input_location,
                const std::string& output_location, WriterType& writer,
                std::vector<size_t>& include, size_t degree, size_t num_length_sections,
                size_t num_tractlets, double width_of_tractlet, double strands_per_acs,
                size_t num_width_sections, bool separate_sets) {
    
    typename InputType::Set selected;
    OutputType converted;
    
    if (include.size())
        input.select(selected, include);
    else
        selected = input;
    
    convert(selected, converted, degree, num_length_sections, num_tractlets, width_of_tractlet,
            strands_per_acs, num_width_sections);
    
    if (separate_sets) {
        
        MR::ProgressBar progress_bar(
                "Selecting from " + class_name(input_location) + " set and outputting into set of "
                + class_name(output_location) + "s, with one fibre per set...", selected.size());
        
        for (size_t elem_i = 0; elem_i < converted.size(); ++elem_i) {
            
            OutputType output_set;
            
            output_set.push_back(converted[elem_i]);
            
            writer.append(output_set);
            
            ++progress_bar;
            
        }
        
    } else {
        
        MR::ProgressBar progress_bar(
                "Selecting from " + class_name(input_location)
                + " set and outputting into a singleton set of " + class_name(output_location)
                + "s...");
        
        writer.append(converted);
        
    }
}

EXECUTE {
    
        std::string input_location = argument[0];
        std::string output_location = argument[1];
        
        std::vector<size_t> include;
        std::vector<size_t> bundle_include;
        size_t degree = Fibre::Strand::DEFAULT_DEGREE;
        size_t num_length_sections = Fibre::Track::NUM_LENGTH_SECTIONS_DEFAULT;
        size_t num_tractlets = 0;
        double width_of_tractlet = -1;
        double strands_per_acs = Fibre::Tractlet::STRANDS_PER_AREA_DEFAULT;
        size_t num_width_sections = 0;    //Image::Expected::Buffer::NUM_WIDTH_SECTIONS_DEFAULT;
        bool separate_sets = false;
        
        Options opt = get_options("include");
        if (opt.size())
            include = parse_sequence<size_t>(opt[0][0]);
        
        opt = get_options("bundle_include");
        if (opt.size())
            bundle_include = parse_sequence<size_t>(opt[0][0]);
        
        opt = get_options("num_length_sections");
        if (opt.size())
            num_length_sections = opt[0][0];
        
        opt = get_options("degree");
        if (opt.size())
            degree = opt[0][0];
        
        opt = get_options("num_tractlets");
        if (opt.size())
            num_tractlets = opt[0][0];
        
        opt = get_options("width_of_tractlet");
        if (opt.size())
            width_of_tractlet = opt[0][0];
        
        opt = get_options("num_width_sections");
        if (opt.size()) {
            num_width_sections = opt[0][0];
            strands_per_acs = -1.0;
        }
        
        opt = get_options("strands_per_acs");
        if (opt.size())
            strands_per_acs = opt[0][0];
        
        opt = get_options("separate_sets");
        if (opt.size())
            separate_sets = true;
        
        if (include.size() && bundle_include.size())
            throw Exception("'-include' and '-bundle_include' cannot be used simultaneously.");
        
        //If bundle_include option provided get bundle_index property and map to individual fibres.
        if (bundle_include.size()) {
            if (File::has_or_txt_extension<Fibre::Strand>(input_location) || File::has_or_txt_extension<
                        Fibre::Track>(input_location)) {
                
                //Although the input could also be a Track set, they can be loaded temporarily as strands first.
                Fibre::Strand::Set strands(input_location, (size_t) 1);
                
                if (!strands.has_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP))
                    throw Exception(
                            "'" + Fibre::Track::BUNDLE_INDEX_EPROP
                            + "' was not found in input strand set.");
                
                for (size_t strand_i = 0; strand_i < strands.size(); ++strand_i) {
                    
                    size_t bundle_index = to<size_t>(
                            strands.get_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP,
                                    strand_i));
                    
                    if (std::find(bundle_include.begin(), bundle_include.end(), bundle_index) != bundle_include.end())
                        include.push_back(strand_i);
                    
                }
                
            }
            
        }
        
        if ((strands_per_acs > 0) && num_width_sections)
            throw Exception(
                    "'-num_width_sections' and '-strands_per_acs' option cannot be used in conjunction.");
        
        else if ((strands_per_acs <= 0.0)
                && !num_width_sections
                && (File::has_txt_or_set_extension<Fibre::Tractlet>(input_location) && (File::has_txt_or_set_extension<
                                                                                                Fibre::Strand>(
                                                                                                output_location)
                                                                                        || File::has_txt_or_set_extension<
                                                                                                Fibre::Track>(
                                                                                                output_location))))
            throw Exception(
                    "Either '-num_width_sections' and '-strands_per_acs' must be supplied for conversion from tractlet to strand/track classes.");
        
        MR::ProgressBar progress_bar("Selecting fibres ...");
        
    if ONE_CONVERSIONS (Strand)
    else if
    ONE_CONVERSIONS(Tractlet)
    else if
    ONE_CONVERSIONS(Track)
    else if
    SET_CONVERSIONS(Strand)
    else if
    SET_CONVERSIONS(Tractlet)
    else if
    SET_CONVERSIONS(Track)
    else if
    TXT_SET_CONVERSIONS(Strand)
    else if
    TXT_SET_CONVERSIONS(Tractlet)
    else if
    TXT_SET_CONVERSIONS(Track)
    else
    throw Exception ("Unsupported input file '" + File::extension(input_location) + "'.");

}

