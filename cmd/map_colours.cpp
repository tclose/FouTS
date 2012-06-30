/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, Mar 14, 2011.

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

#include "bts/fibre/tractlet/set.h"
#include "bts/fibre/strand/set.h"
#include "bts/fibre/track/set.h"

#include "bts/inline_functions.h"

using namespace BTS;

SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);


DESCRIPTION = {
  "Map colours to reference set by switching bundle indices to closest match.",
  "",
  NULL
};

ARGUMENTS = {
  Argument ("reference", "Reference fibres the match fibres will be mapped to.").type_file (),
  Argument ("match", "Input fibres that will have their bundle indices mapped.").type_file (),
  Argument ("output", "Output fibres (if not provided will be same as match fibres)").optional().type_file (),
  Argument ("input", "Input fibres (if not provided will be same as match fibres)").optional().type_file (),
  Argument()
};


OPTIONS = {

    Option ("num_length_sections", "Number of samples to take along the strand path")
     + Argument ("num_length_sections", "").type_integer (1, 100, 2000),

    Option ("degree", "The degree of the Strand coefficients used to describe the strands")
     + Argument ("degree", "").type_integer (1, Fibre::Strand::DEFAULT_DEGREE, LARGE_INT),

Option() };


EXECUTE {

  std::string reference_location = argument[0];
  std::string match_location= argument[1];
	std::string output_location = match_location;
  std::string input_location = match_location;

	if (argument.size() == 3)
	  output_location = argument[2].c_str();

	if (argument.size() == 4) {
    output_location = argument[2].c_str();
    input_location = argument[3].c_str();
	}

  size_t degree = Fibre::Strand::DEFAULT_DEGREE;
//  size_t num_length_sections = 0;


  Options opt = get_options("degree");
  if (opt.size())
    degree = opt[0][0];

  //  Options opt = get_options("num_length_sections");
  //  if (opt.size())
  //    num_length_sections = opt[0][0];



  //TODO: Make this into a template function to avoid replication

	if (File::has_extension<Fibre::Tractlet>(match_location)) {

	  Fibre::Tractlet::Set reference(reference_location, degree);

	  if (!reference.size())
	    throw Exception ("Reference set, loaded from '" + reference_location + "' does not contain any fibres.");

	  if (!reference.has_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP)) {

      std::cout << std::endl << "Warning! Reference set does not have bundle_index property, will use strand order." << std::endl;

	    reference.add_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP,"-1");
	    for (size_t tract_i = 0; tract_i < reference.size(); ++tract_i)
	      reference.set_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP,str(tract_i), tract_i);

	  }

	  std::vector<size_t> bundle_indices;

	  for (size_t ref_i = 0; ref_i < reference.size(); ++ref_i)
	    bundle_indices.push_back(to<size_t>(reference.get_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP, ref_i)));

	  size_t new_bundle_index = *std::max_element(bundle_indices.begin(), bundle_indices.end());

	  //Check bundle indices are unique in reference set.
	  size_t num_bundles = bundle_indices.size();
	  std::unique(bundle_indices.begin(), bundle_indices.end());
	  if (bundle_indices.size() != num_bundles)
	    throw Exception ("Bundle indices were not unique in reference set '" + reference_location + "'.");



    Fibre::Tractlet::Set match(match_location);

	  std::vector<size_t> reference_match(match.size());


	  if (reference.size() > match.size()) {

	    MR::ProgressBar progress_bar("Getting closest match between match template and reference set...");

	    match.distance(reference, reference_match);

	  } else {

	    MR::ProgressBar progress_bar("Getting closest match between match template and reference set...");

	    std::vector<size_t> inv_reference_match;

	    reference.distance(match, inv_reference_match);

	    size_t extra_index = reference.size();

	    for (size_t tract_i = 0; tract_i < match.size(); ++tract_i) {

	      std::vector<size_t>::iterator index = std::find(inv_reference_match.begin(), inv_reference_match.end(), tract_i);

	      if (index != inv_reference_match.end())
	        reference_match[tract_i] = index - inv_reference_match.begin();
	      else {
	        reference_match[tract_i] = extra_index++;
	        bundle_indices.push_back(new_bundle_index++);
	      }

	    }


	  }

	  if (File::has_extension<Fibre::Tractlet>(input_location)) {

	    MR::ProgressBar progress_bar("Mapping bundle indices to output set...");

	    Fibre::Tractlet::Set input;

	    if (input_location != match_location) {
	      input.load(input_location);
	      if (input.size() != match.size())
	        throw Exception ("Size of match set (" + str(match.size()) + ") and input set (" + str(input.size()) + ") do not match.");
	    } else
	      input = match;

      if (!match.has_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP))
        input.add_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP, "-1");

      for (size_t tract_i = 0; tract_i < match.size(); ++tract_i) {

        input.set_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP, str(bundle_indices[reference_match[tract_i]]), tract_i);

        ++progress_bar;
      }

      input.save(output_location);


	  } else if (File::has_extension<Fibre::Tractlet::Set>(input_location)) {

      MR::ProgressBar progress_bar("Mapping bundle indices to output set...");

      if (input_location == output_location)
        throw Exception ("Output location ('" + output_location + "') cannot be the same as input location ('" + input_location + "') for tractlet sets.");

	    Fibre::Tractlet::Set::Reader reader(input_location);

	    std::vector<std::string> elem_header = reader.extend_elem_prop_keys();

	    if (std::find(elem_header.begin(), elem_header.end(), Fibre::Track::BUNDLE_INDEX_EPROP) == elem_header.end())
	      elem_header.push_back(Fibre::Track::BUNDLE_INDEX_EPROP);

	    Fibre::Tractlet::Set::Writer writer(output_location, reader, reader.extend_prop_keys(), elem_header, reader.get_extend_props());

	    Fibre::Tractlet::Set tractlets;

	    size_t count = 0;
	    while (reader.next(tractlets)) {

	      if (tractlets.size() != match.size())
	        throw Exception ("Size of tractlet set " + str(count) + " (" + str(tractlets.size()) + ") does not match that of 'match' set (" + str(match.size()) + ")");

	      tractlets.add_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP, "-1");

	      for (size_t tract_i = 0; tract_i < tractlets.size(); ++tract_i)
	        tractlets.set_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP, str(bundle_indices[reference_match[tract_i]]), tract_i);

	      writer.append(tractlets);
	      ++count;
        ++progress_bar;

	    }


	  } else
	    throw Exception ("Usupported output file type ('" + output_location + "') for match of type tractlet.");



	} else if (File::has_extension<Fibre::Strand>(match_location)) {

    Fibre::Strand::Set reference(reference_location, degree);

    if (!reference.size())
      throw Exception ("Reference set, loaded from '" + reference_location + "' does not contain any fibres.");

    if (!reference.has_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP)) {

      std::cout << std::endl << "Warning! Reference set does not have bundle_index property, will use strand order." << std::endl;

      reference.add_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP,"-1");
      for (size_t tract_i = 0; tract_i < reference.size(); ++tract_i)
        reference.set_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP,str(tract_i), tract_i);

    }
    std::vector<size_t> bundle_indices;

    for (size_t ref_i = 0; ref_i < reference.size(); ++ref_i)
      bundle_indices.push_back(to<size_t>(reference.get_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP, ref_i)));

    size_t new_bundle_index = *std::max_element(bundle_indices.begin(), bundle_indices.end());

    //Check bundle indices are unique in reference set.
    size_t num_bundles = bundle_indices.size();
    std::unique(bundle_indices.begin(), bundle_indices.end());
    if (bundle_indices.size() != num_bundles)
      throw Exception ("Bundle indices were not unique in reference set '" + reference_location + "'.");



    Fibre::Strand::Set match(match_location);

    std::vector<size_t> reference_match(match.size());


    if (reference.size() > match.size()) {

      MR::ProgressBar progress_bar("Getting closest match between match template and reference set...");

      match.distance(reference, reference_match);

    } else {

      MR::ProgressBar progress_bar("Getting closest match between match template and reference set...");

      std::vector<size_t> inv_reference_match;

      reference.distance(match, inv_reference_match);

      size_t extra_index = reference.size();

      for (size_t tract_i = 0; tract_i < match.size(); ++tract_i) {

        std::vector<size_t>::iterator index = std::find(inv_reference_match.begin(), inv_reference_match.end(), tract_i);

        if (index != inv_reference_match.end())
          reference_match[tract_i] = index - inv_reference_match.begin();
        else {
          reference_match[tract_i] = extra_index++;
          bundle_indices.push_back(new_bundle_index++);
        }

      }


    }

    if (File::has_extension<Fibre::Strand>(input_location)) {

      MR::ProgressBar progress_bar("Mapping bundle indices to output set...");

      Fibre::Strand::Set input;

      if (input_location != match_location) {
        input.load(input_location);
        if (input.size() != match.size())
          throw Exception ("Size of match set (" + str(match.size()) + ") and input set (" + str(input.size()) + ") do not match.");
      } else
        input = match;

      if (!match.has_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP))
        input.add_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP, "-1");

      for (size_t tract_i = 0; tract_i < match.size(); ++tract_i) {

        input.set_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP, str(bundle_indices[reference_match[tract_i]]), tract_i);

        ++progress_bar;
      }

      input.save(output_location);


    } else if (File::has_extension<Fibre::Strand::Set>(input_location)) {

      MR::ProgressBar progress_bar("Mapping bundle indices to output set...");

      if (input_location == output_location)
        throw Exception ("Output location ('" + output_location + "') cannot be the same as input location ('" + input_location + "') for tractlet sets.");

      Fibre::Strand::Set::Reader reader(input_location);

      std::vector<std::string> elem_header = reader.extend_elem_prop_keys();

      if (std::find(elem_header.begin(), elem_header.end(), Fibre::Track::BUNDLE_INDEX_EPROP) == elem_header.end())
        elem_header.push_back(Fibre::Track::BUNDLE_INDEX_EPROP);

      Fibre::Strand::Set::Writer writer(output_location, reader, reader.extend_prop_keys(), elem_header, reader.get_extend_props());

      Fibre::Strand::Set tractlets;

      size_t count = 0;
      while (reader.next(tractlets)) {

        if (tractlets.size() != match.size())
          throw Exception ("Size of tractlet set " + str(count) + " (" + str(tractlets.size()) + ") does not match that of 'match' set (" + str(match.size()) + ")");

        tractlets.add_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP, "-1");

        for (size_t tract_i = 0; tract_i < tractlets.size(); ++tract_i)
          tractlets.set_extend_elem_prop(Fibre::Track::BUNDLE_INDEX_EPROP, str(bundle_indices[reference_match[tract_i]]), tract_i);

        writer.append(tractlets);
        ++count;
        ++progress_bar;

      }


    } else
      throw Exception ("Usupported output file type ('" + output_location + "') for match of type tractlet.");



  } else
	  throw Exception ("Unsupported extension of match file '" + match_location+ "'.");


}
