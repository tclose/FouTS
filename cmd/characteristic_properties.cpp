/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Created by Tom Close on 13/03/09.

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

#include "bts/common.h"

#include "bts/cmd.h"
#include "bts/file.h"
#include "progressbar.h"

#include "bts/fibre/strand/set.h"
#include "bts/fibre/tractlet/set.h"

#include "bts/inline_functions.h"


using namespace BTS;

template <typename T> void characteristic_properties(const std::string& input_location, const std::string& output_location);

SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);

DESCRIPTION = {
  "Normalises the density of the given Fourier tracts.",
  "Normalises the density of the given Fourier tracts.",
  NULL
};

ARGUMENTS = {
  Argument ("input", "The Fourier tracts to be normalised.").type_file (),
  Argument ("output", "Output location of the normalised tracts (if omitted the origin location will be used)").optional().type_file(),
  Argument()
};


OPTIONS = {

Option() };


EXECUTE {


  std::string input_location = argument[0];
  std::string output_location;

  if (argument.size() > 1)
    output_location = argument[1].c_str();
  else
    output_location = input_location;

  MR::ProgressBar progress_bar ("Calculating characteristic properties...");

  if (File::has_extension<Fibre::Strand>(input_location))
    characteristic_properties<Fibre::Strand::Set>(input_location, output_location);
  else if (File::has_extension<Fibre::Tractlet>(input_location))
    characteristic_properties<Fibre::Tractlet::Set>(input_location, output_location);
  else
    throw Exception("Unrecognised extension to file '" + input_location + "' (can be " + Fibre::Strand::FILE_EXTENSION + "' or '" + Fibre::Tractlet::FILE_EXTENSION + "')");
}


template <typename T> void characteristic_properties(const std::string& input_location, const std::string& output_location) {

  T fibres(input_location);
  std::vector<std::string> keys;
  fibres.append_characteristic_keys(keys);
  for (std::vector<std::string>::iterator key_it = keys.begin(); key_it != keys.end(); ++key_it)
    fibres.add_extend_elem_prop(*key_it, "NAN");
  fibres.set_characteristics();
  fibres.save(output_location);

}





