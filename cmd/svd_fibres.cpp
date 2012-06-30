/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, Jul 26, 2010.

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

#include "bts/file.h"

#include "bts/math/svd.h"



#include "bts/inline_functions.h"

using namespace BTS;


template <typename T> void        svd_fibres(const std::string& input_location, const std::string& output_location, size_t max_num_components);


SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);


DESCRIPTION = {
  "svd_fibres",
  "",
  NULL
};

ARGUMENTS = {
  Argument ("input location","").type_file (),
  Argument ("output location","").type_file (),
  Argument()
};


OPTIONS = {

  Option("max_num_components","The number of principle components to be returned..")
   + Argument("max_num_components","").type_integer(0,UINT_MAX,UINT_MAX),

Option() };


EXECUTE {


  std::string input_location = argument[0];
	std::string output_location = argument[1];

	size_t max_num_components = UINT_MAX;

  Options opt = get_options("max_num_components");
  if (opt.size())
    max_num_components = opt[0][0];

	if (File::extension(input_location) != File::extension(output_location))
	  throw Exception ("Input location ('" + File::extension(input_location) + "') does not match output location ('" + File::extension(output_location) + "').");

  MR::ProgressBar progress_bar ("Performing SVD on fibres ...");

  if (File::has_extension<Fibre::Strand::Set>(input_location)) {

    svd_fibres<Fibre::Strand::Set>(input_location, output_location, max_num_components);

  } else if (File::has_extension<Fibre::Tractlet::Set>(input_location)) {

    svd_fibres<Fibre::Tractlet::Set>(input_location, output_location, max_num_components);

  } else
	  throw Exception ("Unrecognised extension '" + File::extension(input_location) + "'.");


}

template <typename T> void       svd_fibres(const std::string& input_location, const std::string& output_location, size_t max_num_components) {

  typename T::Reader reader(input_location);

  T input_fibres, output_fibres;

  size_t num_cols = 0;

  std::vector< MR::Math::Vector<double> > input_vectors;

  while (reader.next(input_fibres)) {

    MR::Math::Vector<double>& input_vector  = input_fibres;

    if (!num_cols) //First one.
      num_cols = input_vector.size();
    else if (num_cols != input_vector.size())
      throw Exception ("Inconsistent lengths (" + str(num_cols) + " and " + str(input_vector.size()) + ") found in input location.");

    input_vectors.push_back(input_vector);

    //Save the input_fibres structure so that the usv matrix can be converted back into the right format.
    if (!output_fibres.size())
      output_fibres = input_fibres;

  }

  size_t num_rows = input_vectors.size();

  MR::Math::Matrix<double> M (num_rows, num_cols);

  for (size_t row_i = 0; row_i < num_rows; row_i++)
    M.row(row_i) = input_vectors[row_i];


  Math::USV usv = Math::svd(M);

  std::vector<std::string> ext_header;
  ext_header.push_back("singular_value");

  std::map<std::string,std::string> file_props;
  file_props["method"] = "svd_fibres";

  typename T::Writer writer(output_location, reader, ext_header, file_props);

  size_t num_components = min2(max_num_components, num_cols);

  for (size_t comp_i = 0; comp_i < num_components; comp_i++) {

    output_fibres = usv.v.column(comp_i);

    output_fibres.set_extend_prop("singular_value", str(usv.sigma[comp_i]));

    writer.append(output_fibres);

  }

}
















































