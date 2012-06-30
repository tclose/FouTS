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


extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
}


#include "bts/cmd.h"
#include "progressbar.h"

#include "bts/common.h"

#include "bts/fibre/base/writer.h"


#include "bts/coord.h"
#include "bts/fibre/tractlet/section.h"


#include "bts/mcmc/state.h"
//#include "bts/mcmc/state/writer.h"

#include "bts/fibre/properties/extended.h"
#include "bts/fibre/strand.h"


#include "bts/inline_functions.h"

using namespace BTS;

SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);

DESCRIPTION = {
  "Generates random triples and saves them to file.",
  "",
  NULL
};

ARGUMENTS = {
  Argument ("output", "Random triples").type_file(),
  Argument()
};


OPTIONS = {

  Option ("number", "The number of triples to generate")
   + Argument ("number", "").type_integer (1, 100, 1e6),

  Option ("max", "Maximum value of the triple elements.")
   + Argument ("max", "").type_float (-1e9, 1.0, 1e9),

  Option ("min", "Minimum value of the triple elements.")
   + Argument ("min", "").type_float (-1e9, 1.0, 1e9),
  
  Option ("state_size", "The size MCMC states to generate")
   + Argument ("state_size", "").type_integer (1, 2, 1e6),

  Option ("seed", "Seed for the random generation")
   + Argument ("seed", ""),


 
Option() };



EXECUTE {



  std::string output_location = argument[0];

  size_t number = 100;
  double max_value = 1.0;
  double min_value = 0.0;
  size_t state_size = 2;
  size_t seed = time(NULL);
  
  

  
  Options opt = get_options("number");
  if (opt.size())
    number = opt[0][0];

  opt = get_options("max_value");
  if (opt.size())
    max_value = opt[0][0];

  opt = get_options("min_value");
  if (opt.size())
    min_value = opt[0][0];
    
  opt = get_options("seed");
  if (opt.size()) {
    std::string seed_str = opt[0][0];
    seed = to<size_t>(seed_str);
  }
  opt = get_options("state_size");
  if (opt.size())
    state_size = opt[0][0];


  gsl_rng* rand_gen = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rand_gen, seed);




  if (File::has_extension< Coord >(output_location)) {

    MR::ProgressBar progress_bar ("Randomly generating " + str(number) + " Coord's...", number);

    Coord::Writer writer(output_location);
    
    for (size_t trip_i = 0; trip_i < number; ++trip_i) {
    
      Coord triple;
      
      for (size_t dim_i = 0; dim_i < 3; ++dim_i) 
        triple[dim_i] = (double)gsl_ran_flat(rand_gen, min_value, max_value);
    
      writer.append(triple);
      progress_bar++;
    
    }
    
    //MR::ProgressBar::done();
    
  } else if (File::has_extension<Fibre::Strand::Section>(output_location)) {

#ifdef CODE_TESTING
      MR::ProgressBar progress_bar ("Randomly generating " + str(number) + " Fibre::Strand::Section's...", number);

      std::vector<std::string> header;
      header.push_back("intensity");

      Fibre::Base::Writer<Fibre::Strand::Section> writer(output_location);

      for (size_t section_i = 0; section_i < number; ++section_i) {

        Fibre::Strand::Section section;

        for (size_t dim_i = 0; dim_i < 3; ++dim_i) {
          section.position()[dim_i] = (double)gsl_ran_flat(rand_gen, min_value, max_value);
          section.tangent()[dim_i] = (double)gsl_ran_flat(rand_gen, min_value, max_value);
        }

        section.intensity() = (double)gsl_ran_flat(rand_gen, min_value, max_value);

        writer.append(section);
        progress_bar++;

      }

      //MR::ProgressBar::done();
#else
      throw Exception ("Macro CODE_TESTING needs to be defined to generate random Fibre::Strand::Sections");
#endif /* CODE_TESTING */
    
  } else if (File::has_extension<Fibre::Tractlet::Section>(output_location)) {

#ifdef CODE_TESTING
      MR::ProgressBar progress_bar ("Randomly generating " + str(number) + " Fibre::Tractlet::Section's...", number);

      std::vector<std::string> header;
      header.push_back("intensity");

      Fibre::Base::Writer<Fibre::Tractlet::Section> writer(output_location);

      for (size_t section_i = 0; section_i < number; ++section_i) {

        Fibre::Tractlet::Section section;


        for (size_t dim_i = 0; dim_i < 3; ++dim_i) {
          section.position()[dim_i] = (double)gsl_ran_flat(rand_gen, min_value, max_value);
          section.tangent()[dim_i] = (double)gsl_ran_flat(rand_gen, min_value, max_value);
          section.width1()[dim_i] = (double)gsl_ran_flat(rand_gen, min_value, max_value);
          section.width2()[dim_i] = (double)gsl_ran_flat(rand_gen, min_value, max_value);
        }

        section.intensity() = (double)gsl_ran_flat(rand_gen, min_value, max_value);

        writer.append(section);
        progress_bar++;
      
      }
      
      //MR::ProgressBar::done();
#else
      throw Exception ("Macro CODE_TESTING needs to be defined to generate random Fibre::Tractlet::Sections");
#endif /* CODE_TESTING */

  } else if (File::has_extension<MCMC::State>(output_location)) {

      MR::ProgressBar progress_bar ("Randomly generating " + str(number) + " MCMC::State's of size " + str(state_size) + "...", number);

      MCMC::State::Writer writer(output_location);

      for (size_t state_i = 0; state_i < number; ++state_i) {

        MCMC::State state(state_size);

        for (size_t elem_i = 0; elem_i < state_size; elem_i++)
          state[elem_i] = (double)gsl_ran_flat(rand_gen, min_value, max_value);

        writer.append(state);
        progress_bar++;

      }

      //MR::ProgressBar::done();

    
  } else 
    throw Exception ("Unrecognised file extension '" + File::extension(output_location) + "', can be '" + Coord::FILE_EXTENSION + "', '" + Fibre::Strand::Section::FILE_EXTENSION + "', '" + Fibre::Tractlet::Section::FILE_EXTENSION + "' or '" + MCMC::State::FILE_EXTENSION + "'.");
  


}
