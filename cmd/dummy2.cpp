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

//#include "bts/fibre/strand/set.h"

//#include "bts/fibre/strand/set.h"
//#include "bts/fibre/track/set.h"
//#include "bts/fibre/tractlet/set.h"
//#include "bts/fibre/tractlet/section.h"
//#include "bts/fibre/base/set_reader.h"
//#include "bts/fibre/base/set_writer.h"

using namespace BTS;

SET_VERSION_DEFAULT;
SET_AUTHOR ("Thomas G. Close");
SET_COPYRIGHT (NULL);

DESCRIPTION = {
  "Print strands to screen",
  "",
  NULL
};

ARGUMENTS = {
  Argument()
};


OPTIONS = {

Option() };

EXECUTE {

  printf("Hello World\n");
  
}

