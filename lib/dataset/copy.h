/*
    Copyright 2009 Brain Research Institute, Melbourne, Australia

    Written by J-Donald Tournier, 20/10/09.

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

#ifndef __dataset_copy_h__
#define __dataset_copy_h__

#include "debug.h"
#include "dataset/loop.h"

namespace MR
{
  namespace DataSet
  {

    template <class Set, class Set2>
    void copy (Set& destination, Set2& source, size_t from_axis = 0, size_t to_axis = std::numeric_limits<size_t>::max())
    {
      LoopInOrder loop (source, from_axis, to_axis);
      for (loop.start (destination, source); loop.ok(); loop.next (destination, source))
        destination.value() = source.value();
    }



    template <class Set, class Set2>
    void copy_with_progress (Set& destination, Set2& source, size_t from_axis = 0, size_t to_axis = std::numeric_limits<size_t>::max())
    {
      LoopInOrder loop (source, "copying from \"" + shorten (source.name()) + "\" to \"" + shorten (destination.name()) + "\"...", from_axis, to_axis);
      for (loop.start (destination, source); loop.ok(); loop.next (destination, source))
        destination.value() = source.value();
    }


    template <class Set, class Set2>
    void copy_with_progress_message (const std::string& message, Set& destination, Set2& source, size_t from_axis = 0, size_t to_axis = std::numeric_limits<size_t>::max())
    {
      LoopInOrder loop (source, message, from_axis, to_axis);
      for (loop.start (destination, source); loop.ok(); loop.next (destination, source))
        destination.value() = source.value();
    }

  }
}

#endif
