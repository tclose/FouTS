/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, 12/08/2010.

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

#include "bts/utilities/writer.h"
#include "bts/file.h"

namespace BTS {
    
    namespace Utilities {
        
        template<typename T> void Writer<T>::create(const std::string& location) {
            
            fout.open((location).c_str());
            
            if (!fout.good())
                throw Exception("Error opening file '" + location + "'.");
            
            fout << "%END\n";
            fout.flush();
            
        }
        
        template<typename T> void Writer<T>::append(const T& state) {
            
            if (state.size()) {
                
                int64_t start_pos = fout.tellp();
                start_pos -= sizeof("%END\n") - 1;
                fout.seekp(start_pos);
                
                for (size_t elem_i = 0; elem_i < state.size() - 1; ++elem_i)
                    fout << state[elem_i] << " ";
                
                fout << state[state.size() - 1] << "\n%END\n";
                fout.flush();
                
            } else
                throw Exception("No elements found in state.");
        }
    
    }

}
