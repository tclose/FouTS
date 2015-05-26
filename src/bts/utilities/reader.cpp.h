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

#include "bts/utilities/reader.h"
#include "bts/file.h"

namespace FTS {
    
    namespace Utilities {
        
        template<typename T> void Reader<T>::open(const std::string& location) {
            
            if (File::is_file(location)) {
                
                fin.open((location).c_str());
                
                if (!fin.good())
                    throw Exception("Error opening file '" + location + "'.");
                
            } else
                throw Exception("No file at location '" + location + "'");
            
        }
        
        template<typename T> bool Reader<T>::next(T& state, Properties& properties_row) {
            
            if (!fin.good())
                throw Exception("Corrupted utilities file.");
            
            std::string line;
            std::getline(fin, line);
            
            if (line == "%END")
                return false;
            
            size_t count = 0;
            size_t pos = 0;
            
            while (pos != std::string::npos) {
                pos = line.find(' ', pos + 1);
                count++;
            }
            
            state.resize(count);
            
            std::istringstream line_stream(line);
            
            std::string value;
            
            for (size_t elem_i = 0; elem_i < state.size(); elem_i++) {
                std::getline(line_stream, value, ' ');
                state[elem_i] = to<double>(value);
            }
            
            return true;
            
        }
    
    }

}
