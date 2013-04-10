/*
 Copyright 2008 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close, 18/02/2011.

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

#ifndef __bts_version_h__
#define __bts_version_h__

namespace BTS {
    
    const size_t BTS_MAJOR_VERSION = 0;
    const size_t BTS_MINOR_VERSION = 0;
    const size_t BTS_MICRO_VERSION = 3;
    const size_t BTS_NANO_VERSION = 24;
    
    inline std::string version_number_string() {
        return str(BTS_MAJOR_VERSION) + "." + str(BTS_MINOR_VERSION) + "." + str(BTS_MICRO_VERSION)
               + "." + str(BTS_NANO_VERSION);
    }
    
    inline bool before_version(const std::string& version_number_string, size_t major_version,
                               size_t minor_version = 0, size_t micro_version = 0,
                               size_t nano_version = 0) {
        
        if (!version_number_string.size())
            return true;    //If version number is not present assume that it is from an early version (or MRtrix).
            
        bool is_before = false;
        
        size_t major_dot_index = version_number_string.find(".");
        size_t read_major = to<size_t>(version_number_string.substr(0, major_dot_index));
        
        if (read_major > major_version)
            is_before = false;
        
        else if (read_major < major_version)
            is_before = true;
        
        else {
            
            size_t minor_dot_index = version_number_string.find(".", major_dot_index + 1);
            size_t read_minor = to<size_t>(
                    version_number_string.substr(major_dot_index + 1, minor_dot_index));
            
            if (read_minor > minor_version)
                is_before = false;
            else if (read_minor < minor_version)
                is_before = true;
            else {
                
                size_t micro_dot_index = version_number_string.find(".", minor_dot_index + 1);
                size_t read_micro = to<size_t>(
                        version_number_string.substr(minor_dot_index + 1, micro_dot_index));
                
                if (read_micro > micro_version)
                    is_before = false;
                else if (read_micro < micro_version)
                    is_before = true;
                else {
                    
                    size_t nano_dot_index = version_number_string.find(".", micro_dot_index + 1);
                    size_t read_nano = to<size_t>(
                            version_number_string.substr(micro_dot_index + 1, nano_dot_index));
                    
                    if (read_nano < nano_version)
                        is_before = true;
                    else
                        is_before = false;
                    
                }
                
            }
            
        }
        
        return is_before;
        
    }

}

#endif /* __bts_version_h__ */
