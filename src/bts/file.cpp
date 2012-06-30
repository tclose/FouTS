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

#include <algorithm>

#include "bts/common.h"
#include "bts/file.h"


namespace BTS {

  namespace File {
 
    std::vector<std::string>              list_filenames(std::string dir_location, std::string prefix, std::string suffix) {
      
      if (!File::is_dir(dir_location))
        throw MR::Exception("Input location provided " + dir_location + " is not a directory");
        
      File::Dir dir(dir_location);
        
      std::vector<std::string> filenames;
        
      std::string filename;
        
        
      while (!(filename = dir.read_name()).empty()) {    
           
        if (filename[0] != '.' && filename.size() - prefix.size() - suffix.size() >= 0) {      //If the file is not hidden or has a name shorter than the prefix and suffix combined. 
                     
          if (!filename.substr(0,prefix.size()).compare(prefix) && !filename.substr(filename.size() - suffix.size(), filename.size()).compare(suffix))  //If the file matches the prefix and suffix
            filenames.push_back(filename);
        }
      }
        
      sort(filenames.begin(), filenames.end());
    
      return filenames;
    
    }
    
    
    std::pair<std::string, std::string>   split_at_extension (std::string filename) {
    
      std::pair<std::string, std::string> basename_and_ext;
    
      size_t extension_i = filename.find_last_of(".");

      if (extension_i == std::string::npos)
        basename_and_ext.second = "";
      else
        basename_and_ext.second = filename.substr(extension_i);
        
      basename_and_ext.first = filename.substr(0, extension_i);
    
      return basename_and_ext;
    }
    
    
    void                                  clear_path(const std::string& location) {

      if (File::exists(location)) {

        if (OVERWRITE_FILES)
          remove(location.c_str());
        else
          throw Exception ("File '" + location + "' exists, if you want to automatically overwrite it set OVERWRITE_FILES to true.");

      }

    }


    std::string                           extension(std::string location) {

      std::string file_ext;

      size_t extension_i = location.find_last_of(".");

      if (extension_i == std::string::npos || extension_i == location.size()-1)
        file_ext = "";
      else
        file_ext = location.substr(extension_i+1);

      return file_ext;
    }


    std::string                           strip_extension(std::string location) {

      std::string stripped;

      size_t extension_i = location.find_last_of(".");

      if (extension_i == std::string::npos)
        stripped = location;
      else
        stripped = location.substr(0, extension_i);

      return stripped;
    }


  }
  
}
