/*
 Copyright 2009 Brain Research Institute, Melbourne, Australia

 Created by Tom Close on 13/03/09.

 This file is part of Fourier Tract Sampling (FouTS).

 FouTS is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 FouTS is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with FTS.  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef __bts_file_h__
#define __bts_file_h__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

//#include "bts/common.h"

#include "file/path.h"

#include "exception.h"

#include "types.h"

namespace FTS {
    
    namespace File {
        
        const bool OVERWRITE_FILES = true;
        
        typedef MR::Path::Dir Dir;
        
        const std::string TXT_FILE_EXTENSTION = "txt";
        
        inline bool exists(const std::string& path) {
            return MR::Path::exists(path);
        }
        
        inline void remove(const std::string& location) {
            std::remove(location.c_str());
        }
        
        inline bool is_dir(const std::string path) {
            return MR::Path::is_dir(path);
        }
        
        inline bool is_file(const std::string path) {
            return MR::Path::is_file(path);
        }
        
        inline std::string join(const std::string& first, const std::string& second) {
            return MR::Path::join(first, second);
        }
        
        std::pair<std::string, std::string> split_at_extension(std::string filename);
        
        std::string extension(std::string location);
        
        std::string strip_extension(std::string location);
        
        void clear_path(const std::string& location);
        
        inline std::string basename(std::string location) {
            return MR::Path::basename(location);
        }
        
        inline std::string dirname(std::string location) {
            return MR::Path::dirname(location);
        }
        
        inline bool has_extension(const std::string& name, const std::string& suffix) {
            return MR::Path::has_suffix(name, suffix);
        }
        
        inline bool has_txt_extension(const std::string& name, const std::string& suffix) {
            return has_extension(name, TXT_FILE_EXTENSTION) && has_extension(strip_extension(name),
                           suffix);
        }
        
        inline bool has_or_txt_extension(const std::string& name, const std::string& suffix) {
            return has_txt_extension(name, suffix) || has_extension(name, suffix);
        }
        
        template<typename T> bool has_extension(const std::string& location) {
            return MR::Path::has_suffix(location, T::FILE_EXTENSION);
        }
        
        template<typename T> bool has_or_set_extension(const std::string& location) {
            return has_extension<T>(location) || has_extension<typename T::Set>(location);
        }
        
        template<typename T> bool has_txt_extension(const std::string& name) {
            return has_txt_extension(name, T::FILE_EXTENSION);
        }
        
        template<typename T> bool has_or_txt_extension(const std::string& name) {
            return has_or_txt_extension(name, T::FILE_EXTENSION);
        }
        
        template<typename T> bool has_txt_or_set_extension(const std::string& location) {
            return has_or_set_extension<T>(location) || has_or_txt_extension<typename T::Set>(
                           location);
        }
        
        std::vector<std::string> list_filenames(std::string dir_location, std::string prefix =
                std::string(""),
                                                std::string suffix = std::string(""));
        
        template<typename T> std::vector<T>& load_vector(std::string location,
                                                         std::vector<T>& vector) {
            
            vector.clear();
            
            std::ifstream fin(location.c_str());
            
            char line[256];
            T t;
            while (fin.getline(line, 256))
                vector.push_back(to(line, t));
            
            fin.close();
            
            return vector;
            
        }
        
        template<typename T> void save_vector(std::string location, const std::vector<T>& vector) {
            
            std::ofstream fout(location.c_str());
            
            for (size_t i = 0; i < vector.size(); i++)
                fout << vector[i] << std::endl;
            
            fout.close();
            
        }
        
        inline void mkdir(const std::string& location, bool make_parents = false) {
            
            std::string cmd = "mkdir ";
            
            if (make_parents)
                cmd += "-p ";
            
            cmd += location;
            
            int failed = system(cmd.c_str());
            
            if (failed)
                throw MR::Exception("Could not make directory '" + location + "'.");
            
        }
    
    }

}

#endif

