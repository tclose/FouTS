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

#include <cstdlib>
#include <stdexcept>

#include "bts/common.h"

namespace BTS {
    
    void breakpoint() {
        size_t i = 0;
        i++;
    }
    
    double rand_double() {
        
        size_t rand_size_t = rand();
        
        return (double) rand_size_t / (double) RAND_MAX;
        
    }
    
    std::vector<double> parse_double_vector(const std::string& s) {
        
        std::vector<double> v;
        
        size_t start_number = 0;
        size_t end_number;
        
        int i = 0;
        
        start_number = s.find_first_of("0123456789.-");
        
        while (start_number != std::string::npos) {
            
            end_number = s.find_first_not_of("0123456789.-", start_number);
            
            std::string num_string = s.substr(start_number, end_number - start_number);
            
            if (num_string.find_first_of('-', 1) != std::string::npos || num_string.find('.',
                                                                                 num_string.find(
                                                                                         '.')
                                                                                 + 1)
                                                                         != std::string::npos)    // If there is a '-' character anywhere other than the first character, or if there is more than one decimal point then throw an error.
                throw Exception(
                        "Value, " + num_string + " at position: " + str(i + 1) + ", is invalid");
            
            v.push_back(atof(num_string.c_str()));
            
            start_number = s.find_first_of("0123456789.-", end_number);
            i++;
            
        }
        
        return v;
        
    }
    
    std::vector<int> parse_int_vector(const std::string& s) {
        
        std::vector<int> v;
        
        size_t start_number = 0;
        size_t end_number;
        
        int i = 0;
        
        start_number = s.find_first_of("0123456789-");
        
        while (start_number != std::string::npos) {
            
            end_number = s.find_first_not_of("0123456789-", start_number);
            
            std::string num_string = s.substr(start_number, end_number - start_number);
            
            if (num_string.find_first_of('-', 1) != std::string::npos)    // If there is a '-' character anywhere other than the first character then throw an error.
                throw Exception(
                        "Value, " + num_string + " at position: " + str(i + 1) + ", is invalid");
            
            v.push_back(atoi(num_string.c_str()));
            
            start_number = s.find_first_of("0123456789-", end_number);
            i++;
            
        }
        
        return v;
        
    }
    
    //TODO: Make this a template function and improve the parse to include ranges.
    std::vector<size_t> parse_size_t_vector(const std::string& s) {
        
        std::vector<size_t> v;
        
        size_t start_number = 0;
        size_t end_number;
        
        size_t i = 0;
        
        start_number = s.find_first_of("0123456789-");
        
        while (start_number != std::string::npos) {
            
            end_number = s.find_first_not_of("0123456789-", start_number);
            
            std::string num_string = s.substr(start_number, end_number - start_number);
            
            if (num_string.find_first_of('-', 1) != std::string::npos)    // If there is a '-' character anywhere other than the first character then throw an error.
                throw Exception(
                        "Value, " + num_string + " at position: " + str(i + 1) + ", is invalid");
            
            v.push_back(atoi(num_string.c_str()));
            
            start_number = s.find_first_of("0123456789-", end_number);
            i++;
            
        }
        
        return v;
        
    }
    
    int num_dec_places(int number) {
        int num_places, ten_power;
        
        num_places = 1;
        ten_power = 10;
        while ((number / ten_power) != 0) {
            num_places++;
            ten_power *= 10;
        }
        
        return num_places;
    }
    
    int num_precision_places(double number) {
        int num_places = 0;
        
        double remainder = number - floor(number);
        
        while (remainder != 0.0) {
            num_places++;
            remainder *= 10;
            remainder = remainder - floor(remainder);
        }
        
        return num_places;
    }

}

