/*
 Copyright 2010 Brain Research Institute, Melbourne, Australia

 Written by Thomas G Close on Jun 3, 2010.

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

#ifndef __bts_cmd_h__
#define __bts_cmd_h__

#include "app.h"
#include "args.h"
#include "progressbar.h"

namespace FTS {
    
    class Argument: public MR::Argument {
            
        public:
            Argument(const char* name = NULL, const char* description = NULL)
                    : MR::Argument(name, description) {
            }
            
    };
    
    class Option: public MR::Option {
            
        public:
            
            typedef enum _ByDefault {
                Active, Inactive
            } ByDefault;

        protected:
            
            ByDefault by_default;

        public:
            Option()
                    : by_default(Inactive) {
            }
            
            Option(const char* name, const char* description)
                    : MR::Option(name, description), by_default(Inactive) {
            }
            
            Option& active_by_default(bool is_active = true) {
                by_default = is_active ? Active : Inactive;
                return *this;
            }
            
    };
    
    class App: public MR::App {
            
        public:
            App(int argc, char** argv, const char** cmd_desc, const MR::Argument* cmd_args,
                const MR::Option* cmd_opts, const size_t* cmd_version, const char* cmd_author,
                const char* cmd_copyright)
                    : MR::App(argc, argv, cmd_desc, cmd_args, cmd_opts, cmd_version, cmd_author,
                            cmd_copyright) {
            }
            
//      void   print_help () const;
            
            virtual ~App();
    };

//  MR::ActiveByDefault NA       = MR::NA;
//  MR::ActiveByDefault Active   = MR::Active;
//  MR::ActiveByDefault Inactive = MR::Inactive;

}

#endif /* __bts_cmd_h__ */
