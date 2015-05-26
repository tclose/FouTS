#include <gsl/gsl_version.h>

#include "app.h"
#include "debug.h"
#include "progressbar.h"
#include "file/path.h"
#include "file/config.h"

#include "bts/cmd.h"

namespace FTS {
//
//  void App::print_help () const
//  {
//    fprintf (stderr, "%s: part of the MRtrix package\n\n", App::name().c_str());
//    if (command_description[0]) {
//      print_formatted_paragraph ("", command_description[0], HELP_PURPOSE_INDENT);
//      fprintf (stderr, "\n");
//      for (const char** p = command_description+1; *p; p++) {
//        print_formatted_paragraph ("", *p, HELP_PURPOSE_INDENT);
//        fprintf (stderr, "\n");
//      }
//    }
//    else fprintf (stderr, "(no description available)\n\n");
//
//    fprintf (stderr, "SYNTAX: %s [ options ]", App::name().c_str());
//    for (const Argument* arg = command_arguments; *arg; ++arg) {
//
//      if (arg->flags & Optional)
//        fprintf (stderr, " [");
//
//      fprintf (stderr, " %s", arg->id);
//
//      if (arg->flags & AllowMultiple) {
//        if (! (arg->flags & Optional))
//          fprintf (stderr, " [ %s", arg->id);
//        fprintf (stderr, " ...");
//      }
//      if (arg->flags & (Optional | AllowMultiple))
//        fprintf (stderr, " ]");
//    }
//    fprintf (stderr, "\n\n");
//
//
//
//    for (const Argument* arg = command_arguments; *arg; ++arg)
//      print_formatted_paragraph (std::string ("- ") + arg->id + " ", arg->desc, HELP_ARG_INDENT);
//    fprintf (stderr, "\n");
//
//
//    fprintf (stderr, "OPTIONS:\n\n");
//    for (const Option* opt = command_options; *opt; ++opt) {
//      std::string text ("-");
//      text += opt->id;
//
//      for (std::vector<Argument>::const_iterator optarg = opt->args.begin();
//           optarg != opt->args.end(); ++optarg)
//        text += std::string (" ") + optarg->id;
//
//      print_formatted_paragraph (text + " ", opt->desc, HELP_OPTION_INDENT);
//      // TODO: add argument defaults like this:
//      //print_formatted_paragraph (text + " ", opt->desc + opt->arg_defaults(), HELP_OPTION_INDENT);
//
//      fprintf (stderr, "\n");
//    }
//
//    fprintf (stderr, "Standard options:\n");
//    for (const Option* opt = default_options; *opt; ++opt)
//      print_formatted_paragraph (std::string ("-") + opt->id, opt->desc, HELP_OPTION_INDENT);
//    fprintf (stderr, "\n");
//  }

}
