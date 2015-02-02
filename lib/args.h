/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by J-Donald Tournier, 27/06/08.

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

#ifndef __args_h__
#define __args_h__

#include <vector>
#include <limits>
#include "ptr.h"

#ifdef None
# undef None
#endif

namespace MR
{

  /*! \defgroup CmdParse Command-Line Parsing
   * \brief Classes and functions to parse command-line arguments and options.
   *
   * For a detailed description of the command-line parsing interface, see the
   * \ref command_line_parsing page.
   * */

  //! \cond skip
  typedef enum {
    Text,
    Integer,
    Float,
    ArgFile,
    Choice,
    ImageIn,
    ImageOut,
    IntSeq,
    FloatSeq
  } ArgType;

  const char* argtype_description (ArgType type);

  typedef int ArgFlags;
  const ArgFlags None = 0;
  const ArgFlags Optional = 0x1;
  const ArgFlags AllowMultiple = 0x2;
  //! \endcond



  //! \addtogroup CmdParse
  // @{

  //! A class to specify a command-line argument
  /*! Command-line arguments that are accepted by a particular command are
   * specified as an array of Arguments objects, terminated with an empty
   * Argument (constructed using default parameters). Please refer to \ref
   * command_line_parsing for more information.
   *
   * The list of arguments is provided using the ARGUMENTS macro, like this:
   * \code
   * ARGUMENTS = {
   *
   *   Argument ("input", "the input image")
   *     .type_image_in(),
   *
   *   Argument ("parameter",
   *        "the parameter to use during processing. Allowed values are "
   *        "between 0 and 10 (default = 1).")
   *     .type_float (0.0, 1.0, 10.0),
   *
   *   Argument ("output", "the output image")
   *     .type_image_out(),
   *
   *   Argument () // don't forget to terminate with the default Argument
   * };
   * \endcode
   * The example above specifies that the application expects exactly 3
   * arguments, with the first one being an existing image to be used as input,
   * the second one being a floating-point value, and the last one being an
   * image to be created and used as output.
   *
   * There are a number of types that the argument can be specified as. The
   * argument can also be specified as optional (see optional() function), or
   * as multiple (see allow_multiple() function). Note that in this case only
   * one such argument can be optional and/or multiple, since more than one
   * such argument would lead to ambiguities when parsing the command-line.  */
  class Argument
  {
    public:
      //! constructor
      /*! this is used to construct a command-line argument object, with a name
       * and description. If default arguments are used, the object corresponds
       * to the end-of-list specifier, as detailed in \ref command_line_parsing. */
      Argument (const char* name = NULL, const char* description = NULL) :
        id (name), desc (description), type (Text), flags (None) {
        defaults.text = NULL;
      }

      //! the argument name
      const char* id;
      //! the argument description
      const char* desc;
      //! the argument type
      ArgType  type;
      //! the argument flags (AllowMultiple and/or Optional)
      ArgFlags flags;

      //! a structure to store the various parameters of the Argument
      union {
        const char* text;
        struct {
          const char** list;
          int def;
        } choices;
        struct {
          int def, min, max;
        } i;
        struct {
          float def, min, max;
        } f;
      } defaults;


      operator bool () const {
        return (id);
      }

      //! specifies that the argument is optional
      /*! For example:
       * \code
       * ARGUMENTS = {
       *
       *   Argument ("input", "the input image")
       *     .type_image_in()
       *     .optional()
       *     .allow_multiple(),
       *
       *   Argument ()
       * };
       * \endcode
       * \note Only one argument can be specified as optional and/or multiple.
       */
      Argument& optional () {
        flags |= Optional;
        return (*this);
      }

      //! specifies that multiple such arguments can be specified
      /*! See optional() for details. */
      Argument& allow_multiple () {
        flags |= AllowMultiple;
        return (*this);
      }

      //! specifies that the argument should be a text string
      /*! If desired, a default string can be specified using the \a
       * default_text argument. */
      Argument& type_text (const char* default_text = NULL) {
        type = Text;
        defaults.text = default_text;
        return (*this);
      }

      //! specifies that the argument should be an input image
      Argument& type_image_in () {
        type = ImageIn;
        defaults.text = NULL;
        return (*this);
      }

      //! specifies that the argument should be an output image
      Argument& type_image_out () {
        type = ImageOut;
        defaults.text = NULL;
        return (*this);
      }

      //! specifies that the argument should be an integer
      /*! if desired, a default value can be specified, along with a range of
       * allowed values. */
      Argument& type_integer (int min = std::numeric_limits<int>::min(), int def = 0, int max = std::numeric_limits<int>::max()) {
        type = Integer;
        defaults.i.min = min;
        defaults.i.def = def;
        defaults.i.max = max;
        return (*this);
      }

      //! specifies that the argument should be a floating-point value
      /*! if desired, a default value can be specified, along with a range of
       * allowed values. */
      Argument& type_float (float min = -INFINITY, float def = 0.0, float max = INFINITY) {
        type = Float;
        defaults.f.min = min;
        defaults.f.def = def;
        defaults.f.max = max;
        return (*this);
      }

      //! specifies that the argument should be selected from a predefined list
      /*! The list of allowed values must be specified as a NULL-terminated
       * list of C strings.  If desired, a default value can be specified,
       * in the form of an index into the list. Here is an example usage:
       * \code
       * const char* mode_list [] = { "standard", "pedantic", "approx", NULL };
       *
       * ARGUMENTS = {
       *
       *   Argument ("mode", "the mode of operation")
       *     .type_choice (mode_list, 0),
       *
       *   Argument ()
       * };
       * \endcode
       * \note Each string in the list must be supplied in \b lowercase. */
      Argument& type_choice (const char** choices, int default_index = -1) {
        type = Choice;
        defaults.choices.list = choices;
        defaults.choices.def = default_index;
        return (*this);
      }

      //! specifies that the argument should be a file
      Argument& type_file () {
        type = ArgFile;
        defaults.text = NULL;
        return (*this);
      }

      //! specifies that the argument should be a sequence of comma-separated integer values
      Argument& type_sequence_int () {
        type = IntSeq;
        defaults.text = NULL;
        return (*this);
      }

      //! specifies that the argument should be a sequence of comma-separated floating-point values.
      Argument& type_sequence_float () {
        type = FloatSeq;
        defaults.text = NULL;
        return (*this);
      }


      void print () const;
      void print_usage () const;
  };




  //! A class to specify a command-line option
  /*! Command-line options that are accepted by a particular command are
   * specified as an array of Option objects, terminated with an empty
   * Option (constructed using default parameters). Please refer to \ref
   * command_line_parsing for more information.
   *
   * The list of options is provided using the OPTIONS macro, like this:
   * \code
   * OPTIONS = {
   *
   *   Option ("exact", "do not use approximations when processing"),
   *
   *   Option ("mask",
   *        "only perform processing within the voxels contained in "
   *        "the binary image specified")
   *     + Argument ("image").type_image_in(),
   *
   *   Option ("regularisation", "set the regularisation term")
   *     + Argument ("value").type_float (0.0, 1.0, 100.0),
   *
   *   Option ("dump", "dump all intermediate values to file")
   *     + Argument ("file").type_file(),
   *
   *   Option () // don't forget to terminate with the default Option
   * };
   * \endcode
   * The example above specifies that the application accepts four options, in
   * addition to the standard ones (see \ref command_line_parsing for details).
   * The first option is a simple switch: specifying '-exact' on the
   * command line will cause the application to behave differently.
   * The next options all expect an additional argument, supplied using the
   * Argument class. Note that the description field of the Argument class is
   * unused in this case. Multiple additional arguments can be specified in
   * this way using the addition operator.
   *
   * Options can also be specified as required (see required() function), or
   * as multiple (see allow_multiple() function).
   */
  class Option
  {
    public:
      Option () : id (NULL), desc (NULL), flags (Optional) { }

      Option (const char* name, const char* description) :
        id (name), desc (description), flags (Optional) { }

      Option& operator+ (const Argument& arg) {
        args.push_back (arg);
        return (*this);
      }
      operator bool () const {
        return (id);
      }

      //! the option name
      const char* id;
      //! the option description
      const char* desc;
      //! option flags (AllowMultiple and/or Optional)
      ArgFlags flags;

      //! a vector of argument that must be supplied with the option
      std::vector<Argument> args;

      //! specifies that the option is required
      /*! An option specified as required must be supplied on the command line.
        * For example:
       * \code
       * OPTIONS = {
       *
       *   Option ("roi",
       *       "the region of interest over which to perform the processing. "
       *       "Mulitple such regions can be specified")
       *     .required()
       *     .allow_multiple()
       *     + Argument ("image").type_image_in(),
       *
       *   Argument ()
       * };
       * \endcode
       * \note Only one argument can be specified as optional and/or multiple.
       */
      Option& required () {
        flags &= ~Optional;
        return (*this);
      }

      //! specifies that multiple such options can be specified
      /*! See required() for details. */
      Option& allow_multiple () {
        flags |= AllowMultiple;
        return (*this);
      }

      void print () const;
      void print_usage () const;
  };

  // @}

}

#endif

