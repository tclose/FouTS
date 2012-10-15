/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close, Jan 27, 2011.

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

#include "exception.h"

#include "bts/common.h"
#include "bts/file.h"

#include "bts/fibre/base/object.h"
#include "bts/fibre/strand.h"
#include "bts/fibre/tractlet.h"
#include "bts/fibre/strand/set.h"
#include "bts/fibre/tractlet/set.h"
#include "bts/math/common.h"

namespace BTS {

	namespace Fibre {

	  namespace Base {

      const char*                Object::ALPHA_PROP = "acs_sqrt";

      const char*                Object::PROPS_LIST_END = 0;

      // Sorry a bit of a hack. Is passed to the operator[]() when constructing a strand element without the properties.
      std::vector<const char*>   Object::EMPTY_PROPS;


      void                       Object::reset(size_t size, size_t row_size, const std::vector<const char*>& properties, double fill_value) {

        if (!owner)
          throw Exception ("Cannot reset object as it is a view onto a larger object.");

        MR::Math::Vector<double>::resize(size * row_size + properties.size());
        sze = size;

        *props = properties;

        set(fill_value);

      }


      void                       Object::reset(const std::vector<const char*>& properties, double fill_value) {

        if (!owner)
          throw Exception ("Cannot reset object as it is a view onto a larger object.");

        MR::Math::Vector<double>::resize(properties.size(), fill_value);
        sze = 0;

        *props = properties;

        set(fill_value);

      }


      bool                       Object::props_match(const std::vector<const char*>* properties) const {

        bool match;

        if (props == properties) // Point to the same properties
          match = true;
        else if (!props->size() && !properties->size()) // Both contain no properties
          match = true;
        else if (num_props() != properties->size()) // Are of different size
          match = false;
        else // Compare the property pointers (relies on the fact that properties are sorted by order of memory location when they are inserted).
          match = !memcmp(&(props->operator[](0)),&(properties->operator[](0)), num_props());

        return match;

      }



      Object&                    Object::operator=(const Object& b) {

        if (is_owner()) {

          MR::Math::Vector<double>::operator=(b);

          sze = b.sze;
          *props = *b.props;


        } else {

          // Might consider making these Exceptions instead of assertions, but this function probably should be as fast as possible.
          assert(vsize() == b.vsize()); // Sizes need to match when fibre object is not owner of data.
          assert(props_match(b)); // Properties need to match exactly when fibre object is not owner of data.

          const MR::Math::Vector<double>& bvec = b;

          for (size_t i = 0; i < vsize(); ++i)
            MR::Math::Vector<double>::operator[](i) = bvec[i];

        }

        return *this;

      }


      Object&                    Object::operator=(const MR::Math::Vector<double>& v) {

        assert(vsize() == v.size()); // Sizes need to match

        for (size_t i = 0; i < vsize(); ++i)
          MR::Math::Vector<double>::operator[](i) = v[i];

        return *this;

      }


      void                       Object::add_prop(const char* const name, double value) {

        if (!is_owner())
          throw Exception ("Cannot add property to a fibre object that does not own the underlying data (i.e. is a view onto part of a larger structure).");

        //Loop through properties and insert the new property in the appropriate position for alphabetical order
        size_t insert_index = 0;

        for (std::vector<const char*>::iterator prop_it = props->begin(); prop_it != props->end(); ++prop_it) {

          int compare = strcmp(name, *prop_it);

          if (!compare)
            throw Exception ("Property with the same name ('" + str(name) + "') already exists in fibre object");

          // If new property name comes before current property in alphabetical order, insert it and break loop
          if (compare < 0) {
            props->insert(prop_it,name);
            break;
          }

          ++insert_index;

        }

        //Increment state vector to hold new property
        MR::Math::Vector<double>::resize(vsize()+1);

        //If the property wasn't inserted within the existing properties append it to the end of the list
        if (insert_index == props->size()) {

          props->push_back(name);  //Add pointer to property name to properties list

          MR::Math::Vector<double>::operator[](vsize()-1) = value;         //Set property value.

        } else {

          // Shuffle the values of the existing properties after the inserted property up a position
          for (int prop_i = props->size()-1; prop_i > (int)insert_index; --prop_i)
            MR::Math::Vector<double>::operator[](bsize() + prop_i) = MR::Math::Vector<double>::operator[](bsize() + prop_i - 1);

          MR::Math::Vector<double>::operator[](bsize() + insert_index) = value; //Set property value of inserted property

        }

      }


      void                       Object::remove_prop(const char* const name, bool ignore_missing) {

        if (!is_owner())
          throw Exception ("Cannot remove property from a fibre object that does not own the underlying data (i.e. is a view onto part of a larger structure).");

        //Find iterator to property that matches given property name.
        std::vector<char const*>::iterator prop_it = props->begin();

        for (; prop_it != props->end(); ++prop_it)
          if (name == *prop_it)
            break;

        if (prop_it != props->end()) {

          //Shift property values with indices higher than the removed property down a position
          for (size_t i = (prop_it - props->begin()) + 1; i < num_props(); ++i)
            MR::Math::Vector<double>::operator[](bsize() + i-1) = MR::Math::Vector<double>::operator[](bsize() + i);

          // Remove property from list
          props->erase(prop_it);

          // Resize state vector
          MR::Math::Vector<double>::resize(vsize()-1);

        } else if (!ignore_missing)
          throw Exception ("Fibre object does not have property '" + str(name) + "' (use 'ignore_missing' parameter to ignore).");

      }


      void                       Object::clear_props() {

        if (!is_owner())
          throw Exception ("Cannot clear properties from a fibre object that does not own the underlying data (i.e. is a view onto part of a larger structure).");

        // Resize state vector to remove appended properties
        MR::Math::Vector<double>::resize(bsize());

        // Clear properties list.
        props->clear();

      }


      //FIXME: Change this to template function and make it so that only the properties relevant for the current object are copied.
      void                       Object::copy_all_props(const Object& obj) {

        for (size_t prop_i = 0; prop_i < obj.num_props(); ++prop_i) {

          const char* key = obj.prop_key(prop_i);

          if (!has_prop(key))
            add_prop(key);

          prop(key) = obj.prop(prop_i);

        }

      }



      Object&                    Object::clear() {

        //Check to see if has props because if it does they will be retained (and MR::Math::Vector<double>::resize(0)
        //will throw an exception.
        if (num_props()) {

          for (size_t prop_i = 0; prop_i < num_props(); ++prop_i)
            MR::Math::Vector<double>::operator[](prop_i) = MR::Math::Vector<double>::operator[](bsize() + prop_i);

          MR::Math::Vector<double>::resize(num_props());

        } else
          MR::Math::Vector<double>::clear();

        sze = 0;

        return *this;

      }

      Object&                    Object::clear(const std::vector<const char*>& properties) {

        if (!props_match(&properties)) {

          clear_props();
          for (size_t prop_i = 0; prop_i < properties.size(); ++prop_i)
            add_prop(properties[prop_i]);

        }

        clear();

        return *this;

      }



      size_t                     Object::prop_index(const char* const name) const {

        // Find index of property
        size_t prop_index = 0;

        for (; prop_index < num_props(); ++prop_index)
          if (props->operator[](prop_index) == name)
            break;

        //Check to see if property was found.
        if (prop_index == num_props())
          throw Exception ("Fibre object does not have property '" + str(name) + "', it needs to be added by the 'add_prop' command beforehand.\n"
                           "Note that the comparison is performed by comparing pointer addresses, so the same (const char*) needs to be used throughout.");

        return prop_index;

      }


      void                       Object::resize(size_t new_size, size_t row_size, double fill_value) {

        if (!is_owner())
          throw Exception ("Cannot resize fibre object that does not own the underlying data (i.e. is a view onto part of a larger structure).");

        size_t old_bsize = bsize();
        size_t new_bsize = row_size * new_size;

        assert(old_bsize == row_size * sze);

        if (new_bsize > old_bsize) {

          //Extend state vector to hold the new base state data
          MR::Math::Vector<double>::resize(new_bsize + num_props());

          //Shift the property values out to their new positions
          for (int prop_i = num_props()-1; prop_i >= 0; --prop_i)
            MR::Math::Vector<double>::operator[](new_bsize + prop_i) = MR::Math::Vector<double>::operator[](old_bsize + prop_i);

          for (size_t i = old_bsize; i < new_bsize; ++i)
            MR::Math::Vector<double>::operator[](i) = fill_value;

        } else if (new_bsize < old_bsize) {

          //Shift the property values in to their new positions
          for (size_t prop_i = 0; prop_i < num_props(); ++prop_i)
            MR::Math::Vector<double>::operator[](new_bsize + prop_i) = MR::Math::Vector<double>::operator[](old_bsize + prop_i);

          //Resize the underlying vector.
          MR::Math::Vector<double>::resize(new_bsize + num_props());

        }

        sze = new_size;

      }



      Object&                    Object::set(double value) {

        for (uint i = 0; i < vsize(); ++i)
          this->MR::Math::Vector<double>::operator[](i) = value;

        return *this;

      }


      Object&                    Object::push_back(const Coord& c) {

        if (!is_owner())
          throw Exception ("Cannot push back to a fibre object that does not own the underlying data (i.e. is a view onto part of a larger structure).");

        //Add room for the new coord
        resize(size() + 1);

        //Set the newly added coord.
        operator[](size()-1) = c;

        return *this;

      }


      Object&                    Object::invert() {

        for (uint i = 0; i < vsize(); ++i)
          MR::Math::Vector<double>::operator[](i) = -MR::Math::Vector<double>::operator[](i);

        return *this;

      }


      Coord                      Object::left_product(const MR::Math::Vector<double>& row_vector) const {

        assert(row_vector.size() == size());

        Coord product(0,0,0);

        for (size_t i = 0; i < size(); ++i)
          product += operator[](i) * row_vector[i];

        return product;

      }


      std::string                Object::load_matlab_str(const std::string& location, double scale) {

        std::string output;

        if (!location.size())
          output = Math::matlab_str(Fibre::Tractlet::Set()); //FIXME: This will be a little confusing
        else if (File::has_or_txt_extension<Fibre::Strand>(location)) {
          Fibre::Strand::Set strands (location);
          strands *= scale;
          output = strands.matlab_str();

        } else if (File::has_or_txt_extension<Fibre::Tractlet>(location)) {
          Fibre::Tractlet::Set tractlets (location);
          tractlets *= scale;
          output = tractlets.matlab_str();

        } else if (File::extension(location) == "sta") {

          MR::Math::Vector<double> state (location);
          state *= scale;
          output = Math::matlab_str(state);

        } else
          throw Exception ("Unrecognised file extension '" + File::extension(location) + "'.");

        return output;

      }

      std::ostream& operator<< (std::ostream& stream, const Object& fibre) {

        stream << std::endl;

        for (uint prop_i = 0; prop_i < fibre.num_props(); ++prop_i)
          stream << fibre.prop_key(prop_i) << ": " << fibre.prop(prop_i) << std::endl;

        for (size_t i = 0; i < fibre.size(); i++)
          stream << fibre[i] << std::endl;

        return (stream);

      }



	  }

	}

}


