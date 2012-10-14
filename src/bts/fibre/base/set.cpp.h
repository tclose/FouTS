/*
    Copyright 2010 Brain Research Institute, Melbourne, Australia

    Written by Thomas G Close on Jun 3, 2010.

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

#ifndef ___bts_fibre_base_set_cpp_h__
#define ___bts_fibre_base_set_cpp_h__

#include "bts/common.h"

#include "bts/fibre/base/set.h"
#include "bts/coord.h"

#include "bts/fibre/base/writer.h"
#include "bts/fibre/base/reader.h"


namespace BTS {

  namespace Fibre {

    namespace Base {

      template <typename T> const char*                Set<T>::BASE_INTENSITY_PROP = "base_intensity";
      template <typename T> const Coord                Set<T>::FILE_SEPARATOR = Triple<double>(-INFINITY, NAN, INFINITY);
      template <typename T> const std::string          Set<T>::PROPS_FILE_PREAMBLE = "%%% Extended Properties File %%% - keys: ";
      template <typename T> const std::string          Set<T>::BUNDLE_INDEX_EPROP = "bundle_index";
      template <typename T> const std::string          Set<T>::NUM_ELEMS_EXT_PROP = "__num_elems__";
      template <typename T> const std::string          Set<T>::ELEM_SIZE_EXT_PROP = "__elem_degree__";
      template <typename T> const std::string          Set<T>::ELEM_VSIZE_EXT_PROP = "__elem_row_size__";


      template <typename T>                            Set<T>::Set (const Set<T>& set)
         : Object(set),
           vrows(set.vrows),
           rsize(set.rsize),
           elem_dgree(set.elem_dgree),
           row_ends(set.row_ends ? (size_t*)malloc(sizeof(size_t) * set.size()) : 0),
           elem_dgrees(set.elem_dgrees ? (size_t*)malloc(sizeof(size_t) * set.size()) : 0),
           elem_props(new std::vector<const char*>(*set.elem_props)) {

        if (size() && var_elem_degrees()) {

          if (!elem_dgrees)
            throw Exception ("Allocation of elem_dgrees (" + str(sizeof(size_t) * size()) + ") failed." );

          if (!row_ends)
            throw Exception ("Allocation of row_ends (" + str(sizeof(size_t) * size()) + ") failed." );

          memcpy(row_ends,set.row_ends, sizeof(size_t) * size());
          memcpy(elem_dgrees,set.elem_dgrees, sizeof(size_t) * size());

        }

        ext_props = set.ext_props ? new std::map<std::string,std::string>(*set.ext_props) : 0;
        ext_elem_prop_keys = set.ext_elem_prop_keys ? new std::vector<std::string>(*set.ext_elem_prop_keys) : 0;
        ext_elem_prop_defaults = set.ext_elem_prop_defaults ? new std::vector<std::string>(*set.ext_elem_prop_defaults) : 0;
        ext_elem_prop_values = set.ext_elem_prop_values ? new std::vector<std::vector<std::string> >(*set.ext_elem_prop_values) : 0;

      }


      template <typename T> std::vector<std::string>   Set<T>::extend_prop_keys() const {

        std::vector<std::string> prp_keys;

        if (ext_props)
          for (std::map<std::string,std::string>::iterator prop_it = ext_props->begin(); prop_it != ext_props->end(); ++prop_it)
            prp_keys.push_back(prop_it->first);

        return prp_keys;

      }


      template <typename T>  bool                      Set<T>::elem_props_match(const std::vector<const char*>* elem_properties) const {

        bool match;

        if (elem_props == elem_properties) // Point to the same properties
          match = true;
        else if (!elem_props->size() && !elem_properties->size()) // Both contain no properties
          match = true;
        else if (num_elem_props() != elem_properties->size()) // Are of different size
          match = false;
        else // Compare the property pointers (relies on the fact that properties are sorted by order of memory location when they are inserted).
          match = !memcmp(&(elem_props->operator[](0)),&(elem_properties->operator[](0)), num_elem_props());

        return match;

      }



      template <typename T> size_t                     Set<T>::initial_vsize(size_t size, size_t* elem_degrees, size_t elem_row_size, size_t num_elem_props, size_t num_props) {

        size_t vsize = 0;

        for (size_t i = 0; i < size; ++i)
          vsize += elem_degrees[i] * elem_row_size + num_elem_props;

        vsize += num_props;

        return vsize;

      }



      template <typename T> void                       Set<T>::set_extend_props(const std::map<std::string,std::string>& extended_props) {

        if (ext_props)
          delete ext_props;

        ext_props = new std::map<std::string,std::string>(extended_props);

      }


      template <typename T> void                       Set<T>::set_extend_prop(const std::string& prop, const std::string& value) {

        if (!ext_props)
          ext_props = new std::map<std::string,std::string>();

        (*ext_props)[prop] = value;

      }



      template <typename T> Set<T>&                    Set<T>::operator=(const Set<T>& set) {

        // If owner recreate all structures (if not check to see if all intrinsic properties and sizes are the same.
        if (is_owner()) {

          MR::Math::Vector<double>::operator=(set);

          *props = *set.props;
          *elem_props = *set.elem_props;

          // Copy across the variable row ends and elem sizes if present
          if (set.var_elem_degrees() && set.sze) {

            if (!var_elem_degrees()) {
              row_ends = (size_t*)malloc(sizeof(size_t) * set.sze);
              elem_dgrees = (size_t*)malloc(sizeof(size_t) * set.sze);
            } else {
              row_ends = (size_t*)realloc(row_ends,sizeof(size_t) * set.sze);
              elem_dgrees = (size_t*)realloc(elem_dgrees,sizeof(size_t) * set.sze);
            }

            if (!row_ends)
              throw Exception ("Allocation of row_ends (" + str(sizeof(size_t) * set.sze) + ") failed.");

            if (!elem_dgrees)
              throw Exception ("Allocation of elem_dgrees (" + str(sizeof(size_t) * set.sze) + ") failed.");

            memcpy(row_ends, set.row_ends, sizeof(size_t) * set.sze);
            memcpy(elem_dgrees, set.elem_dgrees, sizeof(size_t) * set.sze);

          } else if (var_elem_degrees() && sze) {

            free(row_ends);
            free(elem_dgrees);
            row_ends = 0;
            elem_dgrees = 0;

          }

          //Now that all the (possibly) variable row sizes componenets are copied across copy the fixed row size components.
          sze = set.sze;
          vrows = set.vrows;
          rsize = set.rsize;
          elem_dgree = set.elem_dgree;

          // Copy across the element properties if present
          if (set.ext_elem_prop_keys) {

            assert(set.ext_elem_prop_defaults);
            assert(set.ext_elem_prop_values);

            if (!ext_elem_prop_keys) {

              assert(!ext_elem_prop_defaults);
              assert(!ext_elem_prop_values);

              ext_elem_prop_keys = new std::vector<std::string>(*set.ext_elem_prop_keys);
              ext_elem_prop_defaults = new std::vector<std::string>(*set.ext_elem_prop_defaults);
              ext_elem_prop_values = new std::vector< std::vector<std::string> >(*set.ext_elem_prop_values);

            } else {

              assert(ext_elem_prop_defaults);
              assert(ext_elem_prop_values);

              ext_elem_prop_keys->operator=(*set.ext_elem_prop_keys);
              ext_elem_prop_defaults->operator=(*set.ext_elem_prop_defaults);
              ext_elem_prop_values->operator=(*set.ext_elem_prop_values);

            }

          } else {

            assert(!set.ext_elem_prop_defaults);
            assert(!set.ext_elem_prop_values);

            if (ext_elem_prop_keys) {

              free(ext_elem_prop_keys);
              free(ext_elem_prop_defaults);
              free(ext_elem_prop_values);
              ext_elem_prop_keys = 0;
              ext_elem_prop_defaults = 0;
              ext_elem_prop_values = 0;

            }

          }


          //Copy accross the set-wide properties if present
          if (set.ext_props) {
            if (!ext_props)
              ext_props = new std::map<std::string,std::string>(*set.ext_props);
            else
              *ext_props = *set.ext_props;

          } else if (ext_props) {
              free(ext_props);
              ext_props = 0;
          }



        } else {

          // Might consider making these Exceptions instead of assertions, but this function probably should be as fast as possible.
          assert(size() == set.size()); // Sizes need to match when fibre object is not owner of data.
          assert(props_match(set)); // Properties need to match exactly when fibre object is not owner of data.
          assert(rsize == set.rsize); // Row sizes need to be the same (NB: when variable they should both be 0).
          assert(elem_dgree == set.elem_dgree);
          assert(!var_elem_degrees() || memcmp(row_ends, set.row_ends, size() * sizeof(size_t))); // All row sizes must match for variable row sizes case.
          assert(!var_elem_degrees() || memcmp(elem_dgrees, set.elem_dgrees, size() * sizeof(size_t))); // All row sizes must match for variable row sizes case.
          assert(var_elem_degrees() == set.var_elem_degrees());

          const MR::Math::Vector<double>& setvec = set;

          for (size_t i = 0; i < vsize(); ++i)
            MR::Math::Vector<double>::operator[](i) = setvec[i];

        }

        return *this;

      }



      template <typename T> void                       Set<T>::resize(size_t new_size, double fill_value, size_t new_elem_degree, size_t new_elem_vsize) {

        //FIXME: Needs some restructuring as functionality from Object was moved in here.

        if (!is_owner())
          throw Exception ("Cannot resize fibre object that does not own the underlying data (i.e. is a view onto part of a larger structure).");

        // Preserve the old size and base size.
        size_t old_size = size();
        size_t old_bsize = bsize();

        //Calculate the the new base size
        size_t new_bsize;

        //If new rsize is set (i.e. non-variable row sizes) simply set the bsize to new_size * rsize
        if (!var_elem_degrees()) {

          if (new_elem_degree && (new_elem_degree != elem_dgree))
            throw Exception ("The the new element size does not match that of set, and the set elements are not variable (see free_elem_degrees()).");

          if (new_elem_vsize && (new_elem_vsize != rsize))
            throw Exception ("The the new element vsize does not match that of set, and the set elements are not variable (see free_elem_degrees()).");

          new_bsize = new_size * rsize;

        } else {

          //If variable row and the new size is smaller the previous simply take the row end at that index to be the new bsize().
          if (new_size <= size())
            new_bsize = new_size ? row_end(new_size-1) : 0;

          //If variable row and the new size is greater than the previous take the multiple of the new vsize.
          else
            new_bsize = (size() ? row_end(size()-1) : 0) + (new_size - size()) * new_elem_vsize;

        }

        //Resize the underlying vector
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


        //If variable rows then extend the variable rows and elem sizes variables.
        if (var_elem_degrees()) {

          if (new_size) {

            if (row_ends)
              row_ends = (size_t*)realloc(row_ends, sizeof(size_t) * new_size);
            else
              row_ends = (size_t*)malloc(sizeof(size_t) * new_size);

            if (!row_ends)
              throw Exception ("Failed to allocate array for row_ends (" + str(sizeof(size_t) * new_size) + " bytes requested).");

            if (elem_dgrees)
              elem_dgrees = (size_t*)realloc(elem_dgrees, sizeof(size_t) * new_size);
            else
              elem_dgrees = (size_t*)malloc(sizeof(size_t) * new_size);

            if (!elem_dgrees)
              throw Exception ("Failed to allocate array for elem_dgrees (" + str(sizeof(size_t) * new_size) + " bytes requested).");


            size_t row_start = old_bsize;
            for (size_t row_i = old_size; row_i < new_size; ++row_i)
              row_ends[row_i] = row_start = row_start + new_elem_vsize;

            //If new_size > old_size, fill with the new elem sizes
            for (size_t row_i = old_size; row_i < new_size; ++row_i)
              elem_dgrees[row_i] = new_elem_degree;

          } else {

            if (row_ends) {
              assert(elem_dgrees);

              free(row_ends);
              free(elem_dgrees);
              elem_dgrees = 0;
              row_ends = 0;
            }
          }

        }

        //Update the size of the set
        sze = new_size;

        //Extend the element properties.
        if (num_extend_elem_props()) {

          for (size_t key_i = 0; key_i < num_extend_elem_props(); ++key_i)
            ext_elem_prop_values->operator[](key_i).resize(new_size, ext_elem_prop_defaults->operator[](key_i));

        }

      }



      template <typename T> T                          Set<T>::push_back(const T& elem) {

        rand_double();

        if (!elem.props_match(elem_props))
          throw Exception ("Element properties (" + str(elem.prop_keys()) + ") do not match that of set " + str(elem_props) + ".");

        if (rsize && (elem.vsize() != rsize))
          throw Exception ("vize() of element (" + str(elem.vsize()) + ") does not match that of the other elements in set (" + str(rsize) + "). Call 'free_elem_degree()' on set if this is intentional.");

        resize(size() + 1, NAN, elem.degree(), elem.vsize());

        T new_elem = operator[](size()-1);

        new_elem = elem;

        return new_elem;

      }


      template <typename T> T                          Set<T>::push_back(const T& elem, bool with_props) {

        if (with_props)
          return push_back(elem);

        else {

          if (!var_elem_degrees() && (elem.bsize() != (rsize - num_elem_props())))
            throw Exception ("bize() of element (" + str(elem.bsize()) + ") does not match that of the other elements in set (" + str(rsize - num_elem_props()) + "). Call 'free_elem_degree()' on set if this is intentional.");

          resize(size() + 1, NAN, elem.degree(), elem.bsize() + num_elem_props());

          T set_elem = operator[](size()-1);

          //Only set base vector part of element (the part without properties).
          set_elem.sub(0, elem.bsize()) = elem.bvector();

          return set_elem;
        }

      }


      template <typename T> Set<T>&                    Set<T>::clear() {

        Object::clear();

        if (row_ends) {
          free(row_ends);
          row_ends = 0;
        }

        if (elem_dgrees) {
          free(elem_dgrees);
          elem_dgrees = 0;
        }

        if (ext_props) {
          delete ext_props;
          ext_props = 0;
        }
        if (ext_elem_prop_keys) {
          delete ext_elem_prop_keys;
          delete ext_elem_prop_defaults;
          delete ext_elem_prop_values;
          ext_elem_prop_keys = 0;
          ext_elem_prop_defaults = 0;
          ext_elem_prop_values = 0;
        }

        return *this;
      }


      template <typename T> Set<T>&    Set<T>::clear(const std::vector<const char*>& properties,
                                                     const std::vector<const char*>& elem_properties) {


        if (!Object::props_match(&properties)) {

          clear_props();
          for (size_t prop_i = 0; prop_i < properties.size(); ++prop_i)
            add_prop(properties[prop_i]);

        }

        if (!elem_props_match(&elem_properties)) {

          clear_elem_props();
          for (size_t prop_i = 0; prop_i < elem_properties.size(); ++prop_i)
            add_elem_prop(properties[prop_i]);

        }

        return *this;

      }


      template <typename T> std::vector<std::string>   Set<T>::extend_elem_prop_keys() const {

        std::vector<std::string> header;

        if (ext_elem_prop_keys)
          header = *ext_elem_prop_keys;
        else
          header = std::vector<std::string>();

        std::sort(header.begin(),header.end());

        return header;

      }


      template <typename T> std::vector<const char*>   Set<T>::elem_prop_keys() const {

        std::vector<const char*> header = *elem_props;

        std::sort(header.begin(),header.end(), cstr_sort);

        return header;

      }


      template <typename T> template <typename U> void Set<T>::copy_extend_elem_props(const Set<U>& set)                       {

        if (set.size() != this->size())
          throw Exception ("Size of set (" + str(this->size()) + ") does not match size of source set (" + str(set.size()) + ") for properties copy.");

        clear_extend_elem_props();
        append_extend_elem_props(set);

      }


      template <typename T> template <typename U> void Set<T>::append_extend_elem_props(const Set<U>& set) {

        if (set.size() != this->size())
          throw Exception ("Size of set (" + str(this->size()) + ") does not match size of source set (" + str(set.size()) + ") for properties append.");

        if (set.ext_elem_prop_keys) {

          assert(set.ext_elem_prop_defaults);
          assert(set.ext_elem_prop_values);

          if (!ext_elem_prop_keys) {

            ext_elem_prop_keys = new std::vector<std::string>(*set.ext_elem_prop_keys);
            ext_elem_prop_defaults = new std::vector<std::string>(*set.ext_elem_prop_defaults);
            ext_elem_prop_values = new std::vector< std::vector<std::string> >(*set.ext_elem_prop_values);

          } else {

            assert(ext_elem_prop_defaults);
            assert(ext_elem_prop_values);

            for (size_t prop_i = 0; prop_i < set.ext_elem_prop_keys->size(); ++prop_i) {

              std::string& key = set.ext_elem_prop_keys->operator[](prop_i);

              if (!has_extend_elem_prop(key))
                add_extend_elem_prop(key, set.ext_elem_prop_defaults->operator[](prop_i));

              for (size_t row_i = 0; row_i < size(); ++row_i)
                set_extend_elem_prop(key, set.ext_elem_prop_values->operator[](prop_i)[row_i], row_i);

            }

          }

        } else if (ext_elem_prop_keys) {

          assert(!set.ext_elem_prop_defaults);
          assert(!set.ext_elem_prop_values);
          assert(ext_elem_prop_defaults);
          assert(ext_elem_prop_values);

          free(ext_elem_prop_keys);
          free(ext_elem_prop_defaults);
          free(ext_elem_prop_values);
          ext_elem_prop_keys = 0;
          ext_elem_prop_defaults = 0;
          ext_elem_prop_values = 0;

        }

      }


      template <typename T> std::vector<std::string>&   Set<T>::insert_elem_prop_keys(std::vector<std::string>& header) const {

        for (uint prop_i = 0; prop_i < num_elem_props(); ++prop_i)
          header.push_back(elem_prop_key(prop_i));

        return header;

      }


      template <typename T> Set<T>&                    Set<T>::select(Set<T>& new_set, std::vector<size_t>& indices) const {

        if (!new_set.is_owner())
          throw Exception ("Cannot select fibres into provided set as it is only a view onto a larger object.");

        new_set.reset(*props,*elem_props);

        new_set.add_extend_elem_props(*this);

        for (size_t i = 0; i < indices.size(); i++) {

          if (indices[i] >= this->size())
            throw Exception ("Supplied index, " + str(indices[i]) + ", is out of range.");

          std::map<std::string,std::string> props_row = get_extend_elem_prop_row(indices[i]);

          T new_elem;

          new_elem = operator[](indices[i]);

          new_set.push_back(new_elem, props_row);

        }

        // Freeze element sizes if frozen in this set.
        if (!var_elem_degrees())
          new_set.freeze_elem_degree();

        return new_set;

      }


      template <typename T> void                       Set<T>::load (const std::string& location) {

        if (File::has_extension(location, File::TXT_FILE_EXTENSTION)) {

          typename T::TextReader reader;
          load_tpl<typename T::TextReader>(location, reader);

        } else {

          typename T::Reader reader;
          load_tpl<typename T::Reader>(location, reader);

        }


      }


      template <typename T> void                       Set<T>::clear_extend_elem_props() {

        if (ext_elem_prop_keys) { delete ext_elem_prop_keys; ext_elem_prop_keys = 0; }
        if (ext_elem_prop_defaults) { delete ext_elem_prop_defaults; ext_elem_prop_defaults = 0; }
        if (ext_elem_prop_values) { delete ext_elem_prop_values; ext_elem_prop_values = 0; }

      }


      template <typename T> template <typename V> void Set<T>::load_tpl (const std::string& location, V& reader) {

        if (!is_owner())
          throw Exception ("Cannot load a fibre set that does not own the underlying data (i.e. is a view onto part of a larger structure).");

        reader.open (location);

        //Separate out the intrinsic and extended properties read from the file header
        std::map<std::string,std::string> loaded_extend_props = reader.get_extend_props();
        std::map<const char*, double> loaded_props = extract_props<typename T::Set>(loaded_extend_props);

        // Extract just the keys for the set properties
        std::vector<const char*> loaded_prop_keys;
        for (std::map<const char*, double>::iterator prop_it = loaded_props.begin(); prop_it != loaded_props.end(); ++prop_it)
          loaded_prop_keys.push_back(prop_it->first);

        //Separate out the intrinsic and extended properties read from the element properties file header
        std::vector<std::string> extended_elem_prop_header = reader.extend_prop_keys();
        std::vector<const char*> loaded_elem_prop_keys = reader.prop_keys();

        // If the number of elements and their sizes have been stored in the file header use them to initialise the set.
        if (loaded_extend_props.count(NUM_ELEMS_EXT_PROP) && loaded_extend_props.count(ELEM_SIZE_EXT_PROP)
                                                                   && loaded_extend_props.count(ELEM_VSIZE_EXT_PROP))
          reset(to<size_t>(loaded_extend_props[NUM_ELEMS_EXT_PROP]),
                to<size_t>(loaded_extend_props[ELEM_SIZE_EXT_PROP]),
                to<size_t>(loaded_extend_props[ELEM_VSIZE_EXT_PROP]),
                loaded_prop_keys,
                loaded_elem_prop_keys);

        // Else reset the set to a variable row size set.
        else
          reset(loaded_prop_keys,
                loaded_elem_prop_keys);


        for (std::map<const char*, double>::iterator prop_it = loaded_props.begin(); prop_it != loaded_props.end(); ++prop_it)
          prop(prop_it->first) = prop_it->second;


        //Add the extended properties to the set
        for (std::vector<std::string>::iterator prop_it = extended_elem_prop_header.begin(); prop_it != extended_elem_prop_header.end(); ++prop_it)
          add_extend_elem_prop(*prop_it, "");

        //Loop through the reader and append the loaded elements.
        T elem;
        std::map<std::string,std::string> properties_row;
        size_t count = 0;

        while (reader.next(elem, properties_row)) {

          if (elem.degree() == 0) {
            std::cout << std::endl << "WARNING!! Omitting element " << count << " as it's size is 0." << std::endl;
            continue;
          }

          // If the next element does not fit the fixed element size loaded from the header, remove the remaining unfilled
          // elements and free the element size for the set.
          if (!var_elem_degrees() && (elem.degree() != elem_degree())) {
            resize(count);
            free_elem_degree();
          }


          //If estimated size is smaller than actual size append a new element to the end
          if (count >= size())
            push_back(elem, properties_row);

          //Otherwise just set the already initialised element.
          else {
            operator[](count) = elem;
            set_extend_elem_prop_row(properties_row, count);
          }

          ++count;

        }

        if (!size())
          std::cout << "WARNING!! No elements loaded from file '" << location << "'.";

        //If estimated size is larger than actual size, resize the set to match the actual size.
        if (size() > count)
          resize(count);

        // Attempt to freeze the row sizes for the loaded set, ignoring failure.
        freeze_elem_degree(true);

      }


      template <typename T> void                       Set<T>::save (const std::string& location) const {

        if (File::has_extension(location, File::TXT_FILE_EXTENSTION)) {

          typename T::TextWriter writer;
          save_tpl<typename T::TextWriter>(location, writer);

        } else {

          typename T::Writer writer;
          save_tpl<typename T::Writer>(location, writer);

        }


      }


      template <typename T> template <typename V> void Set<T>::save_tpl (const std::string& location, V& writer) const {

        //Copy set properties and append intrinsic properties to them.
        std::map<std::string,std::string> properties;

        if (ext_props)
          properties = *ext_props;

        //Insert the properties for the set into the header along with the extended properties
        insert_props(properties);

        //Copy element properties and append intrinsic properties to them.
        std::vector<std::string> elem_prop_header;

        if (ext_elem_prop_keys)
          elem_prop_header = *ext_elem_prop_keys;

        writer.create (location, *elem_props, elem_prop_header, properties);

        for (size_t elem_i = 0; elem_i < size(); elem_i++)
          writer.append(operator[](elem_i), get_extend_elem_prop_row(elem_i));


      }


      template <typename T> Set<T>&                    Set<T>::permute(Set<T>& permuted_set, const std::vector<size_t>& indices) const {

        if (indices.size() != this->size())
          throw Exception ("Size of indices (" + str(indices.size()) + ") does not match size of set (" + str(this->size()) + ").");

        if (permuted_set.size() != this->size())
          throw Exception ("Size of output set (" + str(permuted_set.size()) + ") does not match size of set (" + str(this->size()) + ").");

#ifndef NDEBUG
        std::vector<size_t> unique_indices (indices);

        std::unique(unique_indices.begin(), unique_indices.end());

        if (unique_indices.size() != indices.size())
          throw Exception ("Indices provided to Fibre::Base::Set<T>::permute() are not unique");
#endif

        //To avoid a bug if *this is the same as permuted_set;
        Set<T> this_copy(*this);

        for (size_t indices_i = 0; indices_i < indices.size(); indices_i++) {
          permuted_set[indices[indices_i]] = this_copy[indices_i];
        }

        return permuted_set;

      }


      template <typename T> double                     Set<T>::distance(const Set<T>& reference, std::vector<size_t>& matched_indices, BTS::Math::Munkres<double>&munkres_algorithm, MR::Math::Matrix<double>& similarity, double strands_per_acs, bool add_extra) const {


        if (reference.size() < this->size())
          throw Exception( "Size of reference Strand::Set (" + str (reference.size()) + ") must be greater than or equal to that of the Strand::Set of interest  (" + str(this->size()) + ").");

        matched_indices.resize(this->size());

        this->similarity_matrix(reference, similarity, strands_per_acs);

        munkres_algorithm.match(similarity, matched_indices);

        double dist = 0.0;

        for (size_t fibre_i = 0; fibre_i < this->size(); fibre_i++) {

          dist += this->operator[](fibre_i).distance(reference[matched_indices[fibre_i]], strands_per_acs);
  //        dist += similarity_matrix(fibre_i, matched_indices[fibre_i]);

        }

        if (add_extra)
          for (size_t ref_i = 0; ref_i < reference.size(); ++ref_i)
            //FIXME: This should not include the displacement terms only the terms that contribute to its length.
            if (std::find(matched_indices.begin(), matched_indices.end(), ref_i) == matched_indices.end())
              dist += MR::Math::sqrt(reference[ref_i].norm2());


        return dist;


      }


      template <typename T> double                     Set<T>::distance(const Set<T>& reference, double strands_per_acs, bool add_extra) const {

        std::vector<size_t> matched_indices(size());

        return distance(reference, matched_indices, strands_per_acs, add_extra);

      }


      template <typename T> double                     Set<T>::distance(const Set<T>& reference, std::vector<size_t>& matched_indices, double strands_per_acs, bool add_extra) const {

        BTS::Math::Munkres<double> munkres_algorithm(this->size(), reference.size());
        MR::Math::Matrix<double> similarity(this->size(),reference.size());

        return distance(reference, matched_indices, munkres_algorithm, similarity, strands_per_acs, add_extra);

      }


      template <typename T> void                       Set<T>::similarity_matrix(const Set<T>& reference,
                                                                                                  MR::Math::Matrix<double>& similarity,
                                                                                                  double strands_per_acs) const {

        if ((similarity.rows() != size()) || (similarity.columns() != reference.size()))
          throw Exception ("Dimensions of similarity matrix, " + str(similarity.rows()) + " X " + str(similarity.columns()) + " do not match Strand::Sets, " + str(size()) + ", " + str(reference.size()) + ".");

        for (size_t this_i = 0; this_i < size(); this_i++) {
          for (size_t reference_i = 0; reference_i < reference.size(); reference_i++)
            similarity(this_i, reference_i) = operator[](this_i).distance(reference[reference_i], strands_per_acs);
        }

      }


      template <typename T> Set<T>&                    Set<T>::smallest_distance_set(const Set<T>& reference, Set<T>& smallest) const {

        if (size() != reference.size())
          throw Exception ("Size of set (" + str(size()) + ") does not match that of reference set (" + str(reference.size()) + ").");

        std::vector<size_t> matched_indices(size());

        distance(reference, matched_indices);

        smallest = *this;

        smallest.permute(smallest, matched_indices);

        for (size_t fibre_i = 0; fibre_i < smallest.size(); ++fibre_i)
          smallest[fibre_i] = smallest[fibre_i].smallest_distance_set(reference[fibre_i]);

        return smallest;

      }


      template <typename T> void                       Set<T>::reset_bundle_indices() {

        if (!has_extend_elem_prop(BUNDLE_INDEX_EPROP))
          add_extend_elem_prop(BUNDLE_INDEX_EPROP, "-1");

        for (size_t elem_i = 0; elem_i < size(); ++elem_i)
          set_extend_elem_prop(BUNDLE_INDEX_EPROP, str(elem_i), elem_i);

      }


      template <typename T> void                       Set<T>::add_extend_elem_prop(const std::string& key, const std::string& default_value) {

        this->init_extend_elem_props();

        //If key is not already present.
        if (!has_extend_elem_prop(key)) {

          ext_elem_prop_keys->push_back(key);
          ext_elem_prop_defaults->push_back(default_value);

          ext_elem_prop_values->push_back ( std::vector<std::string>(size(), default_value) );

        } else
          ext_elem_prop_defaults->operator[](key_index(key)) = default_value;

      }


      template <typename T> template <typename U> void Set<T>::add_extend_elem_props(const Set<U>& set) {

        //If set has extended properties
        if (set.ext_elem_prop_keys) {

          for (size_t key_i = 0; key_i < set.ext_elem_prop_keys->size(); key_i++)
            add_extend_elem_prop(set.ext_elem_prop_keys->operator[](key_i), set.ext_elem_prop_defaults->operator[](key_i));
        }

      }


      template <typename T> void                       Set<T>::remove_extend_elem_prop(const std::string& key) {


        if (has_extend_elem_prop(key)) {

          size_t index = key_index(key);

          ext_elem_prop_keys->erase(ext_elem_prop_keys->begin() + index);
          ext_elem_prop_defaults->erase(ext_elem_prop_defaults->begin() + index);
          ext_elem_prop_values->erase(ext_elem_prop_values->begin() + index);


        }

      }


      template <typename T> void                       Set<T>::set_extend_elem_prop(std::string key, std::string value, size_t elem_index) {

        if (elem_index >= size())
          throw Exception ("Element index (" + str(elem_index) + ") is out of range.");

        size_t key_i = key_index(key);

        ext_elem_prop_values->operator[](key_i)[elem_index] = value;


      }


      template <typename T> std::string                Set<T>::get_extend_elem_prop(std::string key, size_t elem_index) const {


        if (elem_index >= size())
          throw Exception ("Element index (" + str(elem_index) + ") is out of range.");

        size_t key_i = key_index(key);

        //Can't use operator[] for key_index since 'this' is const. NB: Key is guaranteed to exist in 'key_index' from the 'check_key()' function.
        std::string value = ext_elem_prop_values->operator[](key_i)[elem_index];

        return value;

      }


      template <typename T> bool                       Set<T>::has_extend_elem_prop(const std::string& key) const {

        bool found = false;

        if (ext_elem_prop_keys) {
          for (size_t key_i = 0; key_i < ext_elem_prop_keys->size(); ++key_i)
            if (ext_elem_prop_keys->operator[](key_i) == key) {
              found = true;
              break;
            }
        }

        return found;

      }


      template <typename T> std::map<std::string, std::string> Set<T>::get_extend_elem_prop_row(size_t elem_index) const {

        if (elem_index >= size())
          throw Exception ("Element index (" + str(elem_index) + ") is out of range.");

        std::map<std::string,std::string> row;

        if (ext_elem_prop_keys)
          for (size_t key_i = 0; key_i < ext_elem_prop_values->size(); ++key_i)
            row[ext_elem_prop_keys->operator[](key_i)] = ext_elem_prop_values->operator[](key_i)[elem_index];

        return row;

      }


      template <typename T> void                       Set<T>::set_extend_elem_prop_row(const std::map<std::string,std::string>& properties, size_t row_index, bool ignore_missing) {

        if (ext_elem_prop_keys) {

          std::vector<std::string> row(ext_elem_prop_keys->size());

          for (size_t key_i = 0; key_i < ext_elem_prop_keys->size(); key_i++) {

            std::map<std::string,std::string>::const_iterator key_it = properties.find(ext_elem_prop_keys->operator[](key_i));

            std::string value;

            if (key_it == properties.end()) {
              if (ignore_missing)
                value = ext_elem_prop_defaults->operator[](key_i);
              else
                throw Exception ("Property '" + ext_elem_prop_keys->operator[](key_i) + "' not found in added properties row.");

            } else
              value = key_it->second;

            //Can't use operator[] for properties since it is constant.
            this->ext_elem_prop_values->operator[](key_i)[row_index] = value;

          }

        }

      }


//      template <typename T> template <typename U> void Set<T>::copy_extend_elem_props(const Set<U>& set) {
//
//        clear_extend_elem_props();
//
//        append_extend_elem_props(set);
//
//      }
//
//
//      template <typename T> template <typename U> void Set<T>::append_extend_elem_props(const Set<U>& set) {
//
//        add_extend_elem_props(set);
//
//        for (size_t elem_i = 0; elem_i < size(); ++elem_i)
//          copy_extend_elem_prop_row(set,elem_i);
//
//      }



      template <typename T> void                        Set<T>::resize(size_t new_size, const T& elem) {

        assert(elem.props_match(elem_props));

        size_t old_size = size();

        resize(new_size, NAN, elem.degree(), elem.vsize());

        for (size_t elem_i = old_size; elem_i < new_size; ++elem_i)
          operator[](elem_i) = elem;

      }


      template <typename T> void                    	  Set<T>::reset(size_t num_elems,
                                                                  size_t elem_degree,
                                                                  size_t elem_vsize,
                                                                  const std::vector<const char*>& properties,
                                                                  const std::vector<const char*>& elem_properties,
                                                                  double fill_value) {

        if (!owner)
          throw Exception ("Cannot reset set as it is part of a larger object.");

        Object::reset(num_elems, elem_vsize, properties, fill_value);

        elem_dgree = elem_degree;
        rsize = elem_vsize;
        vrows = false;

        if (row_ends) {
          free(row_ends);
          row_ends = 0;
        }

        if (elem_dgrees) {
          free(elem_dgrees);
          elem_dgrees = 0;
        }

        *elem_props = elem_properties;

        if (ext_props) {
          free(ext_props);
          ext_props = 0;
        }

        if (ext_elem_prop_keys) {
          delete ext_elem_prop_keys;
          delete ext_elem_prop_defaults;
          delete ext_elem_prop_values;
          ext_elem_prop_keys = 0;
          ext_elem_prop_defaults = 0;
          ext_elem_prop_values = 0;
        }

      }


      template <typename T> void                    	 Set<T>::reset(const std::vector<const char*>& properties,
                                                                    const std::vector<const char*>& elem_properties,
                                                                    double fill_value) {

        if (!owner)
          throw Exception ("Cannot reset set as it is part of a larger object.");

        Object::reset(properties, fill_value);

        free_elem_degree();

        *elem_props = elem_properties;

        if (ext_props) {
          free(ext_props);
          ext_props = 0;
        }

        if (ext_elem_prop_keys) {
          delete ext_elem_prop_keys;
          delete ext_elem_prop_defaults;
          delete ext_elem_prop_values;
          ext_elem_prop_keys = 0;
          ext_elem_prop_defaults = 0;
          ext_elem_prop_values = 0;
        }

      }





      template <typename T> void                       Set<T>::erase(size_t index) {

        if (!is_owner())
          throw Exception ("Cannot erase element from fibre object that does not own the underlying data (i.e. is a view onto part of a larger structure).");

        if (index > sze)
          throw Exception ("Index (" + str(index) + ") out of range (" + str(size()) + ").");

        //Calculate the new base size without the row that is to be erased.
        size_t new_bsize = bsize() - row_size(index);

        //Shift the base data into its new positions
        for (size_t i = row_start(index); i < new_bsize; ++i)
          MR::Math::Vector<double>::operator[](i) = MR::Math::Vector<double>::operator[](i + row_size(index));

        //Shift the property values in to their new positions
        for (size_t prop_i = 0; prop_i < num_props(); ++prop_i)
          MR::Math::Vector<double>::operator[](new_bsize + prop_i) = MR::Math::Vector<double>::operator[](bsize() + prop_i);

        //Resize the underlying vector.
        MR::Math::Vector<double>::resize(new_bsize + num_props());

        //deincrement the size.
        --sze;

        //If variable row sizes remove the entry for the erased row.
        if (var_elem_degrees()) {

          for (size_t i = index; i < size(); ++i)
            row_ends[i] = row_start(i) + row_size(i+1);

          if (size()) {
            row_ends = (size_t*)realloc(row_ends, sizeof(size_t) * size());

            if (!row_ends)
              throw Exception ("Failed to allocate size (" + str(sizeof(size_t) * size()) + ") for row_ends array.");

            for (size_t i = index; i < size(); ++i)
              elem_dgrees[i] = elem_dgrees[i+1];

            elem_dgrees = (size_t*)realloc(elem_dgrees, sizeof(size_t) * size());

            if (!elem_dgrees)
              throw Exception ("Failed to allocate array for elem_dgrees (" + str(sizeof(size_t) * size()) + " bytes requested).");

          } else {
            free(row_ends);
            free(elem_dgrees);
            row_ends = 0;
            elem_dgrees = 0;
          }
        }

        //Extended properties are present also remove the corresponding row.
        if (num_extend_elem_props()) {

          for (size_t key_i = 0; key_i < num_extend_elem_props(); ++key_i) {
            std::vector<std::string>& value_row = ext_elem_prop_values->operator[](key_i);
            value_row.erase(value_row.begin() + index);
          }

        }

      }


      template <typename T> T                          Set<T>::insert(const T& element, size_t index) {

        if (!is_owner())
          throw Exception ("Cannot erase element from fibre object that does not own the underlying data (i.e. is a view onto part of a larger structure).");

        if (index > sze)
          throw Exception ("Index (" + str(index) + ") out of range (" + str(size()) + ").");

        if (!element.props_match(elem_props))
          throw Exception ("Element properties do not match that of set.");

        if (elem_dgree && elem_dgree != element.degree())
          throw Exception ("Element size does not match remainder of set (NB: Use Set::free_elem_degree() function first if this is intentional).");

        size_t old_bsize = bsize();

        //Calculate the new base size without the row that is to be erased.
        size_t new_bsize = bsize() + element.vsize();

        //Resize the underlying vector.
        MR::Math::Vector<double>::resize(new_bsize + num_props());

        //Shift the property values out to their new positions
        for (int prop_i = num_props()-1; prop_i >= 0; --prop_i)
          MR::Math::Vector<double>::operator[](new_bsize + prop_i) = MR::Math::Vector<double>::operator[](old_bsize + prop_i);

        //Shift the base data into its new positions
        for (int i = old_bsize-1; i >= (int)row_start(index); --i)
          MR::Math::Vector<double>::operator[](i + element.vsize()) = MR::Math::Vector<double>::operator[](i);

        //Increment the size.
        ++sze;

        //If variable row sizes remove the entry for the erased row.
        if (var_elem_degrees()) {

          if (size()) {
            if (row_ends)
              row_ends = (size_t*)realloc(row_ends, sizeof(size_t) * size());
            else
              row_ends = (size_t*)malloc(sizeof(size_t) * size());

            if (!row_ends)
              throw Exception ("Failed to allocate size (" + str(sizeof(size_t) * size()) + ") for row_ends array.");

            for (int i = (int)size()-1; i >= (int)index; --i)
              row_ends[i] = element.vsize() + (i ? row_ends[i-1] : 0);


            if (elem_dgrees)
              elem_dgrees = (size_t*)realloc(elem_dgrees, sizeof(size_t) * size());
            else
              elem_dgrees = (size_t*)malloc(sizeof(size_t) * size());

            if (!elem_dgrees)
              throw Exception ("Failed to allocate array for elem_dgrees (" + str(sizeof(size_t) * size()) + " bytes requested).");

            for (int i = (int)size()-1; i > (int)index; --i)
              elem_dgrees[i] = elem_dgrees[i-1];

          } else {

            if (row_ends) {
              free(row_ends);
              row_ends = 0;
            }

            if (elem_dgrees) {
              free(elem_dgrees);
              elem_dgrees = 0;
            }

          }
          elem_dgrees[index] = element.degree();

        }

        //If extended properties are present, populate the corresponding row, with the default.
        if (num_extend_elem_props()) {

          for (size_t key_i = 0; key_i < num_extend_elem_props(); ++key_i) {
            std::vector<std::string>& value_row = ext_elem_prop_values->operator[](key_i);
            value_row.insert(value_row.begin() + index, ext_elem_prop_defaults->operator[](key_i));
          }

        }

        //Finally assign the element to the space created for it.
        operator[](index) = element;

        return operator[](index);

      }


      template <typename T> size_t         	           Set<T>::key_index (const std::string& key) const {

        if (!ext_elem_prop_keys)
          throw Exception ("Property '" + key + "' was not found (set has no extended properties).");

        size_t key_i;

        for (key_i = 0; key_i < ext_elem_prop_keys->size(); ++key_i)
          if (ext_elem_prop_keys->operator[](key_i) == key)
            return key_i;

        throw Exception ("Property '" + key + "' was not found.");

      }


      template <typename T> void                       Set<T>::init_extend_elem_props() {

        if (!ext_elem_prop_keys) {

          assert(!ext_elem_prop_defaults);
          assert(!ext_elem_prop_values);

          ext_elem_prop_keys = new std::vector<std::string>();
          ext_elem_prop_defaults = new std::vector<std::string>();
          ext_elem_prop_values = new std::vector<std::vector<std::string> >();

        }

      }


      template <typename T> std::ostream&              operator<< (std::ostream& stream, const Set<T>& set) {

        for (size_t prop_i = 0; prop_i < set.num_props(); ++prop_i)
          stream << set.prop_key(prop_i) << ": " << set.prop(prop_i) << std::endl;

        for (size_t elem_i = 0; elem_i < set.size(); elem_i++) {

          stream << std::endl;
          stream << "--------------------" << std::endl;
          stream << "        " << elem_i << std::endl;
          stream << "--------------------" << std::endl;



          std::vector<std::string> header = set.extend_elem_prop_keys();

          if (set.num_extend_elem_props())
            stream << "{";

          for (size_t prop_i = 0; prop_i < header.size(); prop_i++) {
            stream << header[prop_i] << ": " << set.get_extend_elem_prop(header[prop_i], elem_i);

            if (prop_i != header.size()-1)
              stream << ", ";
          }

          if (set.num_extend_elem_props())
            stream << "}" << std::endl;

          stream << set[elem_i];

          stream << std::endl;

        }

        return (stream);
      }


      template <typename T> void                       Set<T>::add_elem_prop(const char* const name, double value, bool ignore_present) {


        if (!is_owner())
          throw Exception ("Cannot add element property to a fibre object that does not own the underlying data (i.e. is a view onto part of a larger structure).");

        //Loop through properties and insert the new property in the appropriate position for alphabetical order
        size_t insert_index = 0;

        for (std::vector<const char*>::iterator prop_it = elem_props->begin(); prop_it != elem_props->end(); ++prop_it) {

          int compare = strcmp(name, *prop_it);

          //If property is already present.
          if (!compare) {

            //If not 'ignore_present' throw and exception
            if (!ignore_present)
              throw Exception ("Element property with the same name ('" + str(name) + "') already exists in fibre object");

            //Otherwise if value was supplied set the property to that value, and then return.
            if (!isnan(value))
              for (size_t elem_i = 0; elem_i < size(); ++elem_i)
                operator[](elem_i).prop(insert_index) = value;
            return;
          }

          // If new property name comes before current property in alphabetical order, insert it and break loop
          if (compare < 0) {
            elem_props->insert(prop_it,name);
            break;
          }

          ++insert_index;

        }

        //If the property wasn't inserted withinness the existing properties append it to the end of the list
        if (insert_index == elem_props->size())
          elem_props->push_back(name);  //Add pointer to property name to properties list

        //Increment state vector to hold new property
        MR::Math::Vector<double>::resize(vsize() + size());

        //Move the set properties to their new location
        for (int prop_i = num_props()-1; prop_i >=0; --prop_i)
          MR::Math::Vector<double>::operator[](bsize() + prop_i) = MR::Math::Vector<double>::operator[](bsize() + prop_i - size());

        //Adjust the row size to account for the added property
        ++rsize;

        //The offset that accounts for the accumulated increase in row size over the whole set due to the added property
        //per row.
        size_t accum_offset = size();

        for (int elem_i = size()-1; elem_i >=0; --elem_i)
          for (int elem_elem_i = rsize-1; elem_elem_i >= 0; --elem_elem_i)

            if ((size_t)elem_elem_i == insert_index + rsize - elem_props->size()) {
              MR::Math::Vector<double>::operator[](elem_i * rsize + elem_elem_i) = value;
              --accum_offset;
            } else
              MR::Math::Vector<double>::operator[](elem_i * rsize + elem_elem_i) = MR::Math::Vector<double>::operator[](elem_i * rsize + elem_elem_i - accum_offset);

      }


      template <typename T> void                       Set<T>::remove_elem_prop(const char* const name, bool ignore_missing) {

        if (!is_owner())
          throw Exception ("Cannot remove element property from a fibre object that does not own the underlying data (i.e. is a view onto part of a larger structure).");

        //Find iterator to property that matches given property name.
        std::vector<char const*>::iterator prop_it = elem_props->begin();

        for (; prop_it != elem_props->end(); ++prop_it)
          if (name == *prop_it)
            break;

        if (prop_it != elem_props->end()) {

          size_t remove_index = prop_it - elem_props->begin();

          // Remove property from list
          elem_props->erase(prop_it);

          //Adjust the row size(s) to account for the removed property
          if (var_elem_degrees()) {
            for (size_t row_i = 0; row_i < size(); ++row_i)
              row_ends[row_i] -= row_i+1;
          } else
            --rsize;


          //The offset that accounts for the accumulated increase in row size over the whole set due to the added property
          //per row.
          size_t accum_offset = 0;

          for (size_t elem_i = 0; elem_i < size(); ++elem_i) {
            for (size_t elem_elem_i = 0; elem_elem_i < row_size(elem_i); ++elem_elem_i) {

              //When the elem_elem_i index passes the index of the removed property, increment the accumulated offset
              if (elem_elem_i == remove_index + row_size(elem_i) - num_elem_props())
                ++accum_offset;

              MR::Math::Vector<double>::operator[](row_start(elem_i) + elem_elem_i) = MR::Math::Vector<double>::operator[](row_start(elem_i) + elem_elem_i + accum_offset);

            }

            //If the property is the last in the list the accumulated offset is not incremented in the above loop and needs
            // to be done here.
            if (remove_index == num_elem_props())
              ++accum_offset;

          }

          //Move the set properties to their new location
          for (size_t prop_i = 0; prop_i < num_props(); ++prop_i)
            MR::Math::Vector<double>::operator[](bsize() + prop_i - size()) = MR::Math::Vector<double>::operator[](bsize() + prop_i);


          // Resize state vector
          MR::Math::Vector<double>::resize(vsize() - size());

        } else if (!ignore_missing)
          throw Exception ("Fibre object does not have element property '" + str(name) + "' (use 'ignore_missing' parameter to ignore).");


      }


      template <typename T> void                       Set<T>::clear_elem_props() {

        if (!is_owner())
          throw Exception ("Cannot clear properties from a fibre object that does not own the underlying data (i.e. is a view onto part of a larger structure).");

        //Adjust the row size to account for the removed property
        rsize -= elem_props->size();

        //The offset that accounts for the accumulated increase in row size over the whole set due to the added property
        //per row.
        size_t accum_offset = 0;

        for (size_t elem_i = 0; elem_i < size(); ++elem_i) {
          for (size_t elem_elem_i = 0; elem_elem_i < rsize; ++elem_elem_i)
            MR::Math::Vector<double>::operator[](elem_i * rsize + elem_elem_i) = MR::Math::Vector<double>::operator[](elem_i * rsize + elem_elem_i + accum_offset);

          accum_offset += elem_props->size();
        }


        //Move the set properties to their new location
        for (size_t prop_i = 0; prop_i < num_props(); ++prop_i)
          MR::Math::Vector<double>::operator[](bsize() + prop_i - elem_props->size() * size()) = MR::Math::Vector<double>::operator[](bsize() + prop_i);

        // Resize state vector to remove appended properties
        MR::Math::Vector<double>::resize(vsize() - elem_props->size() * size());

        // Clear properties list.
        elem_props->clear();

      }


      template <typename T> template <typename U> void Set<T>::copy_elem_props(const Set<U>& set) {

        if (size() != set.size())
          throw Exception ("Set sizes must be the same (" + str(size()) + " and " + str(set.size()) + ") to copy element properties");

        for (size_t prop_i = 0; prop_i < set.num_elem_props(); ++prop_i) {

          const char* key = set.elem_prop_key(prop_i);

          if (!has_elem_prop(key))
            add_elem_prop(key);

          for (size_t elem_i = 0; elem_i < size(); ++elem_i)
            operator[](elem_i).prop(key) = set[elem_i].prop(prop_i);

        }

      }


      template <typename T> size_t                     Set<T>::elem_prop_index(const char* const name) const {

        // Find index of property
        size_t prop_index = 0;

        for (; prop_index < num_elem_props(); ++prop_index)
          if (elem_props->operator[](prop_index) == name)
            break;

        //Check to see if property was found.
        if (prop_index == num_elem_props())
          throw Exception ("Fibre set does not have element property '" + str(name) + "', it needs to be added by the 'add_elem_prop' command beforehand. "
              "             Note that the comparison is performed by comparing pointer addresses, so the same (const char*) needs to be used throughout.");

        return prop_index;

      }


      template <typename T> void                       Set<T>::elem_resize(size_t new_elem_degree, size_t new_elem_vsize, double fill_value) {

        if (!is_owner())
          throw Exception ("Cannot resize fibre set that does not own the underlying data (i.e. is a view onto part of a larger structure).");

        if (!new_elem_degree)
          throw Exception ("Cannot change fixed element size to 0. However, if row size is 'freed' using 'free_elem_degrees' their size can be 0.");

          //Buffer the current set
        Set<T> orig(*this);

        //Get Vector version of buffer to access its MR::Math::Vector<double>::operator[]().
        MR::Math::Vector<double>& orig_vector = orig;

        //Fix the row size by removing the variable row_ends

        if (row_ends) {

          if (row_ends) {
            free(row_ends);
            row_ends = 0;
          }

          if (elem_dgrees) {
            free(elem_dgrees);
            elem_dgrees = 0;
          }

        }

        //Set the new elem and row size including element properties
        elem_dgree = new_elem_degree;
        rsize = new_elem_vsize;

        //Resize the underlying vector class
        MR::Math::Vector<double>::resize(rsize * size() + num_props());

        //Loop through each element and set the new row
        for (size_t elem_i = 0; elem_i < size(); ++elem_i) {

          //Copy across the data from the buffer that is to be included in the resized set
          for (size_t elem_elem_i = 0; elem_elem_i < std::min(rsize, orig.row_size(elem_i)) - num_elem_props(); ++elem_elem_i)
            MR::Math::Vector<double>::operator[](row_start(elem_i) + elem_elem_i) = orig_vector[orig.row_start(elem_i) + elem_elem_i];

          //Fill in the extra data in the resized set if required
          for (size_t elem_elem_i = orig.row_size(elem_i) - num_elem_props(); elem_elem_i < (rsize - num_elem_props()); ++elem_elem_i)
            MR::Math::Vector<double>::operator[](row_start(elem_i) + elem_elem_i) = fill_value;

          //Copy across the element properties
          for (size_t prop_i = 0; prop_i < num_elem_props(); ++prop_i)
            MR::Math::Vector<double>::operator[](row_end(elem_i) - num_elem_props() + prop_i) = orig_vector[orig.row_end(elem_i) - num_elem_props() + prop_i];

        }

        for (size_t prop_i = 0; prop_i < num_props(); ++prop_i)
          MR::Math::Vector<double>::operator[](bsize() + prop_i) = orig_vector[orig.bsize() + prop_i];



      }


      template <typename T> std::vector<std::string>   Set<T>::parse_extend_elem_props_line(const std::string& line) {

        std::istringstream line_stream(line);

        std::string value;
        std::vector<std::string> values_row;

        while(std::getline(line_stream, value, '\t'))
          values_row.push_back(value);

        return values_row;

      }


      template <typename T> std::vector<std::string>   Set<T>::read_extend_elem_prop_header(std::ifstream& file_in) {

        check_good(file_in);

        char file_intro[PROPS_FILE_PREAMBLE.size()+2];

        file_in.get(file_intro, PROPS_FILE_PREAMBLE.size()+1);

        if (strcmp(file_intro, PROPS_FILE_PREAMBLE.c_str()))
          throw Exception ("Required file intro '" + PROPS_FILE_PREAMBLE + "' was not found at start of extended properties file (found '" + std::string(file_intro) +"' instead).");

        std::vector<std::string> header;

        std::string key_line;

        if (!std::getline(file_in, key_line))
          throw Exception ("No ext_elem_prop_keys found in extended properties file.");

        std::istringstream line_stream(key_line);

        std::string key;

        header.clear();

        while(std::getline(line_stream, key, '\t'))
          header.push_back(key);

        return header;

      }


      template <typename T> void                       Set<T>::load_extend_elem_props(std::string location) {

        clear();

        std::ifstream in;
        in.open(location.c_str());

        *ext_elem_prop_keys = read_extend_elem_prop_header(in);

        std::string line;

        try {
          while (std::getline(in, line)) {
            std::vector<std::string> values = parse_extend_elem_props_line(line);

            for (size_t key_i = 0; key_i < num_extend_elem_props(); key_i++)
              ext_elem_prop_values->operator[](key_i).push_back(values[key_i]);

          }

        } catch (Exception& e) {
          throw Exception ("Error reading extended properties from " + location + ": " + e[0] + ".");
        }




      }


      template <typename T> void                       Set<T>::save_extend_elem_props(std::string location) {

        std::ofstream out;

        out.open(location.c_str());

        out << PROPS_FILE_PREAMBLE;

        for (std::vector<std::string>::const_iterator key_it = ext_elem_prop_keys->begin(); key_it != ext_elem_prop_keys->end(); key_it++)
          out << *key_it << '\t';

        out << std::endl;

        for (size_t row_i = 0; row_i < size(); row_i++) {

          for (size_t key_i = 0; key_i < ext_elem_prop_keys->size(); key_i++)

            out << ext_elem_prop_values->operator[](key_i)[row_i] << '\t';

          out << std::endl;

        }

      }


      template <typename T> void                       Set<T>::add_extend_elem_prop_row(const std::vector<std::string>& values_row) {

        if (num_extend_elem_props() != values_row.size())
          throw Exception ("Number of supplied values " + str(values_row.size()) + " does not match extended element properties (" + str(num_extend_elem_props()) + ").");

        for (size_t key_i = 0; key_i < num_extend_elem_props(); key_i++)
          ext_elem_prop_values->operator[](key_i).push_back(values_row[key_i]);

      }


      template <typename T> void                       Set<T>::add_extend_elem_prop_row(const std::map<std::string,std::string>& properties) {

        std::vector<std::string> row(num_extend_elem_props());

        for (size_t key_i = 0; key_i < num_extend_elem_props(); key_i++) {

          std::map<std::string,std::string>::const_iterator key_it = properties.find(ext_elem_prop_keys->operator[](key_i));

          if (key_it == properties.end())
            throw Exception ("Property '" + ext_elem_prop_keys->operator[](key_i)+ "' not found in added properties row.");

          //Can't use operator[] for properties since properties are const.
          row[key_i] = key_it->second;

        }

        for (size_t key_i = 0; key_i < num_extend_elem_props(); key_i++)
          ext_elem_prop_values->operator[](key_i).push_back(row[key_i]);

      }


      template <typename T> size_t                     Set<T>::freeze_elem_degree(bool ignore_fail) {

        //Check to see if it is already fixed
        if (var_elem_degrees()) {

          bool success = true;
          size_t fix_elem_dgree = 0;
          size_t fix_rsize = 0;

          if (size()) {
            fix_rsize = row_size(0);
            fix_elem_dgree = elem_degree(0);

            for (size_t row_i = 0; row_i < size(); ++row_i)
              if (fix_rsize != row_size(row_i) || fix_elem_dgree != elem_degree(row_i)) {
                success = false;
                break;
              }

            if (success) {
              free(row_ends);
              free(elem_dgrees);
              row_ends = 0;
              elem_dgrees = 0;
            }

          }

          if (success) {
            elem_dgree = fix_elem_dgree;
            rsize = fix_rsize;
            vrows = false;
          } else if (!ignore_fail)
            throw Exception ("Attempt to fix elem sizes failed as multiple values exist in the fibre object.");

        }

        return elem_dgree;

      }


      /*! Frees the restriction on row sizes being equal to a single member variable. This excludes the use of some
       * functions, such as 'bmatrix()' that require a fixed row size and will also take up an extra sizeof(void*)
       * of memory for each element
       *
       */
      template <typename T> void                       Set<T>::free_elem_degree() {

        //Check to see if it is already free if not allocate the vector to hold the individual element sizes
        if (!var_elem_degrees()) {

          if (size()) {

            size_t alloc_size = sizeof(size_t) * size();
            if (alloc_size)
              row_ends = (size_t*)malloc(alloc_size);

            if (!row_ends)
              throw Exception ("Could not allocate memory for variable row sizes (" + str(alloc_size) + " bytes requested).");

            //Allocate an array to hold the variable element sizes
            elem_dgrees = (size_t*)malloc(sizeof(size_t) * size());

            if (!elem_dgrees)
              throw Exception ("Could not allocate memory for variable elem sizes (" + str(sizeof(size_t) * size()) + " bytes requested).");

            //Set the variable element sizes to the current elemnt size before setting it to zero.
            size_t row_start = 0;
            for (size_t row_i = 0; row_i < size(); ++row_i) {
              row_ends[row_i] = row_start + rsize;
              row_start = row_ends[row_i];
              elem_dgrees[row_i] = elem_dgree;
            }

          } else {
            elem_dgrees = 0;
            row_ends = 0;
          }

          rsize = 0;
          elem_dgree = 0;
          vrows = true;

        }

      }



      template <typename T> template <typename U> void  Set<T>::copy_relevant_elem_props(const Base::Set<U>& source) {

        if (size() != source.size())
          throw Exception ("Source and destination set sizes must be equal to copy element props (" + str(source.size()) + " and " + str(size()) + " respectively).");

        if (source.num_elem_props()) {

          // Declare the variables used withinness the loop
          size_t prop_i = 0;
          const char* prp_key;

          while ((prp_key = T::PROPS_LIST[prop_i]) != PROPS_LIST_END) {

            std::vector<const char*>::const_iterator prop_it = std::find(source.props->begin(), source.props->end(), prp_key);

            if (prop_it != source.props->end()) {
              if (!has_prop(prp_key))
                add_elem_prop(prp_key);

              for (size_t elem_i = 0; elem_i < size(); ++elem_i)
                operator[](elem_i).prop(prp_key) = source[elem_i].prop(prp_key);

            }

            ++prop_i;
          }

        }

      }

    }

  }

}


#endif /* ___bts_fibre_set_icpp_h__ */
