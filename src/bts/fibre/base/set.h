/*
 Copyright 2010 Brain Research Institute/National ICT Australia (NICTA), Melbourne, Australia

 Written by Thomas G Close, 5/05/09.

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

#ifndef __bts_fibre_base_set_h__
#define __bts_fibre_base_set_h__

#include <assert.h>

#include "math/vector.h"

#include "bts/math/munkres.h"

#include "bts/fibre/base/object.h"
#include "bts/coord.h"
#include "bts/fibre/properties.h"

#define BASE_SET_FUNCTIONS(Derived) \
  Derived&                      resize(size_t size, double fill_value = NAN) \
    { assert(!var_elem_degrees()); Base::Set<Derived::Element>::resize(size, fill_value, 0, 0); return *this; } \
\
  Derived&                      resize(size_t size, const Derived::Element& elem) \
    { Base::Set<Derived::Element>::resize(size,elem); return *this; } \
\
  void                          select(Derived& set, const std::vector<size_t>& indices) const \
  { set.reset(*props, *elem_props); Object::copy_props(set, *this); Base::Set<Derived::Element>::select(set, indices);  } \


namespace BTS {
    
    namespace Fibre {
        
        namespace Base {
            
            template<typename T> class Set: public Object {
                    
                    template<typename U> friend class Set;

                    //Nested classes and typedefs
                public:
                    
                    typedef T Element;

                    template<typename U> class Tensor;

                    //Public const static members.
                public:
                    
                    const static Coord FILE_SEPARATOR;
                    const static char* BASE_INTENSITY_PROP;

                    const static std::string PROPS_FILE_PREAMBLE;

                    const static std::string BUNDLE_INDEX_EPROP;

                    const static std::string NUM_ELEMS_EXT_PROP;
                    const static std::string ELEM_SIZE_EXT_PROP;
                    const static std::string ELEM_VSIZE_EXT_PROP;

                protected:
                    
                    size_t initial_vsize(size_t size, size_t* elem_degrees, size_t elem_row_size,
                                         size_t num_elem_props, size_t num_props);

                    //Protected member variables;
                protected:
                    
                    bool vrows; /**< A flag that determines whether the row sizes and elem
                     sizes vary or whether they are determined by the 'rsize' and elem_dgree variablse alone */
                    size_t rsize; /**< The vsize of each element if element sizes are fixed. @see Object::vsize()*/
                    size_t elem_dgree; /**< The vsize of each element if element sizes are fixed*/
                    size_t* row_ends; /**< The end index of each element if element sizes vary. */
                    size_t* elem_dgrees; /**< The size of each element if element sizes vary. */
                    
                    std::vector<const char*>* elem_props;

                    std::map<std::string, std::string>* ext_props;

                    std::vector<std::string>* ext_elem_prop_keys;
                    std::vector<std::string>* ext_elem_prop_defaults;
                    std::vector<std::vector<std::string> >* ext_elem_prop_values;

                    //Public member functions
                public:
                    
                    /*! Initialises an empty fibre base set with a given base_size and variable row sizes
                     *
                     * @param row_size The size of the fibre object's elements
                     * @param size Size of the base fibre object
                     * @param default_value The default value with which to fill the fibre object
                     * @param properties A list of properties to be appended to the state vector.
                     */
                    Set(const std::vector<const char*>& props,
                        const std::vector<const char*>& elem_props,
                        const std::map<std::string, std::string>& extended_props)
                            : Object(0, props.size(), props), vrows(true), rsize(0), elem_dgree(0), row_ends(
                                      0), elem_dgrees(0), elem_props(
                                      new std::vector<const char*>(elem_props)), ext_props(0), ext_elem_prop_keys(
                                      0), ext_elem_prop_defaults(0), ext_elem_prop_values(0) {
                        
                        if (extended_props.size())
                            ext_props = new std::map<std::string, std::string>(extended_props);
                        
                    }
                    
                    /*! Initialises a fibre base set with a given base_size and fixed row sizes.
                     *
                     * @param row_size The size of the fibre object's elements
                     * @param size Size of the base fibre object
                     * @param default_value The default value with which to fill the fibre object
                     * @param properties A list of properties to be appended to the state vector.
                     */
                    Set(size_t size, size_t elem_degree, size_t elem_vsize,
                        const std::vector<const char*>& props,
                        const std::vector<const char*>& elem_props,
                        const std::map<std::string, std::string>& extended_props = Properties())
                            : Object(size, size * elem_vsize + props.size(), props), vrows(false), rsize(
                                      elem_vsize), elem_dgree(elem_degree), row_ends(0), elem_dgrees(
                                      0), elem_props(new std::vector<const char*>(elem_props)), ext_props(
                                      0), ext_elem_prop_keys(0), ext_elem_prop_defaults(0), ext_elem_prop_values(
                                      0) {
                        
                        assert((elem_vsize - elem_props.size()) % 3 == 0);
                        
                        if (extended_props.size())
                            ext_props = new std::map<std::string, std::string>(extended_props);
                        
                    }
                    
                    /*! Initialises a fibre base set from a view onto a larger set or tensor with fixed row sizes
                     *
                     * @param size Size of the new set
                     * @param elem_degree Fixed size of the elements of the set
                     * @param elem_row_size Row size of the elements of the set
                     * @param view View on to the fibre objects that will form the set
                     * @param props Pointer to the set properties
                     * @param elem_props Pointer to the element properties
                     */
                    explicit Set(size_t size, size_t elem_degree, size_t elem_vsize,
                                 const MR::Math::Vector<double>::View& view,
                                 std::vector<const char*>* props,
                                 std::vector<const char*>* elem_props)
                            : Object(size, view, props), vrows(false), rsize(elem_vsize), elem_dgree(
                                      elem_degree), row_ends(0), elem_dgrees(0), elem_props(
                                      elem_props), ext_props(0), ext_elem_prop_keys(0), ext_elem_prop_defaults(
                                      0), ext_elem_prop_values(0) {
                        
                        assert(size * elem_vsize + props->size() == view.size());
                        assert((elem_vsize - elem_props->size()) % 3 == 0);
                        
                    }
                    
                    /*! Initialises a fibre base set from a view onto a larger set or tensor with variable row sizes
                     *
                     * @param size Size of the new set
                     * @param elem_degrees Variable row sizes of the elements of the set
                     * @param elem_row_size Row size of the elements of the set
                     * @param view View on to the fibre objects that will form the set
                     * @param props Pointer to the set properties
                     * @param elem_props Pointer to the element properties
                     */
                    explicit Set(size_t size, size_t* element_sizes, size_t* element_ends,
                                 const MR::Math::Vector<double>::View& view,
                                 std::vector<const char*>* props,
                                 std::vector<const char*>* elem_props)
                            : Object(size, view, props), vrows(true), rsize(0), elem_dgree(0), row_ends(
                                      size ? (size_t*) malloc(sizeof(size_t) * size) : 0), elem_dgrees(
                                      size ? (size_t*) malloc(sizeof(size_t) * size) : 0), elem_props(
                                      elem_props), ext_props(0), ext_elem_prop_keys(0), ext_elem_prop_defaults(
                                      0), ext_elem_prop_values(0) {
                        
                        assert(element_ends[size-1] == view.size());
                        
                        if (size) {
                            
                            if (!row_ends)
                                throw Exception(
                                        "Failed to allocate array for row_ends (" + str(
                                                sizeof(size_t) * size)
                                        + " bytes requested).");
                            
                            if (!elem_dgrees)
                                throw Exception(
                                        "Failed to allocate array for elem_dgrees (" + str(
                                                sizeof(size_t) * size)
                                        + " bytes requested).");
                            
                            memcpy(row_ends, element_ends, sizeof(size_t) * size);
                            memcpy(elem_dgrees, element_sizes, sizeof(size_t) * size);
                            
                        }
                        
                    }
                    
                    Set(const Set& set);

                    ~Set() {
                        
                        if (is_owner()) {
                            delete elem_props;
                            
                            if (ext_props)
                                delete ext_props;
                            if (ext_elem_prop_keys)
                                delete ext_elem_prop_keys;
                            if (ext_elem_prop_defaults)
                                delete ext_elem_prop_defaults;
                            if (ext_elem_prop_values)
                                delete ext_elem_prop_values;
                            
                        }
                        
                        if (row_ends)
                            free(row_ends);
                        
                        if (elem_dgrees)
                            free(elem_dgrees);
                        
                    }
                    
                    Set& operator=(const Set& set);

                    void reset(size_t num_elems, size_t elem_degree, size_t elem_row_size,
                               const std::vector<const char*>& props,
                               const std::vector<const char*>& elem_props, double default_value =
                                       NAN);

                    void reset(const std::vector<const char*>& props,
                               const std::vector<const char*>& elem_props, double fill_value = NAN);

                    T operator[](size_t idx) {
                        return T(elem_degree(idx), sub(row_start(idx), row_end(idx)), elem_props,
                                this);
                    }
                    
                    const T operator[](size_t idx) const {
                        return T(elem_degree(idx), sub(row_start(idx), row_end(idx)), elem_props,
                                this);
                    }
                    
                    T elem(size_t idx) {
                        return operator[](idx);
                    }
                    
                    const T elem(size_t idx) const {
                        return operator[](idx);
                    }
                    
                    double base_intensity() const {
                        return has_prop(BASE_INTENSITY_PROP) ? prop(BASE_INTENSITY_PROP) : 1.0;
                    }
                    
                    /*! Checks to see if the (intrinsic) element properties match between two fibre set objects, i.e. whether they have the same
                     * property keys. Is used before some mathematical operations
                     *
                     * @param props the props to match to
                     * @return true if fibre object matches the properties exactly.
                     */
                    bool elem_props_match(const std::vector<const char*>* props) const;

                    /*! Checks to see if the (intrinsic) element properties match between two fibre set objects, i.e. whether they have the same
                     * property keys. Is used before some mathematical operations
                     *
                     * @param props the props to match to
                     * @return true if fibre object matches the properties exactly.
                     */
                    bool props_match(const std::vector<const char*>* properties,
                                     const std::vector<const char*>* elem_properties) const {
                        return Object::props_match(properties) && elem_props_match(elem_properties);
                    }
                    
                    bool props_match(const Set<T>& set) const {
                        return props_match(set.props, set.elem_props);
                    }
                    
                    void set_base_intensity(double base_intens) {
                        if (!has_prop(BASE_INTENSITY_PROP))
                            add_prop(BASE_INTENSITY_PROP, base_intens);
                        else
                            prop(BASE_INTENSITY_PROP) = base_intens;
                    }
                    
                    void remove_base_intensity() {
                        remove_prop(BASE_INTENSITY_PROP);
                    }
                    
                    bool has_var_base_intensity() const {
                        return has_prop(BASE_INTENSITY_PROP);
                    }
                    
                    double& var_base_intensity() {
                        assert(has_var_base_intensity());
                        return prop(BASE_INTENSITY_PROP);
                    }
                    
                    //! Returns either the fixed elem size if present or the variable element size if that is present.
                    size_t elem_degree(size_t idx) const {
                        return elem_dgree ? elem_dgree : (elem_dgrees ? elem_dgrees[idx] : 0);
                    }
                    
                    size_t elem_degree() const {
                        assert(!var_elem_degrees());
                        return elem_dgree;
                    }
                    
                    void elem_resize(size_t new_elem_degree, size_t new_elem_vsize,
                                     double default_value);

                    void push_back(const T& elem);

                    void push_back(const T& elem,
                                   const std::map<std::string, std::string>& elem_properties_row) {
                        push_back(elem);
                        set_extend_elem_prop_row(elem_properties_row, size() - 1, true);
                    }
                    
                    /*! Erase the element at the given index
                     *
                     * @param index The index of the element to erase
                     */
                    void erase(size_t index);

                    /*! Inserts an object at the given index
                     *
                     * @param index
                     */
                    T insert(const T& element, size_t index);

                    void append(const Set& set) {
                        for (size_t elem_i = 0; elem_i < set.size(); ++elem_i)
                            push_back(set[elem_i], set.get_extend_elem_prop_row(elem_i));
                    }
                    
                    Set<T>& clear();

                    Set<T>& clear(const std::vector<const char*>& properties,
                                  const std::vector<const char*>& elem_properties);

                    double distance(const Set& reference, double strands_per_acs = 0.0,
                                    bool add_extra = false) const;

                    /*! Returns the indices of current set in 'matched_indices' that form the closest match with the reference set.
                     *
                     */
                    double distance(const Set& reference, std::vector<size_t>& matched_indices,
                                    double strands_per_acs = 0.0, bool add_extra = false) const;

                    //TODO: create another overload that allows a custom similarity matrix to be passed but doesn't need to pass the munkres algorithm.
                    double distance(const Set& reference, std::vector<size_t>& indices,
                                    BTS::Math::Munkres<double>&munkres_algorithm,
                                    MR::Math::Matrix<double>& similarity, double strands_per_acs =
                                            0.0,
                                    bool add_extra = false) const;

                    Set& smallest_distance_set(const Set& reference, Set& smallest) const;

                    MR::Math::Matrix<double> similarity_matrix(const Set& reference,
                                                               double strands_per_acs = 0.0) const {
                        MR::Math::Matrix<double> matrix(this->size(), reference.size());
                        similarity_matrix(reference, matrix, strands_per_acs);
                        return matrix;
                    }
                    
                    void similarity_matrix(const Set& reference,
                                           MR::Math::Matrix<double>& similarity,
                                           double strands_per_acs = 0.0) const;

                    // Methods to get and manipulate element properties.
                    
                    bool has_extend_prop(const std::string& key) {
                        return ext_props && ext_props->count(key);
                    }
                    
                    std::string get_extend_prop(const std::string& key) {
                        return (ext_props && ext_props->count(key)) ? (*ext_props)[key] :
                                                                      throw Exception(
                                                                              "'" + key
                                                                              + "' key was not found in extended properties.");
                    }
                    
                    std::map<std::string, std::string> get_extend_props() const {
                        return ext_props ? *ext_props : std::map<std::string, std::string>();
                    }
                    
                    std::vector<std::string> extend_prop_keys() const;

                    void set_extend_props(const std::map<std::string, std::string>& extended_props);

                    void set_extend_prop(const std::string& prop, const std::string& value);

                    size_t num_extend_elem_props() const {
                        return ext_elem_prop_keys ? ext_elem_prop_keys->size() : 0;
                    }
                    
                    std::vector<std::string> extend_elem_prop_keys() const;

                    /*! Return a copy of the element property pointers */
                    std::vector<const char*> elem_prop_keys() const;

                    void add_extend_elem_prop(const std::string& key,
                                              const std::string& default_value);

                    template<typename U> void add_extend_elem_props(const Set<U>& set);

                    void clear_extend_elem_props();

                    void remove_extend_elem_prop(const std::string& key);

                    template<typename U> void copy_relevant_elem_props(const Base::Set<U>& set);

                    void set_extend_elem_prop(std::string key, std::string value,
                                              size_t elem_index);

                    std::string get_extend_elem_prop(std::string key, size_t elem_index) const;

                    template<typename U> U get_extend_elem_prop(std::string key,
                                                                size_t elem_index) const {
                        return this->get_extend_elem_prop<U>(key, elem_index);
                    }
                    
                    bool has_extend_elem_prop(const std::string& key) const;

                    std::map<std::string, std::string> get_extend_elem_prop_row(size_t index) const;

                    void set_extend_elem_prop_row(
                            const std::map<std::string, std::string>& properties, size_t row_index,
                            bool ignore_missing = false);

                    template<typename U> inline void copy_extend_elem_prop_row(const Set<U>& set,
                                                                               size_t elem_index,
                                                                               size_t this_index) {
                        this->set_extend_elem_prop_row(set.get_extend_elem_prop_row(elem_index),
                                this_index, true);
                    }
                    
                    template<typename U> inline void copy_extend_elem_prop_row(const Set<U>& set,
                                                                               size_t elem_index) {
                        this->copy_extend_elem_prop_row(set, elem_index, elem_index);
                    }
                    
                    template<typename U> inline void copy_extend_elem_props(const Set<U>& set);

                    template<typename U> inline void append_extend_elem_props(const Set<U>& set);

                    Set& permute(Set& permuted, const std::vector<size_t>& indices) const;

                    Set& select(Set& new_set, const std::vector<size_t>& indices) const;

                    void load(const std::string& location);

                    void save(const std::string& location) const;

                    void reset_bundle_indices();

//        protected:
                    
                    template<typename V> void load_tpl(const std::string& location, V& reader);

                    template<typename V> void save_tpl(const std::string& location,
                                                       V& writer) const;

//        protected:
                    
                    size_t key_index(const std::string& key) const;

                    void init_extend_elem_props();

                    static void check_good(std::ifstream& in) {
                        if (!in.good())
                            throw Exception("Could not read from properties file.");
                    }
                    
                    static void read_props_preamble(std::ifstream& in);

                    void add_extend_elem_prop_row(const std::vector<std::string>& values_row);

                    void add_extend_elem_prop_row(
                            const std::map<std::string, std::string>& properties);

                    std::vector<std::string> parse_extend_elem_props_line(const std::string& line);

                    std::vector<std::string> read_extend_elem_prop_header(std::ifstream& file_in);

                    void load_extend_elem_props(std::string location);

                    void save_extend_elem_props(std::string location);

                    //! Returns the number of elem_properties the object has.
                    size_t num_elem_props() const {
                        return elem_props->size();
                    }
                    
                    /*! Adds a new element property to the object.
                     *  Can only be used if the object is not part of a larger object (eg. a set)
                     *
                     * @param[in] name Pointer to the element property string (NB: comparisons are performed on the pointer itself not its
                     *                 contents so must be consitent application wide. It is implemented as a c-string instead of
                     *                 an enum so for input/output purposes only)
                     * @param[in] value The value assigned to the element property
                     */
                    void add_elem_prop(const char* const name, double value = NAN,
                                       bool ignore_present = true);

                    /*! Adds a new element property to the object.
                     *  Can only be used if the object is not part of a larger object (eg. a set)
                     *
                     * @param[in] name Pointer to the element property string (NB: comparisons are performed on the pointer itself not its
                     *                 contents so must be consitent application wide. It is implemented as a c-string instead of
                     *                 an enum so for input/output purposes only)
                     */
                    bool has_elem_prop(const char* name) const {
                        for (size_t elem_prop_i = 0; elem_prop_i < num_elem_props();
                                ++elem_prop_i) {
                            if (name == elem_props->operator[](elem_prop_i))
                                return true;
                        }
                        return false;
                    }
                    
                    /*! Adds a new element property to the object.
                     *   Can only be used if the object is not part of a larger object (eg. a set)
                     *
                     * @param[in] name The element property name
                     * @param[in] ignore_missing Set to true if it doesn't matter whether the element property is present
                     */
                    void remove_elem_prop(const char* const name, bool ignore_missing = true);

                    //! Removes all elem_properties.
                    void clear_elem_props();

                    /*! Copy all element properties from @arg set to the current set
                     *
                     * @param set
                     */
                    template<typename U> void copy_elem_props(const Set<U>& set);

                    /*! Returns the index of the given element property relative the other elem_properties
                     *
                     * @param name
                     * @return index of element property for the given objecct
                     */
                    size_t elem_prop_index(const char* const name) const;

                    /*! Access the element property at the given index for const object.
                     *
                     * @param index of element property
                     * @return reference of corresponding element property
                     */
                    const char* elem_prop_key(size_t index) const {
                        return elem_props->operator[](index);
                    }
                    
                    std::vector<std::string>& insert_elem_prop_keys(
                            std::vector<std::string>& header) const;

                    /*! Attempts to fix the row size to a single member variable (stored in 'rsize') instead of many row sizes
                     *  (stored in 'row_ends'). This allows the use of some functions, such as 'bmatrix()' that require a fixed
                     *  row size. If the number of rows are variable it will throw an exception unless the 'ignore_fail' flag is set.
                     *
                     * @param ignore_fail Ignore failure to fix the row size (eg. for use in load functions)
                     * @return Fixed row size
                     */
                    size_t freeze_elem_degree(bool ignore_fail = false);

                    /*! Frees the restriction on row sizes being equal to a single member variable. This excludes the use of some
                     * functions, such as 'bmatrix()' that require a fixed row size and will also take up an extra sizeof(void*)
                     * of memory for each element
                     *
                     */
                    void free_elem_degree();

                    /*! Resize the number of @see Triple<double> objects in the fibre object (does not effect the properties)
                     *
                     * @param size The number of @see Triple<double> objects in the new fibre object
                     * @param fill_value The value assigned to the any empty values.
                     * @param new_elem_degree The size() of the new element if the row size is variable (otherwise it is not required and will be ignored).
                     * @param new_elem_vsize The vsize() of the new element if the row size is variable (otherwise it is not required and will be ignored).
                     */
                    void resize(size_t size, double fill_value = NAN, size_t new_elem_degree = 0,
                                size_t new_vsize = 0);

                    void resize(size_t size, const T& elem);

                    //! Return a MR::Math::Matrix<double> window to the underlying state vector interpreted as a matrix (without properties).
                    MR::Math::Matrix<double> bmatrix() {
                        assert(!var_elem_degrees());
                        return MR::Math::Matrix<double>(data, size(), rsize, stride());
                    }
                    
                    /*! Used in SetReader<T> to load up the element without the properties which are set afterwards
                     *
                     * @param elem The element to push to the back of the set.
                     * @param without_props Flag to to not push back the properties of the element, which are loaded separately in the reader.
                     * @return
                     */
                    void push_back(const T& elem, bool without_props);

                protected:
                    
                    /*! Returns the row size for the given row index. If the row size is fixed it is simply returned otherwise the
                     * variable row size is looked up.
                     *
                     * @param row_index Index of the row to be looked up
                     * @return The size of the specified row
                     */
                    size_t row_size(size_t row_index) const {
                        return !var_elem_degrees() ? rsize :
                                                     (row_ends ? (!row_index ? row_ends[0] :
                                                                               row_ends[row_index] - row_ends[row_index
                                                                                       - 1]) :
                                                                 0);
                    }
                    
                    /*! Returns the start index for the given row. If 'rsize' is defined then it is simply the row index x the row
                     *  size ('rsize'). If 'row_ends' is defined then the index of the previous end is returned.
                     *
                     * @param row_index Index of the row to be looked up
                     * @return The start index of the specified row
                     */
                    size_t row_start(size_t row_index) const {
                        return !var_elem_degrees() ? row_index * rsize :
                                                     (row_ends ? (row_index ? row_ends[row_index - 1] :
                                                                              0) :
                                                                 0);
                    }
                    
                    /*! Returns the end index for the given row. If 'rsize' is defined then it is simply the row index x the row size ('rsize')
                     *
                     * @param row_index Index of the row to be looked up
                     * @return The end index of the specified row
                     */
                    size_t row_end(size_t row_index) const {
                        return !var_elem_degrees() ? (row_index + 1) * rsize :
                                                     (row_ends ? row_ends[row_index] : 0);
                    }
                    
                    /*! Test whether the object has variable row sizes or not
                     *
                     * @return True if row size is variable
                     */
                    bool var_elem_degrees() const {
                        return vrows;
                    }
                    
                    /*! Used in constructors called in-turn from Set<T> constructors to calculate the row ends from variable element
                     *  sizes and a fixed elem_row_size
                     *
                     * @param elem_degrees
                     * @param elem_row_size
                     */
                    void calculate_row_ends(size_t* elem_degrees, size_t elem_row_size,
                                            size_t num_elem_props) {
                        
                        size_t row_start = 0;
                        
                        for (size_t row_i = 0; row_i < size(); ++row_i)
                            row_ends[row_i] = row_start + elem_degrees[row_i] * elem_row_size
                                              + num_elem_props;
                        
                    }
                    
            };
            
            template<typename T> std::ostream& operator<<(std::ostream& stream, const Set<T>& set);
            
            typedef std::map<std::string, std::string> Properties;
        
        }
    
    }

}

#include "bts/fibre/base/set.cpp.h"

#endif
