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

#ifndef __bts_fibre_base_h__
#define __bts_fibre_base_h__

#include "math/vector.h"
#include "math/matrix.h"
#include "bts/triple.h"
#include "bts/coord.h"

#define LOOP(op) \
  for (size_t i = 0; i < vsize(); ++i) { \
    op \
  }

//  Derived                      base()
//  { assert(!EMPTY_PROPS.size()); return Derived(sze, bvector(), &EMPTY_PROPS); }

//! Redefines general functions inherited from Fibre::Object so they return references to the derived class
#define BASE_GENERAL_FUNCTIONS(Derived) \
\
  Derived&                      set(double value) \
    { Base::Object::set(value); return *this;} \
\
  void                          copy_props(const Base::Object& source) \
    { Base::Object::copy_props(*this, source); } \
\
  std::map<std::string,std::string> extract_and_set_props(std::map<std::string,std::string>& combined_props) \
    { return Base::Object::extract_and_set_props<Derived>(combined_props); } \
\
  Derived&                      operator=(double c) \
    { MR::Math::Vector<double>::operator=(c); return *this; } \
\
  Derived&                      operator=(const MR::Math::Vector<double>& v) \
    { Base::Object::operator=(v); return *this; }

//! Redefines mulitiplication and division functions inherited from Fibre::Object so they return references to the derived class
#define BASE_INVERSION_FUNCTIONS(Derived) \
  Derived                       invert() \
    { invert(); return *this; } \
  \
  Derived                       operator-() const \
    { Derived answer (*this); answer.invert(); return answer; } \


//! Redefines mulitiplication and division functions inherited from Fibre::Object so they return references to the derived class
#define BASE_MULT_DIVIDE_FUNCTIONS(Derived) \
  Derived&                      operator*=  (const Derived& f) \
    { assert(props_match(f)); MR::Math::Vector<double>::operator*=(f); return *this; } \
    \
  Derived&                      operator*=  (double c) \
    { MR::Math::Vector<double>::operator*=(c); return *this; } \
    \
  Derived&                      operator/=  (const Derived& f) \
    { assert(props_match(f)); MR::Math::Vector<double>::operator/=(f); return *this; } \
    \
  Derived&                      operator/=  (double c) \
    { MR::Math::Vector<double>::operator/=(c); return *this; } \
 \
  Derived                       operator* (const Derived& f) const \
    { Derived answer (*this); answer *= f; return answer;} \
    \
  Derived                       operator* (double c) const \
    { Derived answer (*this); answer *= c; return answer;} \
   \
  Derived                       operator/ (const Derived& f) const \
    { Derived answer (*this); answer /= f; return answer;} \
   \
  Derived                       operator/ (double c) const \
    { Derived answer (*this); answer /= c; return answer;} \


//! Redefines addition and subtraction functions inherited from Fibre::Object so they return references to the derived class
#define BASE_ADD_SUBTRACT_FUNCTIONS(Derived) \
  Derived&                      operator+=  (const Derived& f) \
    { assert(props_match(f)); MR::Math::Vector<double>::operator+=(f); return *this; } \
  \
  Derived&                      operator-=  (const Derived& f) \
    { assert(props_match(f)); MR::Math::Vector<double>::operator-=(f); return *this; } \
  \
  Derived                       operator+ (const Derived& f) const \
    { Derived answer (*this); answer += f; return answer;} \
  \
  Derived                       operator- (const Derived& f) const \
    { Derived answer (*this); answer -= f; return answer;} \


namespace BTS {
    
    namespace Fibre {
        
        class Strand;
        class Track;
        class Tractlet;
        
        namespace Base {
            
            template<typename T> class Set;
            template<typename T> class Tensor;
            template<typename T> class Reader;
            template<typename T> class Writer;
            
            /*! Forms the base of all fibre objects (with the exception of tensors).
             *
             * Objectd on MR::Math::Vector<double> it can be treated as a one dimensional vector by statically casting it to
             * MR::Math::Vector<double>. However, note that the operator[]() is overidden in Fibre::Object to return a
             * BTS::Triple<double>.
             *
             */
            class Object: public MR::Math::Vector<double> {
                    
                    friend class Tensor<Object> ;

                public:
                    
                    const static char* ALPHA_PROP; /**< Human-readable name of the Apparrent Connection Strength (ACS) property used for printing and saving.
                     NB: Existence checks for properties use the pointer address not
                     the character strings they refer to, so be sure to use exactly this char*
                     when using the acs property (same for all properties).  */
                    
                    //! Used to signify the end of a props array.
                    const static char* PROPS_LIST_END;

                    /*! Passed to the operator[]() when constructing a fibre object view without the properties.
                     * Would like to make it const but the properties member variable of an object needs to be non const
                     * for the cases when it owns the properties. */
                    static std::vector<const char*> EMPTY_PROPS;

                public:
                    
                    /*! Loads a BTS data file and converts the loaded fibre object into a string that is readable by the MATLAB
                     *  software package. Is used to store fibre objects used as parameters (such as step sizes) withinness fibre object
                     *  file headers.
                     *
                     * @param location Location of the file containing the fibre object data
                     * @param scale An optional scalar that scales the loaded values
                     */
                    static std::string load_matlab_str(const std::string& location, double scale =
                            1.0);

                    /*! Extracts properties for the template type T (i.e. those present in T::PROPS_LIST) from a map containing
                     *  both template properties and extended properties, such as those loaded from a properties file, and
                     *  returns them in a map.
                     *
                     * @param combined_properties map containing string versions of property keys and values amongst extended properties
                     * @return map of present property keys and their values as doubles
                     */
                    template<typename T> static std::map<const char*, double>
                    extract_props(std::map<std::string, std::string>& combined_properties);

                    /*! Extracts propert keys for the template type T (i.e. those present in T::PROPS_LIST) from a vector containing
                     *  both template property keys and extended property keys, such as those loaded from a properties file, and
                     *  returns them in a vector.
                     *
                     * @param combined_header vector containing string versions of property keys amongst extended property keys
                     * @return vector of present property keys
                     */
                    template<typename T> static std::vector<const char*> extract_props(
                            std::vector<std::string>& combined_header);

                    /*! Selects propert keys for the template type T (i.e. those present in T::PROPS_LIST) from a vector containing
                     *  both template property keys and extended property keys, such as those loaded from a properties file, and
                     *  returns them in a vector. Differs from extract_props in that it doesn't delete the properties after they
                     *  are selected. @see extract_props
                     *
                     * @param combined_header vector containing string versions of property keys amongst extended property keys
                     * @return vector of present property keys
                     */
                    template<typename T> static std::vector<const char*>
                    select_props(const std::vector<const char*>& combined_header);

                    /*! Copies properties present in the source object into the destination if they appear in the PROPS_LIST
                     *  of the destination
                     *
                     * @param destination Destination where the properties are copied to.
                     * @param source Source where the properties are copied from
                     */
                    template<typename T> void copy_props(T& destination, const Object& source);

                public:
                    
                    //TODO: Make all fibre objects implement a 'template_match' function that determines if all their fields are equivalent.
                    
                protected:
                    
                    size_t sze; /**< The size of the fibre object */
                    std::vector<const char*>* props; /**< List of keys for properties that the fibre object has (the values are appended to the state vector). */
                    
                public:
                    
                    /*! Initialises a fibre base object with a given base_size, predefind props
                     *
                     * @param size Size of the base fibre object
                     * @param v_size the total "vector size" of the object
                     * @param props A list of props to be appended to the state vector.
                     */
                    Object(size_t size, size_t v_size, const std::vector<const char*>& properties)
                            : MR::Math::Vector<double>(v_size), sze(size), props(
                                      new std::vector<const char*>(properties)) {
                        
                        assert(props);
                    }
                    
                protected:
                    
                    /*! Initialises a fibre base object 'view' from an existing fibre object. As distinct
                     *  from the first constructor the fibre object return does not 'own' the underlying data and therefore
                     *  it cannot be resized and properties cannot be added or removed. However, the existing data fields can be
                     *  freely accessed and changed see base MR::Math::Vector<double> class.
                     *
                     * @param size The number of elements in the fibre object (not including properties)
                     * @param view MR::Math::Vector<double>::View onto the underlying data
                     * @param props A list of props to be appended to the state vector.
                     * @param dummy Only needed when creating a Fibre::Base::Set<Object> template functions
                     * (which is only used for testing)
                     *
                     * @see MR::Math::Vector<T>::View
                     */
                    Object(size_t size, const MR::Math::Vector<double>::View& view,
                           std::vector<const char*>* props, const void* dummy = 0)
                            : MR::Math::Vector<double>(view), sze(size), props(props) {
                        
                        assert(props);
                    }
                    
                public:
                    
                    //! Copy constructor.
                    Object(const Object& b)
                            : MR::Math::Vector<double>(b), sze(b.sze), props(
                                      new std::vector<char const*>(*b.props)) {
                    }
                    
                    //! Deletes properties if owner of the data
                    ~Object() {
                        
                        if (owner)
                            delete props;
                        
                    }
                    
                    /*! Sets the LHS fibre object to the RHS value. If the LHS object does not own the data (i.e. is a view onto
                     * part of a larger object) both the size and properties must match
                     *
                     * @param b The RHS fibre object
                     * @return Reference to LHS object
                     */
                    Object& operator=(const Object& b);

                    /*! Sets the values of the object to the vector values in @arg v
                     *
                     * @param v Vector containing the values of the object
                     * @return Reference to LHS object
                     */
                    Object& operator=(const MR::Math::Vector<double>& v);

                    /*! Clears the fibre object and resizes it out to the provided size filling with the @arg fill_value
                     *
                     * @param size Size of the new object
                     * @param row_size Row size of the new object
                     * @param props Properties of the new object
                     * @param fill_value Value used to fill the new object
                     */
                    void reset(size_t size, size_t row_size, const std::vector<const char*>& props,
                               double fill_value = NAN);

                    /*! Clears the fibre object and frees the row size
                     *
                     * @param props Properties of the new object
                     */
                    void reset(const std::vector<const char*>& props, double fill_value = NAN);

                    /*! Returns the value of the underlying MR::Math::Vector<double>::owner field to specify if the fibre object is the
                     *  owner of the data or does not own the data (i.e. is a view onto part of a larger structure
                     *
                     * @return Whether the fibre object owns the data
                     */
                    bool is_owner() const {
                        return this->MR::Math::Vector<double>::owner;
                    }
                    
                    //! Returns the "size" of the object, i.e. the number of @see BTS::Triple<double> objects it contains
                    size_t size() const {
                        return sze;
                    }
                    
                    //! Returns the "vector size" of the object, i.e. the size of the underlying MR::Math::Vector<double> object
                    size_t vsize() const {
                        return MR::Math::Vector<double>::size();
                    }
                    
                    /*! Returns the "base size" of the object, i.e. the vector size of the underlying MR::Math::Vector<double> object
                     *  minus the size of the properties.
                     */
                    size_t bsize() const {
                        assert(MR::Math::Vector<double>::size() >= num_props());
                        return MR::Math::Vector<double>::size() - num_props();
                    }
                    
                    //! Returns the number of properties the object has.
                    size_t num_props() const {
                        return props->size();
                    }
                    
                    /*! Checks to see if the (intrinsic) properties match between two fibre objects, i.e. whether they have the same
                     * property keys. Is used before some mathematical operations
                     *
                     * @param props the props to match to
                     * @return true if fibre object matches the properties exactly.
                     */
                    bool props_match(const std::vector<const char*>* props) const;

                    /*! Checks to see if the (intrinsic) properties match between the current object and a second fibre object.
                     *  @see props_match(const std::vecor<const char*>& o)
                     *
                     * @param o Second fibre object
                     * @return true if fibre object matches the properties of the second fibre object exactly.
                     */
                    bool props_match(const Object& o) const {
                        return props_match(o.props);
                    }
                    
                    /*! Adds a new property to the object.
                     *  Can only be used if the object is not part of a larger object (eg. a set)
                     *
                     * @param[in] name Pointer to the property string (NB: comparisons are performed on the pointer itself not its
                     *                 contents so must be consitent application wide. It is implemented as a c-string instead of
                     *                 an enum so for input/output purposes only)
                     * @param[in] value The value assigned to the property
                     */
                    void add_prop(const char* const name, double value = NAN);

                    /*! Adds a new property to the object.
                     *  Can only be used if the object is not part of a larger object (eg. a set)
                     *
                     * @param[in] name Pointer to the property string (NB: comparisons are performed on the pointer itself not its
                     *                 contents so must be consitent application wide. It is implemented as a c-string instead of
                     *                 an enum so for input/output purposes only)
                     */
                    bool has_prop(const char* name) const {
                        for (size_t prop_i = 0; prop_i < num_props(); ++prop_i) {
                            if (name == props->operator[](prop_i))
                                return true;
                        }
                        return false;
                    }
                    
                    /*! Adds a new property to the object.
                     *   Can only be used if the object is not part of a larger object (eg. a set)
                     *
                     * @param[in] name The property name
                     * @param[in] ignore_missing Set to true if it doesn't matter whether the property is present
                     */
                    void remove_prop(const char* const name, bool ignore_missing = false);

                    //! Removes all properties.
                    void clear_props();

                    /*! Copies all properties from the object @arg obj to the current object
                     *
                     * @param obj
                     */
                    void copy_all_props(const Object& obj);

                    /*! Returns the index of the given property relative the other properties
                     *
                     * @param name
                     * @return index of property for the given objecct
                     */
                    size_t prop_index(const char* const name) const;

                    /*! Access the property at the given index.
                     *
                     * @param index of property
                     * @return reference of corresponding property
                     */
                    double& prop(size_t index) {
                        assert(index < num_props());
                        return MR::Math::Vector<double>::operator[](bsize() + index);
                    }
                    
                    /*! Access the property at the given index for const object.
                     *
                     * @param index of property
                     * @return reference of corresponding property
                     */
                    const double& prop(size_t index) const {
                        assert(index < num_props());
                        return MR::Math::Vector<double>::operator[](bsize() + index);
                    }
                    
                    /*! Return the name of the property at the given index
                     *
                     * @param index of property
                     * @return key for the current index
                     */
                    const char* prop_key(size_t index) const {
                        assert(index < num_props());
                        return props->operator[](index);
                    }
                    
                    /*! Return a copy of the property names vector
                     *
                     * @return a copy of the property names vector
                     */
                    std::vector<const char*> prop_keys() const {
                        return *props;
                    }
                    
                    /*! Access the property at the given index.
                     *
                     * @param name of property
                     * @return reference of corresponding property
                     */
                    double& prop(const char* const name) {
                        return prop(prop_index(name));
                    }
                    
                    /*! Access the property at the given index for const object
                     *
                     * @param name of property
                     * @return reference of corresponding property
                     */
                    const double& prop(const char* const name) const {
                        return prop(prop_index(name));
                    }
                    
                    /*! Sets all values (including intrinsic properties) to the given value.
                     *
                     * @param value The double to set all values of the fibre object to
                     */
                    Object& set(double value);

                    /*! Accesses the element at the given index
                     *
                     * @param idx the index of the element
                     * @return The fibre object at the given index
                     */
                    Coord operator[](size_t idx) {
                        return Coord(sub(idx * 3, (idx + 1) * 3));
                    }
                    
                    /*! Accesses the element at the given index
                     *
                     * @param idx the index of the element
                     * @return The fibre object at the given index
                     */
                    const Coord operator[](size_t idx) const {
                        return Coord(sub(idx * 3, (idx + 1) * 3));
                    }
                    
                    /*! Accesses the element at the given index (synonym for operator[])
                     *
                     * @param idx the index of the element
                     * @return The fibre object at the given index
                     */
                    Coord elem(size_t idx) {
                        return operator[](idx);
                    }
                    
                    /*! Accesses the element at the given index (synonym for operator[] const)
                     *
                     * @param idx the index of the element
                     * @return The fibre object at the given index
                     */
                    const Coord elem(size_t idx) const {
                        return operator[](idx);
                    }
                    
                    /*! Append another @see Coord to the the fibre object
                     *
                     * @param c A coord to set the fibre object base values to.
                     */
                    Object& push_back(const Coord& c);

                    //! Remove all base values from the fibre object (preserves properties)
                    Object& clear();

                    /*! Remove all base values from the fibre object and reset properties to those provided
                     *
                     * @param props New properties to be used for the object
                     * @return
                     */
                    Object& clear(const std::vector<const char*>& props);

                    //! Sets all values to zero
                    Object& zero() {
                        set(0.0);
                        return *this;
                    }
                    
                    //! Sets all values to Not-a-Number (NaN)
                    Object& invalidate() {
                        set(NAN);
                        return *this;
                    }
                    
                    /*! Inverts all values in a copy of the fibre object
                     *
                     * @return A copy of the fibre object that has been inverted
                     */
                    Object& invert();

                    /*! Inverts all values in a copy of the fibre object
                     *
                     * @return A copy of the fibre object that has been inverted
                     */
                    Object operator-() const {
                        Object answer(*this);
                        answer.invert();
                        return answer;
                    }
                    
                    /*! Calculates the square of the 2-norm
                     *
                     * @return The squared-2-norm
                     */
                    double norm2() const {
                        return MR::Math::norm2(*this);
                    }
                    
                    /*! Calculates the 2-norm
                     *
                     * @return The squared-2-norm
                     */
                    double norm() const {
                        return MR::Math::sqrt(norm2());
                    }
                    
                    /*! Calculates the inner_product (dot) between two fibres
                     *
                     * @param fibre The fibre object to find the inner product with
                     * @return The inner product
                     */
                    double inner_product(const Object& fibre) {
                        assert(props_match(fibre));
                        return MR::Math::dot(*this, fibre);
                    }
                    
                    /*! Treats the fibre object as a Mx3 matrix and then pre-multiplies it with a 1xM vector to produce a 1X3 triple.
                     *
                     * @param row_vector A 1xM vector (where M is the size() of the fibre object
                     * @return The left product
                     */
                    Coord left_product(const MR::Math::Vector<double>& row_vector) const;

                    BASE_MULT_DIVIDE_FUNCTIONS(Object)
                    ;

                    BASE_ADD_SUBTRACT_FUNCTIONS(Object)
                    ;

                    /*! Inserts properties present in the current object into the combined properties map provided
                     *
                     * @param props The map into which the properties of the current object will be inserted
                     * @return A reference to the passed map
                     */
                    std::map<std::string, std::string>& insert_props(
                            std::map<std::string, std::string>& props) const;

                    /*! Inserts property keys present in the current object into the combined property keys vector provided
                     *
                     * @param props The other properties to be combined with this objects properties.
                     * @return The vector containing the combined property keys
                     */
                    std::vector<std::string> insert_prop_keys(
                            const std::vector<std::string>& props) const;

                    /*! Extracts properties from combined properties that are present in T::PROPS_LIST provided and set them in
                     * this object.
                     *
                     * @param combined_properties Map containing both properties and extended properties.
                     * @return A reference to the passed properties
                     */
                    template<typename T> std::map<std::string, std::string>&
                    extract_and_set_props(std::map<std::string, std::string>& combined_properties);

                    //! Return a MR::Math::Vector<double> window onto the underlying state vector with properties
                    MR::Math::Vector<double> vector() {
                        return MR::Math::Vector<double>(data, vsize(), stride());
                    }
                    
                    //! Return a MR::Math::Vector<double> window onto the underlying state vector without properties
                    MR::Math::Vector<double> bvector() {
                        return sub(0, bsize());
                    }
                    
                    //! Return a MR::Math::Vector<double> window onto the underlying state vector without properties
                    const MR::Math::Vector<double> bvector() const {
                        return sub(0, bsize());
                    }
                    
                    //! Return a MR::Math::Matrix<double> window to the underlying state vector interpreted as a matrix (without properties).
                    MR::Math::Matrix<double> bmatrix(size_t row_size) {
                        return MR::Math::Matrix<double>(data, size(), row_size, stride());
                    }
                    
                    //! Return a MR::Math::Matrix<double> window to the underlying state vector interpreted as a matrix (without properties).
                    const MR::Math::Matrix<double> bmatrix(size_t row_size) const {
                        return MR::Math::Matrix<double>(data, size(), row_size, stride());
                    }
                    
                    /*! Resize the number of elements in the fibre object (does not effect the properties).
                     * More general form is provided here for use in Base::Set<T>
                     *
                     * @param size The number of elements in the new fibre object
                     * @param fill_value The value assigned to the any empty values.
                     */
                    void resize(size_t size, double fill_value = NAN) {
                        resize(size, 3, fill_value);
                    }
                    
                    template<typename T> friend class Set;

                    //! Only included to get template function to work (a strand will never have element properties.
                    Object& push_back(const Coord& c, bool without_properties) {
                        return push_back(c);
                    }
                    
                protected:
                    
                    /*! Resize the number of elements in the fibre object general form is provided here for use in derived classes
                     *  such as Tractlet and Set<T> classes.
                     *
                     * @param size The number of elements in the new fibre object
                     * @param fill_value The value assigned to the any empty values.
                     * @param new_row_size The size of the new rows (only applicable if row size if variable).
                     */
                    void resize(size_t size, size_t row_size, double fill_value);
                    
            };
            
            std::ostream& operator<<(std::ostream& stream, const Object& fibre);
            
            template<typename T> std::map<const char*, double> Object::extract_props(
                    std::map<std::string, std::string>& loaded_props) {
                
                std::map<const char*, double> intrinsic;
                
                // Declare the variables used withinness the loop
                size_t prop_i = 0;
                const char* prp_key;
                
                while ((prp_key = T::PROPS_LIST[prop_i]) != PROPS_LIST_END) {
                    if (loaded_props.count(prp_key)) {
                        intrinsic[prp_key] = to<double>(loaded_props[prp_key]);
                        loaded_props.erase(loaded_props.find(prp_key));
                    }
                    ++prop_i;
                }
                
                return intrinsic;
                
            }
            
            template<typename T> std::vector<const char*> Object::extract_props(
                    std::vector<std::string>& header) {
                
                std::vector<const char*> intrinsic;
                
                // Declare the variables used withinness the loop
                size_t prop_i = 0;
                const char* prp_key;
                std::vector<std::string>::iterator prop_it;
                
                while ((prp_key = T::PROPS_LIST[prop_i]) != PROPS_LIST_END) {
                    
                    std::string prp_key_str = str(prp_key);
                    
                    prop_it = std::find(header.begin(), header.end(), prp_key_str);
                    
                    if (prop_it != header.end()) {
                        intrinsic.push_back(prp_key);
                        header.erase(prop_it);
                    }
                    ++prop_i;
                }
                
                return intrinsic;
                
            }
            
            template<typename T> std::vector<const char*> Object::select_props(
                    const std::vector<const char*>& properties) {
                
                std::vector<const char*> selected;
                
                // Declare the variables used withinness the loop
                size_t prop_i = 0;
                const char* prp_key;
                
                while ((prp_key = T::PROPS_LIST[prop_i]) != PROPS_LIST_END) {
                    
                    if (std::find(properties.begin(), properties.end(), prp_key) != properties.end())
                        selected.push_back(prp_key);
                    
                    ++prop_i;
                }
                
                return selected;
                
            }
            
            template<typename T> void Object::copy_props(T& destination, const Object& source) {
                
                if (source.num_props()) {
                    
                    // Declare the variables used withinness the loop
                    size_t prop_i = 0;
                    const char* prp_key;
                    
                    while ((prp_key = T::PROPS_LIST[prop_i]) != PROPS_LIST_END) {
                        
                        std::vector<const char*>::const_iterator prop_it = std::find(
                                source.props->begin(), source.props->end(), prp_key);
                        
                        if (prop_it != source.props->end()) {
                            if (destination.has_prop(prp_key))
                                destination.prop(prp_key) = source.prop(prp_key);
                            else
                                destination.add_prop(prp_key, source.prop(prp_key));
                        }
                        
                        ++prop_i;
                    }
                    
                }
                
            }
            
            template<typename T> std::map<std::string, std::string>& Object::extract_and_set_props(
                    std::map<std::string, std::string>& combined_properties) {
                
                std::map<const char*, double> props = extract_props<T>(combined_properties);
                
                for (std::map<const char*, double>::iterator prop_it = props.begin();
                        prop_it != props.end(); ++prop_it)
                    if (!has_prop(prop_it->first))
                        add_prop(prop_it->first, prop_it->second);
                    else
                        prop(prop_it->first) = prop_it->second;
                
                return combined_properties;
                
            }
            
            inline std::map<std::string, std::string>& Object::insert_props(
                    std::map<std::string, std::string>& properties) const {
                
                for (uint prop_i = 0; prop_i < num_props(); ++prop_i)
                    properties[prop_key(prop_i)] = str(prop(prop_i));
                
                return properties;
                
            }
        
        }
    
    }

}

#undef LOOP

#endif /* __bts_fibre_object_h__ */
