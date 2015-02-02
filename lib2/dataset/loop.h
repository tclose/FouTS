/*
    Copyright 2009 Brain Research Institute, Melbourne, Australia

    Written by J-Donald Tournier, 16/10/09.

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

#ifndef __dataset_loop_h__
#define __dataset_loop_h__

#include "progressbar.h"
#include "dataset/position.h"
#include "dataset/stride.h"

namespace MR
{
  namespace DataSet
  {

    /** \defgroup loop Looping functions
      @{ */

    //! a class to loop over arbitrary numbers of axes of a DataSet
    /*! This class can be used to loop over any number of axes of one of more
     * DataSets. Its use is best illustrated with the following examples.
     *
     * If \a vox in the following example is a 3D DataSet (i.e. vox.ndim() ==
     * 3), then:
     * \code
     * float sum = 0.0;
     * Loop loop;
     * for (loop.start (vox); loop.ok(); loop.next (vox))
     *   sum += vox.value();
     * \endcode
     * is equivalent to:
     * \code
     * float sum = 0.0;
     * for (vox[2] = 0; vox[2] < vox.dim(2); ++vox[2])
     *   for (vox[1] = 0; vox[1] < vox.dim(1); ++vox[1])
     *     for (vox[0] = 0; vox[0] < vox.dim(0); ++vox[0])
     *       sum += vox.value();
     * \endcode
     * This has the advantage that the dimensionality of the DataSet does not
     * need to be known at compile-time. In other words, if the DataSet was
     * 4-dimensional, the first looping construct would correctly iterate over
     * all voxels, whereas the second one would only process the first 3D
     * volume.
     *
     * \section multiloop Looping over multiple DataSets
     * It is often required to loop over more than one DataSet of the same
     * dimensions. This is done trivially by passing any additional DataSets to
     * be looped over to both the start() and next() member functions. For
     * example, this code snippet will copy the contents of the DataSet \a src
     * into a DataSet \a dest, assumed to have the same dimensions as \a src:
     * \code
     * Loop loop;
     * for (loop.start (dest, src); loop.ok(); loop.next (dest, src))
     *   dest.value() = vox.value();
     * \endcode
     *
     * \section restrictedloop Looping over a specific range of axes
     * It is also possible to explicitly specify the range of axes to be looped
     * over. In the following example, the program will loop over each 3D
     * volume in the DataSet in turn:
     * \code
     * Loop outer (3); // outer loop iterates over axes 3 and above
     * for (outer.start (vox); outer.ok(); outer.next (vox)) {
     *   float sum = 0.0;
     *   Loop inner (0, 3); // inner loop iterates over axes 0 to 2
     *   for (inner.start (vox); inner.ok(); inner.next (vox))
     *     sum += vox.value();
     *   print ("total = " + str (sum) + "\n");
     * }
     * \endcode
     *
     * \section progressloop Displaying progress status
     * The Loop object can also display its progress as it proceeds, using the
     * appropriate constructor. In the following example, the program will
     * display its progress as it averages a DataSet:
     * \code
     * float sum = 0.0;
     *
     * Loop loop ("averaging...");
     * for (loop.start (vox); loop.ok(); loop.next (vox))
     *   sum += vox.value();
     *
     * float average = sum / float (DataSet::voxel_count (vox));
     * print ("average = " + str (average) + "\n");
     * \endcode
     * The output would look something like this:
     * \code
     * myprogram: averaging... 100%
     * average = 23.42
     * \endcode
     *
     * \sa LoopInOrder
     */
    class Loop
    {
      public:
        //! Constructor
        /*! Construct a Loop object to iterate over the axes specified. By
         * default, the Loop will iterate over all axes of the first DataSet
         * supplied to next(). If \a from_axis and \a to_axis are specified,
         * the Loop will iterate from axis \a from_axis up to but \b not
         * including axis \a to_axis. */
        Loop (size_t from_axis = 0, size_t to_axis = std::numeric_limits<size_t>::max()) :
          from_ (from_axis), to_ (to_axis), cont_ (true) { }

        //! Constructor with progress status
        /*! Construct a Loop object to iterate over the axes specified and
         * display the progress status with the specified message. By default,
         * the Loop will iterate over all axes of the first DataSet supplied to
         * next(). If \a from_axis and \a to_axis are specified, the Loop will
         * iterate from axis \a from_axis up to but \b not including axis \a
         * to_axis. */
        Loop (const std::string& message, size_t from_axis = 0, size_t to_axis = std::numeric_limits<size_t>::max()) :
          from_ (from_axis), to_ (to_axis), cont_ (true), progress_ (message, 1) { }


        //! Start the loop to iterate over a single DataSet
        /*! Start the loop by resetting the appropriate coordinates of each of
         * the specified DataSets to zero, and initialising the progress status
         * if appropriate. Note that only those axes specified in the Loop
         * constructor will have their coordinates set to zero; the coordinates
         * of all other axes will be left untouched. */
        template <class Set>
        inline void start (Set& set) {
          cont_ = true;
          for (size_t n = from_; n < max_axis(set); ++n)
            set[n] = 0;
          if (progress_)
            progress_.set_max (voxel_count (set, from_, to_));
        }
        //! Start the loop to iterate over two DataSets
        /*! \copydetails Loop::start(Set&) */
        template <class Set, class Set2>
        inline void start (Set& set, Set2& set2) {
          cont_ = true;
          for (size_t n = from_; n < max_axis(set); ++n) {
            set[n] = 0;
            set2[n] = 0;
          }
          if (progress_)
            progress_.set_max (voxel_count (set, from_, to_));
        }

        //! Start the loop to iterate over three DataSets
        /*! \copydetails Loop::start(Set&) */
        template <class Set, class Set2, class Set3>
        inline void start (Set& set, Set2& set2, Set3& set3) {
          cont_ = true;
          for (size_t n = from_; n < max_axis(set); ++n) {
            set[n] = 0;
            set2[n] = 0;
            set3[n] = 0;
          }
          if (progress_)
            progress_.set_max (voxel_count (set, from_, to_));
        }

        //! Check whether the loop should continue iterating
        /*! \return true if the loop has not completed, false otherwise. */
        bool ok () const {
          return (cont_);
        }

        //! Proceed to next iteration for a single DataSet
        /*! Advance coordinates of all specified DataSets to the next position
         * to be processed, and update the progress status if appropriate. */
        template <class Set>
        void next (Set& set) {
          next_impl (from_, set);
          ++progress_;
        }

        //! Proceed to next iteration for two DataSets
        /*! \copydetails Loop::next(Set&) */
        template <class Set, class Set2>
        void next (Set& set, Set2& set2) {
          next_impl (from_, set, set2);
          ++progress_;
        }

        //! Proceed to next iteration for three DataSets
        /*! \copydetails Loop::next(Set&) */
        template <class Set, class Set2, class Set3>
        void next (Set& set, Set2& set2, Set3& set3) {
          next_impl (from_, set, set2, set3);
          ++progress_;
        }

        //! set position along relevant axes of \a target to that of \a reference
        template <class Set, class Set2>
        void set_position (const Set& reference, Set2& target) const {
          for (size_t i = from_; i < max_axis (reference); ++i) {
            const ssize_t x = reference[i];
            target[i] = x;
          }
        }

        //! set position along relevant axes of \a targets to that of \a reference
        template <class Set, class Set2, class Set3>
        void set_position (const Set& reference, Set2& target, Set3& target2) const {
          for (size_t i = from_; i < max_axis(reference); ++i) {
            const ssize_t x = reference[i];
            target[i] = x;
            target2[i] = x;
          }
        }

        //! set position along relevant axes of \a targets to that of \a reference
        template <class Set, class Set2, class Set3, class Set4>
        void set_position (const Set& reference, Set2& target, Set3& target2, Set4& target3) const {
          for (size_t i = from_; i < max_axis(reference); ++i) {
            const ssize_t x = reference[i];
            target[i] = x;
            target2[i] = x;
            target3[i] = x;
          }
        }


      private:
        const size_t from_, to_;
        bool cont_;
        ProgressBar progress_;


        template <class Set>
        size_t max_axis (const Set& set) const {
          return std::min (set.ndim(), to_);
        }

        template <class Set>
        void next_impl (size_t axis, Set& set) {
          if (axis < max_axis(set)) {
            if (set[axis] + 1 < set.dim (axis)) ++set[axis];
            else {
              if (axis+1 == max_axis(set)) {
                cont_ = false;
              }
              else {
                next_impl (axis+1, set);
                if (cont_) set[axis] = 0;
              }
            }
          }
          else cont_ = false;
        }

        template <class Set, class Set2>
        void next_impl (size_t axis, Set& set, Set2& set2) {
          if (axis < max_axis(set)) {
            if (set[axis] + 1 < set.dim (axis)) {
              ++set[axis];
              ++set2[axis];
            }
            else {
              if (axis+1 == max_axis(set)) {
                cont_ = false;
              }
              else {
                next_impl (axis+1, set, set2);
                if (cont_) {
                  set[axis] = 0;
                  set2[axis] = 0;
                }
              }
            }
          }
          else cont_ = false;
        }

        template <class Set, class Set2, class Set3>
        void next_impl (size_t axis, Set& set, Set2& set2, Set3& set3) {
          if (axis < max_axis(set)) {
            if (set[axis] + 1 < set.dim (axis)) {
              ++set[axis];
              ++set2[axis];
              ++set3[axis];
            }
            else {
              if (axis+1 == max_axis(set)) {
                cont_ = false;
              }
              else {
                next_impl (axis+1, set, set2, set3);
                if (cont_) {
                  set[axis] = 0;
                  set2[axis] = 0;
                  set3[axis] = 0;
                }
              }
            }
          }
          else cont_ = false;
        }

    };








    //! a class to loop over arbitrary numbers and orders of axes of a DataSet
    /*! This class can be used to loop over any number of axes of one of more
     * DataSets, in any specified order. Its use is essentially identical to
     * that of the Loop class, with the difference that axes can now be
     * iterated over in any arbitrary order. This is best illustrated with the
     * following examples.
     *
     * \section strideorderloop Looping with smallest stride first
     * The looping strategy most likely to make most efficient use of the
     * memory infrastructure is one where the innermost loop iterates over the
     * axis with the smallest absolute stride, since voxels along this axis are
     * most likely to be adjacent. This is most likely to optimise both
     * throughput to and from system RAM or disk (which are typically optimised
     * for bursts of contiguous sections of memory), and CPU cache usage.
     *
     * The LoopInOrder class is designed to facilitate this. In the following
     * example, the DataSet of interest is passed as an argument to the
     * constructor, so that its strides can be used to compute the nesting
     * order for the loops over the corresponding axes. Here, we assume that
     * \a vox is a 3D DataSet (i.e. vox.ndim() == 3) with strides [ 2 -1 3 ]:
     * \code
     * float sum = 0.0;
     * LoopInOrder loop (vox);
     * for (loop.start (vox); loop.ok(); loop.next (vox))
     *   sum += vox.value();
     * \endcode
     * This is equivalent to:
     * \code
     * float sum = 0.0;
     * for (vox[2] = 0; vox[2] < vox.dim(2); ++vox[2])
     *   for (vox[0] = 0; vox[0] < vox.dim(0); ++vox[0])
     *     for (vox[1] = 0; vox[1] < vox.dim(1); ++vox[1])
     *       sum += vox.value();
     * \endcode
     *
     * \section restrictedorderloop Looping over a specific range of axes
     * It is also possible to explicitly specify the range of axes to be looped
     * over. In the following example, the program will loop over each 3D
     * volume in the DataSet in turn using the Loop class, and use the
     * LoopInOrder class to iterate over the axes of each volume to ensure
     * efficient memory bandwidth use when each volume is being processed.
     * \code
     * Loop outer (3); // outer loop iterates over axes 3 and above
     * for (outer.start (vox); outer.ok(); outer.next (vox)) {
     *   float sum = 0.0;
     *   LoopInOrder inner (vox, 0, 3); // inner loop iterates over axes 0 to 2
     *   for (inner.start (vox); inner.ok(); inner.next (vox))
     *     sum += vox.value();
     *   print ("total = " + str (sum) + "\n");
     * }
     * \endcode
     *
     * \section arbitraryorderloop Arbitrary order loop
     * It is also possible to specify the looping order explictly, as in the
     * following example:
     * \code
     * float value = 0.0;
     * std::vector<size_t> order (3);
     * order[0] = 1;
     * order[1] = 0;
     * order[2] = 2;
     *
     * LoopInOrder loop (vox, order);
     * for (loop.start (vox); loop.ok(); loop.next (vox))
     *   value = Math::exp (value - vox.value());
     * \endcode
     * This will iterate over the axes in the same order as the first example
     * above, irrespective of the strides of the DataSet.
     *
     * \section multiorderloop Looping over multiple DataSets:
     * As with the Loop class, it is possible to loop over more than one
     * DataSet of the same dimensions, by passing any additional DataSets to
     * be looped over to both the start() and next() member functions. For
     * example, this code snippet will copy the contents of the DataSet \a src
     * into a DataSet \a dest (assumed to have the same dimensions as \a src),
     * with the looping order optimised for the \a dest DataSet:
     * \code
     * LoopInOrder loop (dest);
     * for (loop.start (dest, src); loop.ok(); loop.next (dest, src))
     *   dest.value() = vox.value();
     * \endcode
     *
     * \section progressloopinroder Displaying progress status
     * As in the Loop class, the LoopInOrder object can also display its
     * progress as it proceeds, using the appropriate constructor. In the
     * following example, the program will display its progress as it averages
     * a DataSet:
     * \code
     * float sum = 0.0;
     *
     * LoopInOrder loop (vox, "averaging...");
     * for (loop.start (vox); loop.ok(); loop.next (vox))
     *   sum += vox.value();
     *
     * float average = sum / float (DataSet::voxel_count (vox));
     * print ("average = " + str (average) + "\n");
     * \endcode
     * The output would look something like this:
     * \code
     * myprogram: averaging... 100%
     * average = 23.42
     * \endcode
     *
     * \sa Loop
     */
    class LoopInOrder
    {
      public:
        //! Constructor from axes indices
        /*! Construct a LoopInOrder object to iterate over the axes specified. */
        LoopInOrder (const std::vector<size_t>& axes) :
          axes_ (axes), cont_ (true) { }

        //! Construct from axes indices with progress status
        /*! Construct a LoopInOrder object to iterate over the axes specified and
         * display the progress status with the specified message. */
        LoopInOrder (const std::vector<size_t>& axes, const std::string& message) :
          axes_ (axes), cont_ (true), progress_ (message, 1) { }

        //! Construct from DataSet strides
        /*! Construct a LoopInOrder object to iterate over the axes of \a set
         * in order of smallest stride first. If supplied, the optional
         * arguments \a from_axis and \a to_axis can be used to restrict those
         * axes that will be looped over: the Loop will then iterate from axis
         * \a from_axis up to but \b not including axis \a to_axis. */
        template <class Set>
        LoopInOrder (const Set& set, size_t from_axis = 0, size_t to_axis = std::numeric_limits<size_t>::max()) :
          axes_ (Stride::order (set, from_axis, to_axis)), cont_ (true) { }

        //! Constructor from DataSet strides with progress status
        /*! Construct a LoopInOrder object to iterate over the axes specified
         * in order of smallest stride first, and display the progress status
         * with the specified message. If supplied, the optional arguments \a
         * from_axis and \a to_axis can be used to restrict those axes that
         * will be looped over: the Loop will then iterate from axis \a
         * from_axis up to but \b not including axis \a to_axis. */
        template <class Set>
        LoopInOrder (const Set& set, const std::string& message, size_t from_axis = 0, size_t to_axis = std::numeric_limits<size_t>::max()) :
          axes_ (Stride::order (set, from_axis, to_axis)), cont_ (true), progress_ (message, 1) { }

        //! Start the loop to iterate over a single DataSet
        /*! Start the loop by resetting the appropriate coordinates of each of
         * the specified DataSets to zero, and initialising the progress status
         * if appropriate. Note that only those axes specified in the
         * LoopInOrder constructor will have their coordinates set to zero; the
         * coordinates of all other axes will be left untouched. */
        template <class Set>
        inline void start (Set& set) {
          cont_ = true;
          for (size_t n = 0; n < axes_.size(); ++n)
            set[axes_[n]] = 0;
          if (progress_)
            progress_.set_max (voxel_count (set, axes_));
        }
        //! Start the loop to iterate over two DataSets
        /*! \copydetails LoopInOrder::start(Set&) */
        template <class Set, class Set2>
        inline void start (Set& set, Set2& set2) {
          cont_ = true;
          for (size_t n = 0; n < axes_.size(); ++n) {
            set[axes_[n]] = 0;
            set2[axes_[n]] = 0;
          }
          if (progress_)
            progress_.set_max (voxel_count (set, axes_));
        }

        //! Start the loop to iterate over three DataSets
        /*! \copydetails LoopInOrder::start(Set&) */
        template <class Set, class Set2, class Set3>
        inline void start (Set& set, Set2& set2, Set3& set3) {
          cont_ = true;
          for (size_t n = 0; n < axes_.size(); ++n) {
            set[axes_[n]] = 0;
            set2[axes_[n]] = 0;
            set3[axes_[n]] = 0;
          }
          if (progress_)
            progress_.set_max (voxel_count (set, axes_));
        }

        //! Check whether the loop should continue iterating
        /*! \return true if the loop has not completed, false otherwise. */
        bool ok () const {
          return (cont_);
        }

        //! Proceed to next iteration for a single DataSet
        /*! Advance coordinates of all specified DataSets to the next position
         * to be processed, and update the progress status if appropriate. */
        template <class Set>
        void next (Set& set) {
          next_impl (0, set);
          ++progress_;
        }

        //! Proceed to next iteration for two DataSets
        /*! \copydetails LoopInOrder::next(Set&) */
        template <class Set, class Set2>
        void next (Set& set, Set2& set2) {
          next_impl (0, set, set2);
          ++progress_;
        }

        //! Proceed to next iteration for three DataSets
        /*! \copydetails LoopInOrder::next(Set&) */
        template <class Set, class Set2, class Set3>
        void next (Set& set, Set2& set2, Set3& set3) {
          next_impl (0, set, set2, set3);
          ++progress_;
        }

        //! set position along relevant axes of \a target to that of \a reference
        template <class Set, class Set2>
        void set_position (const Set& reference, Set2& target) const {
          for (size_t i = 0; i < axes_.size(); ++i) {
            const size_t a = axes_[i];
            const ssize_t x = reference[a];
            target[a] = x;
          }
        }

        //! set position along relevant axes of \a targets to that of \a reference
        template <class Set, class Set2, class Set3>
        void set_position (const Set& reference, Set2& target, Set3& target2) const {
          for (size_t i = 0; i < axes_.size(); ++i) {
            const size_t a = axes_[i];
            const ssize_t x = reference[a];
            target[a] = x;
            target2[a] = x;
          }
        }

        //! set position along relevant axes of \a targets to that of \a reference
        template <class Set, class Set2, class Set3, class Set4>
        void set_position (const Set& reference, Set2& target, Set3& target2, Set4& target3) const {
          for (size_t i = 0; i < axes_.size(); ++i) {
            const size_t a = axes_[i];
            const ssize_t x = reference[a];
            target[a] = x;
            target2[a] = x;
            target3[a] = x;
          }
        }

        template <class Set>
        size_t max_axis (const Set& set) const {
          size_t a = 0;
          for (size_t i = 0; i < axes_.size(); ++i)
            if (axes_[i] > a) 
              a = axes_[i];
          return a;
        }

      private:
        const std::vector<size_t> axes_;
        bool cont_;
        ProgressBar progress_;

        template <class Set>
        void next_impl (size_t axis, Set& set) {
          size_t a = axes_[axis];
          if (set[a] + 1 < set.dim (a)) ++set[a];
          else {
            if (axis+1 == axes_.size()) {
              cont_ = false;
            }
            else {
              next_impl (axis+1, set);
              if (cont_) set[a] = 0;
            }
          }
        }

        template <class Set, class Set2>
        void next_impl (size_t axis, Set& set, Set2& set2) {
          size_t a = axes_[axis];
          if (set[a] + 1 < set.dim (a)) {
            ++set[a];
            ++set2[a];
          }
          else {
            if (axis+1 == axes_.size()) {
              cont_ = false;
            }
            else {
              next_impl (axis+1, set, set2);
              if (cont_) {
                set[a] = 0;
                set2[a] = 0;
              }
            }
          }
        }

        template <class Set, class Set2, class Set3>
        void next_impl (size_t axis, Set& set, Set2& set2, Set3& set3) {
          size_t a = axes_[axis];
          if (set[a] + 1 < set.dim (a)) {
            ++set[a];
            ++set2[a];
            ++set3[a];
          }
          else {
            if (axis+1 == axes_.size()) {
              cont_ = false;
            }
            else {
              next_impl (axis+1, set, set2, set3);
              if (cont_) {
                set[a] = 0;
                set2[a] = 0;
                set3[a] = 0;
              }
            }
          }
        }

    };


    //! @}
  }
}

#endif


