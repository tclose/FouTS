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


#ifndef __bts_math_munkres_h__
#define __bts_math_munkres_h__

namespace BTS {
  namespace Math {

    template <typename T> class Munkres;
  }
}

#include <set>


#include "math/matrix.h"
#include "math/vector.h"
#include "bts/common.h"

namespace BTS {

  namespace Math {

    class Index {

      public:

        size_t row, column;

        Index () : row(99999999), column(99999999) {}
        Index (size_t row, size_t column) : row(row), column(column) {}

        friend std::ostream& operator<< (std::ostream& stream, const Index& idx)
        {
          stream << "(" << idx.row << ", " << idx.column << ")";

          return stream;
        }

    };

    template <typename T> class Munkres {

      protected:

        const static int    UNASSIGNED;
        const static int    NOT_PROPOSED;

        size_t nrows, ncolumns;

        MR::Math::Matrix<T> Costs;


        std::set<size_t>      all_rows;
        std::set<size_t>      all_cols;

        std::set<size_t>      open_rows;
        std::set<size_t>      open_cols;

        std::vector<int>    assigned_col;
        std::vector<int>    assigned_row;

        std::vector<int>    proposed_col;

      public:


        Munkres(size_t nrows, size_t ncolumns) : nrows(nrows), ncolumns(ncolumns), Costs(nrows,ncolumns), assigned_col(nrows), assigned_row(ncolumns), proposed_col(nrows) {

          if (ncolumns < nrows)
            throw Exception ("Number of columns (" + str(ncolumns) + ") must exceed number of rows (" + str(nrows) + ").");

          for (size_t row_i = 0; row_i < nrows; row_i++) all_rows.insert(row_i);
          for (size_t col_i = 0; col_i < ncolumns; col_i++) all_cols.insert(col_i);

        }

        ~Munkres()                  {}


        void                match(MR::Math::Matrix<T>& Costs, std::vector<size_t>& matched_indices);

      protected:

        void                  subtract_row_min();

        T                     min_open_value();

        void                  initial_assignment();

        bool                  found_solution();

        Index                 new_match();

        void                  increment_costs();

        void                  add_to_closed_rows(T value);

        void                  subtract_from_open_columns(T value);

        void                  reassign_matches(Index new_match);

        void                  assign(size_t row, size_t col);

        void                  unassign(size_t row, size_t col);

        void                  propose(size_t row, size_t col);

        void                  clear_assignments();

        void                  clear_proposals();

        void                  open_all_rows();

        void                  open_all_columns();

      public:

        friend std::ostream& operator<< (std::ostream& stream, const Munkres<T>& H)
        {

          stream << std::endl << "----------------" << std::endl << "Cost: " << std::endl << H.Costs << std::endl;

          stream << "Assigned cols: ";
          for (std::vector<int>::const_iterator it = H.assigned_col.begin(); it != H.assigned_col.end(); it++)
            stream << *it << ", ";
          stream << std::endl;

          stream << "Assigned rows: ";
          for (std::vector<int>::const_iterator it = H.assigned_row.begin(); it != H.assigned_row.end(); it++)
            stream << *it << ", ";
          stream << std::endl;

          stream << "proposed cols: ";
          for (std::vector<int>::const_iterator it = H.proposed_col.begin(); it != H.proposed_col.end(); it++)
            stream << *it << ", ";
          stream << std::endl;

          stream << "Open rows: ";
          for (std::set<size_t>::const_iterator row_it = H.open_rows.begin(); row_it != H.open_rows.end(); row_it++)
            stream << *row_it << ", ";
          stream << std::endl;

          std::cout << "Open columns: ";
          for (std::set<size_t>::const_iterator col_it = H.open_cols.begin(); col_it != H.open_cols.end(); col_it++)
            stream << *col_it << ", ";
          stream << std::endl;

          stream << "-------end------" << std::endl << std::endl;

          return (stream);
        }


    };


    template <typename T> const int            Munkres<T>::UNASSIGNED = -1;
    template <typename T> const int            Munkres<T>::NOT_PROPOSED = -1;


    template <typename T> void                 Munkres<T>::match(MR::Math::Matrix<T>& Costs, std::vector<size_t>& matched_indices) {


      if ((Costs.rows() != nrows) || (Costs.columns() != ncolumns) || (matched_indices.size() != nrows))
        throw Exception ("Size of cost matrix " + str(Costs.rows()) + "x" + str(Costs.columns()) + " or matched_indices vector " + str(matched_indices.size()) + " does not match Munkres object " + str(nrows) + "x" + str(ncolumns) + ".");

      this->Costs = Costs;

      open_rows = all_rows;
      open_cols = all_cols;

      clear_assignments();

      subtract_row_min();

      initial_assignment();

      while (!found_solution()) {

        clear_proposals();

        Index proposed_match = new_match();

        reassign_matches(proposed_match);

      }

      //Assign assigned cols to output vector reference (and cast to unsigned ints)
      for (size_t i = 0; i < nrows; i++)
        matched_indices[i] = (size_t)assigned_col[i];

    }


    template <typename T> void                  Munkres<T>::assign(size_t row, size_t col) {


      assigned_row[col]  = row;
      assigned_col[row]  = col;

#ifdef _DEBUG_VERBOSE_
      std::cout << "assigned: " << row << ", " << col << *this << std::endl;
#endif

    }


    template <typename T> void                  Munkres<T>::unassign(size_t row, size_t col) {



      assigned_row[col]  = UNASSIGNED;
      assigned_col[row]  = UNASSIGNED;

#ifdef _DEBUG_VERBOSE_
      std::cout << "unassigned: " << row << ", " << col << *this << std::endl;
#endif

    }


    template <typename T> void                  Munkres<T>::propose(size_t row, size_t col) {


      proposed_col[row] = col;

#ifdef _DEBUG_VERBOSE_
      std::cout << "proposed: " << row << ", " << col << *this << std::endl;
#endif

    }


    template <typename T> void                  Munkres<T>::open_all_rows() {

      open_rows = all_rows;

#ifdef _DEBUG_VERBOSE_
      std::cout << "open all rows: " << *this << std::endl;
#endif

    }


    template <typename T> void                  Munkres<T>::open_all_columns() {

      open_cols = all_cols;

#ifdef _DEBUG_VERBOSE_
      std::cout << "open all cols: " << *this << std::endl;
#endif

    }


    template <typename T> void                  Munkres<T>::clear_assignments() {


      fill(assigned_col.begin(), assigned_col.end(), UNASSIGNED);
      fill(assigned_row.begin(), assigned_row.end(), UNASSIGNED);

#ifdef _DEBUG_VERBOSE_
      std::cout << "cleared assignments" << *this << std::endl;
#endif

    }


    template <typename T> void                  Munkres<T>::clear_proposals() {

      fill(proposed_col.begin(), proposed_col.end(), NOT_PROPOSED);

#ifdef _DEBUG_VERBOSE_
      std::cout << "cleared proposals" << *this << std::endl;
#endif
    }


    template <typename T> void                  Munkres<T>::subtract_row_min() {


      for (size_t row_i = 0; row_i < nrows; row_i++) {
        MR::Math::Vector<T> row = Costs.row(row_i);
        size_t index;
        T min_val = min (row, index);
        row -= min_val;
      }

#ifdef _DEBUG_VERBOSE_
      std::cout << "subtracted min cost in row" << *this << std::endl;
#endif


    }


    template <typename T> void                  Munkres<T>::initial_assignment() {

      std::set<size_t>::iterator next_row = open_rows.begin();

      for (std::set<size_t>::iterator row_it = next_row; row_it != open_rows.end(); row_it = next_row) {

        ++next_row;

        for (std::set<size_t>::iterator col_it = open_cols.begin(); col_it != open_cols.end(); col_it++) {

          if ( Costs (*row_it, *col_it) == 0 ) {

            assign(*row_it, *col_it);
            open_rows.erase(row_it);
            open_cols.erase(col_it);
            break;
          }
        }

      }

      open_all_rows();
      open_all_columns();

#ifdef _DEBUG_VERBOSE_
      std::cout << "Initial assignment: " <<  *this << std::endl;
#endif


    }



// Close each column containing a starred zero.  If 'n' columns are closed, the starred zeros describe a complete set of unique
// assignments.

    template <typename T> bool                  Munkres<T>::found_solution() {

      size_t num_assigned = 0;

      for (size_t col_i = 0; col_i < ncolumns; col_i++) {

        if (assigned_row[col_i] != UNASSIGNED) {
          open_cols.erase(col_i);
          num_assigned++;
        }

      }

      bool found_sol = false;
      if (num_assigned == nrows)
        found_sol = true;

#ifdef _DEBUG_VERBOSE_
      std::cout << "Covered columns: " << *this << std::endl;
#endif

      return found_sol;

    }



// Find a open zero and propose it.  If there is no assigned zero in the row containing this proposed zero, reassign the new chain
// of matches. Otherwise, close this row and open the column containing the assigned zero. Continue in this manner until there are
// no open zeros left.

    template <typename T> Index                 Munkres<T>::new_match() {

      while (true) {

        std::set<size_t>::iterator next_row = open_rows.begin();

        for (std::set<size_t>::iterator row_it = next_row; row_it != open_rows.end(); row_it = next_row) {

          ++next_row;

          for (std::set<size_t>::iterator col_it = open_cols.begin(); col_it != open_cols.end(); col_it++) {

            if ( Costs (*row_it, *col_it) == 0) {

              propose(*row_it, *col_it);

              if ( assigned_col[*row_it] == UNASSIGNED ) {

#ifdef _DEBUG_VERBOSE_
                std::cout << "proposed new zero, (no assignment returning...) (" << *row_it << ", " << *col_it << ")"  << *this << std::endl;
#endif
                return Index(*row_it, *col_it);
              } else {
                open_cols.insert(assigned_col[*row_it]);
                open_rows.erase(row_it);

#ifdef _DEBUG_VERBOSE_
                std::cout << "proposed new zero, (found assignment continuing...)" << *this << std::endl;
#endif

                break;
              }

            }

          }

        }

        increment_costs();

      }

    }



// Add the minimum open value to every element of each closed row, and subtract it from every element of each opened column.
// This step uses the smallest opened value in the cost matrix to modify the matrix. It may seem that since the values in the
// cost matrix are being altered, we would lose sight of the original problem.  However, we are only changing certain values
/// that have already been tested and found not to be elements of the minimal assignment.  Also we are only
// changing the values by an amount equal to the smallest value in the cost matrix, so we will not jump over the optimal
// (i.e. minimal assignment) with this change.


    template <typename T> void                  Munkres<T>::increment_costs() {

        T min_value = min_open_value();

        add_to_closed_rows(min_value);

        subtract_from_open_columns(min_value);


    }


    template <typename T> T                     Munkres<T>::min_open_value() {


      T min_value = INFINITY;

      for (std::set<size_t>::iterator row_it = open_rows.begin(); row_it != open_rows.end(); row_it++) {

        for (std::set<size_t>::iterator col_it = open_cols.begin(); col_it != open_cols.end(); col_it++) {

          if ( Costs(*row_it, *col_it) < min_value )
            min_value = Costs(*row_it, *col_it);

        }

      }

#ifdef _DEBUG_VERBOSE_
      std::cout << "min open value: " << min_value << *this << std::endl;
#endif

      return min_value;

    }


    template <typename T> void                  Munkres<T>::add_to_closed_rows(T value) {


      for (size_t row_i = 0; row_i < nrows; row_i++) {

        //If row is closed.
        if (open_rows.find(row_i) == open_rows.end()) {

          for (size_t col_i = 0; col_i < ncolumns; col_i++)
            Costs (row_i, col_i) += value;

        }

      }

    }


    template <typename T> void                  Munkres<T>::subtract_from_open_columns(T value) {

      for (size_t row_i = 0; row_i < nrows; row_i++) {

        for (std::set<size_t>::iterator col_it = open_cols.begin(); col_it != open_cols.end(); col_it++) {

          Costs( row_i, *col_it) -= value;

        }
      }
    }



// Construct a chain of alternating proposed and assigned zeros as follows.  Let Z0 represent the unclosed proposed zero
// found in Step 4.  Let Z1 denote the assigned zero in the column of Z0 (if any). Let Z2 denote the proposed zero in the
// row of Z1 (there will always be one).  Continue until the series terminates at a proposed zero that has no assigned zero
// in its column.  Unassign each assigned zero of the series, assign each proposed zero of the series, remove all proposals
// and open all rows and columns.  This step is very similar to the augmenting location algorithm (for solving the maximal
// matching problem)

    template <typename T> void                  Munkres<T>::reassign_matches(Index new_match) {

      std::vector<Index> unassign_indices;
      std::vector<Index> assign_indices;

      assign_indices.push_back(new_match);

      size_t col_i = new_match.column;

      while (assigned_row[col_i] != UNASSIGNED) {


        Index unassign_idx (assigned_row[col_i], col_i);
        Index assign_idx (assigned_row[col_i], proposed_col[assigned_row[col_i]]);

        assign_indices.push_back(assign_idx);
        unassign_indices.push_back(unassign_idx);

        col_i = proposed_col[assigned_row[col_i]];

      }

#ifdef _DEBUG_VERBOSE_
      std::cout << "reassigning matches" << *this << std::endl;
#endif


      for (std::vector<Index>::iterator it = unassign_indices.begin(); it != unassign_indices.end(); it++)
        unassign(it->row, it->column);

      for (std::vector<Index>::iterator it = assign_indices.begin(); it != assign_indices.end(); it++)
        assign(it->row, it->column);

      assign_indices.clear();
      unassign_indices.clear();

      open_all_rows();
      open_all_columns();


    }




  }

}


#endif

