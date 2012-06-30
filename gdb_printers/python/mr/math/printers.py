# Pretty-printers for 'BTS::Fibre::*'.

# Copyright (C) 2008, 2009 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import gdb
import itertools
import re


class MathVectorPrinter:
    "Print Math::Vector"

    class _iterator:
        def __init__ (self, size, stride, data_ptr):
            self.size = size
            self.stride = stride
            self.data_ptr = data_ptr
            self.count = 0

        def __iter__(self):
            return self

        def next(self):
            if self.count == self.size:
                raise StopIteration
            count = self.count
            elt = self.data_ptr.dereference()
            self.data_ptr = self.data_ptr + self.stride
            self.count = self.count + 1
            return ('[%d]' % count, elt)

    def __init__(self, typename, val):
        self.typename = typename
        self.val = val

    def children(self):
        return self._iterator(self.val['size'],
                              self.val['stride'],
                              self.val['data'])

    def to_string(self):
        return ("%s: size %d" % (self.typename, self.val['size']))


    def display_hint(self):
        return 'array'

class RowPrinter:
  
  def __init__(self):
    print 'done'
    

class MathMatrixPrinter:
    "Print Math::Matrix"

    class _iterator:
        def __init__ (self, num_rows, num_cols, row_size, data_ptr):
            self.num_rows = num_rows
            self.num_cols = num_cols
            self.row_size = row_size
            self.data_ptr = data_ptr
            self.row_count = 0
            self.col_count = 0

        def __iter__(self):
            return self

        def next(self):
            if self.row_count == self.num_rows:
                raise StopIteration
              
            row_count = self.row_count
            col_count = self.col_count
            
            if col_count == 0:
                print '\n'
            
            elt = self.data_ptr.dereference()
            
            self.col_count += 1
            if self.col_count == self.num_cols:
                self.row_count += 1
                self.col_count = 0
                self.data_ptr = self.data_ptr + (self.row_size - self.num_cols + 1)
            else:
                self.data_ptr = self.data_ptr + 1
              
            return ('(%d,%d)' % (row_count, col_count), elt)

    def __init__(self, typename, val):
        self.typename = typename
        self.val = val

    def children(self):
        return self._iterator(self.val['size1'],
                              self.val['size2'],
                              self.val['tda'],
                              self.val['data'])

    def to_string(self):
        return ("%s: rows %d, cols %d" % (self.typename, self.val['size1'], self.val['size2']))


    def display_hint(self):
        return 'array'


def register_mr_math_printers (obj):
    "Register MRtrix (MR) Math pretty-printers with objfile Obj."

    if obj == None:
        obj = gdb

    obj.pretty_printers.append (lookup_function)

def lookup_function (val):
    "Look-up and return a pretty-printer that can print val."

    # Get the type.
    type = val.type

    # If it points to a reference, get the reference.
    if type.code == gdb.TYPE_CODE_REF:
        type = type.target ()

    # Get the unqualified type, stripped of typedefs.
    type = type.unqualified ().strip_typedefs ()

    # Get the type name.    
    typename = type.tag
    if typename == None:
        return None

    # Iterate over local dictionary of types to determine
    # if a printer is registered for that type.  Return an
    # instantiation of the printer if found.
    for function in pretty_printers_dict:
        if function.search (typename):
            return pretty_printers_dict[function] (val)
        
    # Cannot find a pretty printer.  Return None.
    return None

def build_mr_math_dictionary ():
    # bts objects requiring pretty-printing.
    # In order from:
    pretty_printers_dict[re.compile('^MR::Math::Vector<.*>$')] = lambda val: MathVectorPrinter('MR::Math::Vector', val)
    pretty_printers_dict[re.compile('^MR::Math::Matrix<.*>$')] = lambda val: MathMatrixPrinter('MR::Math::Matrix', val)
    pretty_printers_dict[re.compile('^BTS::Fibre::Strand::Set::Tensor$')] = lambda val: MathMatrixPrinter('BTS::Fibre::Strand::Set::Tensor', val)
    pretty_printers_dict[re.compile('^BTS::Fibre::Tract::Set::Tensor$')] = lambda val: MathMatrixPrinter('BTS::Fibre::Tract::Set::Tensor', val)  
    pretty_printers_dict[re.compile('^BTS::MCMC::State$')] = lambda val: MathVectorPrinter('BTS::MCMC::State', val)
    pretty_printers_dict[re.compile('^BTS::MCMC::State::Tensor$')] = lambda val: MathMatrixPrinter('BTS::MCMC::State::Tensor', val) 
            
pretty_printers_dict = {}

build_mr_math_dictionary ()
