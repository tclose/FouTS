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

class TriplePrinter:
    "Print Triple<T>"

    def __init__(self, val):
        self.val = val

    def to_string(self):
        return self.val['p']

    def display_hint(self):
        return 'array'


class CoordPrinter:
    "Print Coord"

    def __init__(self, val):
      self.val = val

    def to_string(self):
      x = self.val['data'].dereference()
      y = (self.val['data'] + self.val['stride']).dereference()
      z = (self.val['data'] + 2 * self.val['stride']).dereference()
      return "[%g, %g, %g]" % (x, y, z)

    def display_hint(self):
        return 'array'


class FibreStrandTrackPrinter:
    "Print Fibre::Strand and Fibre::Track objects"

    def __init__(self, typename, val):
        self.typename = typename
        self.val = val
        self.stride = val['stride']
        self.data_ptr = val['data']
        self.size = self.val['sze']

    def children(self):
      for i in xrange(self.size):
        x = self.next_elem()
        y = self.next_elem()
        z = self.next_elem()
        yield '%d' % i, '[%g, %g, %g]' % (x, y, z)

    def to_string(self):
        if self.typename == 'Fibre::Strand':
            size_name = 'degree'
        else:
            size_name = 'size'
        return ("%s (%s=%d):" % (self.typename, size_name, self.size))

    def next_elem(self):
      elem = self.data_ptr.dereference()
      self.data_ptr = self.data_ptr + self.stride
      return elem


class FibreTractletPrinter:
    "Print Fibre::Strand and Fibre::Track objects"

    def __init__(self, val):
        self.val = val
        self.stride = val['stride']
        self.data_ptr = val['data']
        self.size = self.val['sze']
        self.degree = self.val['dgree']
        self.vsize = self.val['size']
        self.bsize = self.degree * 9
        self.num_props = self.vsize - self.bsize


    def children(self):
      for i in xrange(self.degree):
        yield '%d' % i, '[%g, %g, %g | %g, %g, %g | %g, %g, %g ]' % (self.get_elem(i, 0, 0), self.get_elem(i, 0, 1), self.get_elem(i, 0, 2),
                                                                     self.get_elem(i, 1, 0), self.get_elem(i, 1, 1), self.get_elem(i, 1, 2),
                                                                     self.get_elem(i, 2, 0), self.get_elem(i, 2, 1), self.get_elem(i, 2, 2))

    def to_string(self):
        alpha = 1
        if self.num_props:
            alpha = self.get_prop(0)
        return ("Fibre::Tractlet (degree={}, alpha={}): ".format(self.degree, alpha))

    def get_elem(self, degree_i, ax_i, dim_i):
      address = self.data_ptr + ((ax_i * self.degree + degree_i) * 3 + dim_i) * self.stride
      return address.dereference()

    def get_prop(self, prop_index):
        assert prop_index < self.num_props
        address = self.data_ptr + (self.bsize + prop_index) * self.stride
        return address.dereference()
class FibreSectionPrinter:
    "Print Fibre::Strand::Section"

    class _iterator:
        def __init__ (self, section, stride, data):
            self.field = 'position'
            self.stride = stride
            self.data = data

        def __iter__(self):
            return self

        def next(self):
            if self.field == 'end':
                raise StopIteration
            elif self.field == 'position':
                self.field = 'tangent'
                x = str(self.data_ptr.dereference())
                y = str((self.data_ptr + self.stride * 1).dereference())
                z = str((self.data_ptr + self.stride * 2).dereference())
                return ('position', gdb.Value("[" + x + ", " + y + ", " + z + "]"))
            elif self.field == 'tangent':
                self.field = 'intensity'
                x = str((self.data_ptr + self.stride * 3).dereference())
                y = str((self.data_ptr + self.stride * 4).dereference())
                z = str((self.data_ptr + self.stride * 5).dereference())
                return ('tangent', gdb.Value("[" + x + ", " + y + ", " + z + "]"))
            elif self.field == 'intensity':
                self.field = 'end'
                intensity = (self.data + 6 * self.stride).dereference()
                return ('intensity', gdb.Value(intensity))

            else:
                raise StopIteration


    def __init__(self, typename, val):
        self.typename = typename
        self.val = val

    def children(self):
        return self._iterator(self.val['stride'],
                              self.val['data'])

    def to_string(self):
        return "%s" % self.typename


    class _iterator:
        def __init__ (self, size, degree, stride, data_ptr):
            self.size = size
            self.degree = degree
            self.stride = stride
            self.data_ptr = data_ptr
            self.count = 0
            self.x1_offset = 0
            self.y1_offset = stride
            self.z1_offset = stride * 2
            self.x2_offset = stride * degree * 3
            self.y2_offset = stride * (degree * 3 + 1)
            self.z2_offset = stride * (degree * 3 + 2)
            self.x3_offset = stride * (degree * 6)
            self.y3_offset = stride * (degree * 6 + 1)
            self.z3_offset = stride * (degree * 6 + 2)

        def __iter__(self):
            return self

        def next(self):
            if self.count == self.size:
                raise StopIteration
            count = self.count
            x1 = (self.data_ptr + self.x1_offset).dereference()
            y1 = (self.data_ptr + self.y1_offset).dereference()
            z1 = (self.data_ptr + self.z1_offset).dereference()
            x2 = (self.data_ptr + self.x2_offset).dereference()
            y2 = (self.data_ptr + self.y2_offset).dereference()
            z2 = (self.data_ptr + self.z2_offset).dereference()
            x3 = (self.data_ptr + self.x3_offset).dereference()
            y3 = (self.data_ptr + self.y3_offset).dereference()
            z3 = (self.data_ptr + self.z3_offset).dereference()
            self.data_ptr = self.data_ptr + self.stride * 3
            self.count = self.count + 1
            return ('[%d]' % count, gdb.Value("%f.4 %f.4 %f.4 | %f.4 %f.4 %f.4 | %f.4 %f.4 %f.4 " % (x1, y1, z1, x2, y2, z2, x3, y3, z3)))

    def __init__(self, typename, val):
        self.typename = typename
        self.val = val



    def to_string(self):
        return ("%s: size %d" % (self.typename, self.val['sze']))


class FibreSetPrinter:
    "Print Fibre::Set"

    class _iterator:
        def __init__ (self, set, size):
          self.set = set
          self.size = size
          self.count = 0

        def __iter__(self):
            return self

        def next(self):
            if self.count == self.size:
                raise StopIteration
            count = self.count
            self.count = self.count + 1
            elem = self.set['elem'](count)
            return ('[%d]' % count, elem)

    def __init__(self, typename, val):
        self.typename = typename
        self.val = val

    def children(self):
        return self._iterator(self.val,
                              self.val['sze'])

    def to_string(self):
        return ("%s: size %d" % (self.typename, self.val['sze']))



def register_bts_fibre_printers (obj):
    "Register Bayesian Tractlet Sampling (BTS) Fibre pretty-printers with objfile Obj."

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

def build_bts_fibre_dictionary ():
    # bts objects requiring pretty-printing.
    # In order from:
    pretty_printers_dict[re.compile('^CoordView$')] = lambda val: CoordViewPrinter(val)
    pretty_printers_dict[re.compile('^BTS::Triple<.*>$')] = lambda val: TriplePrinter(val)
    pretty_printers_dict[re.compile('^BTS::Coord$')] = lambda val: CoordPrinter(val)
    pretty_printers_dict[re.compile('^BTS::Fibre::Strand$')] = lambda val: FibreStrandTrackPrinter("Fibre::Strand", val)
    pretty_printers_dict[re.compile('^BTS::Fibre::Tractlet$')] = lambda val: FibreTractletPrinter(val)
    pretty_printers_dict[re.compile('^BTS::Fibre::Track$')] = lambda val: FibreStrandTrackPrinter("Fibre::Track", val)
#    pretty_printers_dict[re.compile('^BTS::Fibre::Strand::Set$')] = lambda val: FibreSetPrinter("Fibre::Strand::Set", val)
#    pretty_printers_dict[re.compile('^BTS::Fibre::Tractlet::Set$')] = lambda val: FibreSetPrinter("Fibre::Tractlet::Set", val)
#    pretty_printers_dict[re.compile('^BTS::Fibre::Set<.*>$')] = lambda val: FibreSetPrinter("Fibre::Set", val)
#    pretty_printers_dict[re.compile('^BTS::Fibre::Track::Set$')] = lambda val: FibreSetPrinter("Fibre::Track::Set", val)
#    pretty_printers_dict[re.compile('^BTS::Fibre::Strand::Tensor$')] = lambda val: FibreStrandPrinter("Fibre::Strand::Tensor", val)
#    pretty_printers_dict[re.compile('^BTS::Fibre::Tractlet::Tensor$')] = lambda val: FibreTractletPrinter("Fibre::Tractlet::Tensor", val)
#    pretty_printers_dict[re.compile('^BTS::Fibre::Track::Tensor$')] = lambda val: FibreTrackPrinter("Fibre::Track::Tensor", val)
#    pretty_printers_dict[re.compile('^BTS::Fibre::Strand::Set::Tensor$')] = lambda val: FibreSetPrinter("Fibre::Strand::Set::Tensor", val)
#    pretty_printers_dict[re.compile('^BTS::Fibre::Tractlet::Set::Tensor$')] = lambda val: FibreSetPrinter("Fibre::Tractlet::Set::Tensor", val)
    pretty_printers_dict[re.compile('^BTS::Fibre::Tractlet::Section$')] = lambda val: FibreSectionPrinter("Fibre::Tractlet::Section", val)
    pretty_printers_dict[re.compile('^BTS::Fibre::Strand::BasicSection$')] = lambda val: FibreSectionPrinter("Fibre::Strand::BasicSection", val)
    pretty_printers_dict[re.compile('^BTS::Fibre::Strand::Section$')] = lambda val: FibreSectionPrinter("Fibre::Strand::Section", val)
#    pretty_printers_dict[re.compile('^BTS::Fibre::Tractlet::Section::Tensor$')] = lambda val: FibreTractletSectionPrinter("Fibre::Tractlet::Section::Tensor", val)
#    pretty_printers_dict[re.compile('^BTS::Fibre::Strand::Section::Tensor$')] = lambda val: FibreStrandSectionPrinter("Fibre::Strand::Section::Tensor", val)
#    pretty_printers_dict[re.compile('^BTS::Fibre::Strand::BasicSection::Tensor$')] = lambda val: FibreStrandSectionPrinter("Fibre::Strand::BasicSection::Tensor", val)


pretty_printers_dict = {}

build_bts_fibre_dictionary ()
