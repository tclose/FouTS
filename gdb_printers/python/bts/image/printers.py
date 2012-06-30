# Pretty-printers for 'BTS::Image::*'.

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


class ImageVoxelPrinter:
    "Print Image::*::Voxel"

    class _iterator:
        def __init__ (self, start, finish):
            self.item = start
            self.finish = finish
            self.count = 0

        def __iter__(self):
            return self

        def next(self):
            if self.item == self.finish:
                raise StopIteration
            count = self.count
            self.count = self.count + 1
            elt = self.item.dereference()
            self.item = self.item + 1
            return ('[%d]' % count, elt)

    def __init__(self, typename, template_type, val):
        self.typename = typename
        self.template_type = template_type
        self.base_val = val.cast(gdb.lookup_type('BTS::Image::Voxel<%s>' % template_type))


    def children(self):
        return self._iterator(self.base_val['intensities']['_M_impl']['_M_start'],
                              self.base_val['intensities']['_M_impl']['_M_finish'])

    def to_string(self):
        return self.typename

    def display_hint(self):
        if self.template_type == 'double':
        	  return 'array'
        else:
            return ''
          
          
          
class ImageDensityVoxelPrinter:
    "Print Image::Density::Voxel"

    class _iterator:
        def __init__ (self, val):
            self.field = 'coordinate'
            self.val = val

        def __iter__(self):
            return self

        def next(self):
            if self.field == 'end':
                raise StopIteration
            elif self.field == 'coordinate':
                self.field = 'centre_point'
                return ('coordinate', self.val['coordinate'])
            elif self.field == 'centre_point':
                self.field = 'tangent'
                return ('centre_point', self.val['centre_point'])
            elif self.field == 'dens':
                self.field = 'end'
                return ('dens', self.val['dens'])


    def __init__(self, val):
        self.val = val

    def to_string(self):
        return 'Image::Density::Voxel'
          
          

class VoxelIterator:
    def __init__(self, rbtree):
        self.size = rbtree['_M_t']['_M_impl']['_M_node_count']
        self.node = rbtree['_M_t']['_M_impl']['_M_header']['_M_left']
        self.count = 0

    def __iter__(self):
        return self

    def __len__(self):
        return int (self.size)

    def next(self):
        if self.count == self.size:
            raise StopIteration
        result = self.node
        self.count = self.count + 1
        if self.count < self.size:
            # Compute the next node.
            node = self.node
            if node.dereference()['_M_right']:
                node = node.dereference()['_M_right']
                while node.dereference()['_M_left']:
                    node = node.dereference()['_M_left']
            else:
                parent = node.dereference()['_M_parent']
                while node == parent.dereference()['_M_right']:
                    node = parent
                    parent = parent.dereference()['_M_parent']
                if node.dereference()['_M_right'] != parent:
                    node = parent
            self.node = node
        return result



class ImageBufferPrinter:
    "Print Image::*::Buffer"

    # Turn an VoxelIterator into a pretty-print iterator.
    class _iter:
        def __init__(self, rbiter, type):
            self.rbiter = rbiter
            self.count = 0
            self.type = type

        def __iter__(self):
            return self

        def next(self):
            if self.count % 2 == 0:
                n = self.rbiter.next()
                n = n.cast(self.type).dereference()['_M_value_field']
                self.pair = n
                item = n['first']
            else:
                item = self.pair['second']
            result = ('[%d]' % self.count, item)
            self.count = self.count + 1
            return result

    def __init__ (self, typename, val):
        
        self.val = val
        
        if typename == 'BTS::Image::Buffer':
            self.typename = 'BTS::Image::Buffer'
            self.voxel_typename = 'BTS::Image::Voxel<double>'
        if typename == 'BTS::Image::Density::Buffer':
            self.typename = 'BTS::Image::Density::Buffer'
            self.voxel_typename = 'BTS::Image::Density::Voxel'
        elif typename == 'BTS::Image::Container::Buffer':
            self.typename = 'BTS::Image::Container::Buffer< %s >' % val.type.template_argument(0)
            self.voxel_typename = 'BTS::Image::Container::Voxel< %s >' % val.type.template_argument(0)
        elif typename == 'BTS::Image::Expected::Buffer_tpl':
            self.typename = 'BTS::Image::Expected::Buffer_tpl< %s >' % val.type.template_argument(0)
            self.voxel_typename = '%s' % val.type.template_argument(0)            
        elif typename == 'BTS::Image::Observed::Buffer_tpl':
            self.typename = 'BTS::Image::Observed::Buffer_tpl< %s >' % val.type.template_argument(0)
            self.voxel_typename = '%s' % val.type.template_argument(0)     
        else:
            self.typename = typename
            self.voxel_typename = '%s::Voxel' % '::'.join(typename.split('::')[:-1])

          


    def to_string (self):
        return '%s' % (self.typename)

    def children (self):
        nodetype = gdb.lookup_type('std::_Rb_tree_node< std::pair< const BTS::Image::Coord, %s > >' % self.voxel_typename)
        nodetype = nodetype.pointer()
        return self._iter (VoxelIterator (self.val['voxels']), nodetype)

    def display_hint (self):
        return 'map'

class ImageReferencePrinter:
    "Print Image::Reference"

    # Turn an VoxelIterator into a pretty-print iterator.
    class _iter:
        def __init__(self, rbiter, type):
            self.rbiter = rbiter
            self.count = 0
            self.type = type

        def __iter__(self):
            return self

        def next(self):
            if self.count % 2 == 0:
                n = self.rbiter.next()
                n = n.cast(self.type).dereference()['_M_value_field']
                self.pair = n
                item = n['first']
            else:
                item = self.pair['second']
            result = ('[%d]' % self.count, item)
            self.count = self.count + 1
            return result

    def __init__ (self, typename, val):
        
        self.val = val
        self.typename = typename
        self.voxel_typename = 'std::__debug::vector<%s*, std::allocator<%s*> >' % (val.type.template_argument(0), val.type.template_argument(0))

    def to_string (self):
        return '%s' % (self.typename)

    def children (self):
        nodetype = gdb.lookup_type('std::_Rb_tree_node< std::pair< const BTS::Image::Coord, %s > >' % self.voxel_typename)
        nodetype = nodetype.pointer()
        return self._iter (VoxelIterator (self.val['voxels']), nodetype)

    def display_hint (self):
        return 'map'


class ImagePrinter:
    "Print Image::Observed"

    def __init__(self, val):
        self.val = val

    def to_string(self):
        if self.val['type'].string() == 'trilinear':
            return self.val.cast(gdb.lookup_type('BTS::Image::Expected::Trilinear::Buffer'))

        elif self.val['type'].string() == 'top_hat':
            return self.val.cast(gdb.lookup_type('BTS::Image::Expected::TopHat::Buffer'))
            
        elif self.val['type'].string() == 'gaussian':
            return self.val.cast(gdb.lookup_type('BTS::Image::Expected::Gaussian::Buffer'))
    
        elif self.val['type'].string() == 'quartic':
            return self.val.cast(gdb.lookup_type('BTS::Image::Expected::Quartic::Buffer'))
    
        elif self.val['type'].string() == 'realistic':
            return self.val.cast(gdb.lookup_type('BTS::Image::Expected::Realistic::Buffer'))
                        
        elif self.val['type'].string() == 'sinc':
            return self.val.cast(gdb.lookup_type('BTS::Image::Expected::Sinc::Buffer'))
   
        elif self.val['type'].string() == 'reverse_sqrt':
            return self.val.cast(gdb.lookup_type('BTS::Image::Expected::ReverseSqrt::Buffer'))            
                    
        else:
            return ("Unknown observed/expected image type '%s'." % self.val['type'].string())


def register_bts_image_printers (obj):
    "Register Bayesian Tractlet Sampling (BTS) Image pretty-printers with objfile Obj."

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

def build_bts_image_dictionary ():
    # bts objects requiring pretty-printing.
    # In order from:
    pretty_printers_dict[re.compile('^BTS::Image::Buffer$')] = lambda val: ImageBufferPrinter("BTS::Image::Buffer", val)        
    pretty_printers_dict[re.compile('^BTS::Image::Expected::Quartic::Buffer$')] = lambda val: ImageBufferPrinter("BTS::Image::Expected::Quartic::Buffer", val)    
    pretty_printers_dict[re.compile('^BTS::Image::Expected::Gaussian::Buffer$')] = lambda val: ImageBufferPrinter("BTS::Image::Expected::Gaussian::Buffer", val)    
    pretty_printers_dict[re.compile('^BTS::Image::Expected::Trilinear::Buffer$')] = lambda val: ImageBufferPrinter("BTS::Image::Expected::Trilinear::Buffer", val)    
    pretty_printers_dict[re.compile('^BTS::Image::Observed::Buffer$')] = lambda val: ImageBufferPrinter("BTS::Image::Observed::Buffer", val)    
    pretty_printers_dict[re.compile('^BTS::Image::Density::Buffer$')] = lambda val: ImageBufferPrinter("BTS::Image::Density::Buffer", val)    
    pretty_printers_dict[re.compile('^BTS::Image::Container::Buffer<.*>$')] = lambda val: ImageBufferPrinter("BTS::Image::Container::Buffer", val)     
    pretty_printers_dict[re.compile('^BTS::Image::Reference::Buffer<.*>$')] = lambda val: ImageReferencePrinter("BTS::Image::Reference::Buffer", val)     
    pretty_printers_dict[re.compile('^BTS::Image::Expected::Buffer$')] = lambda val: ImagePrinter(val)  
    pretty_printers_dict[re.compile('^BTS::Image::Expected::Buffer_tpl<.*>$')] = lambda val: ImageBufferPrinter('BTS::Image::Expected::Buffer_tpl', val)
    pretty_printers_dict[re.compile('^BTS::Image::Observed::Buffer_tpl<.*>$')] = lambda val: ImageBufferPrinter('BTS::Image::Observed::Buffer_tpl', val)
    pretty_printers_dict[re.compile('^BTS::Image::Voxel<double>$')] = lambda val: ImageVoxelPrinter("BTS::Image::Voxel<double>", 'double', val)
    pretty_printers_dict[re.compile('^BTS::Image::Expected::Quartic::Voxel$')] = lambda val: ImageVoxelPrinter("BTS::Image::Expected::Quartic::Voxel", 'double', val)
    pretty_printers_dict[re.compile('^BTS::Image::Expected::Gaussian::Voxel$')] = lambda val: ImageVoxelPrinter("BTS::Image::Expected::Gaussian::Voxel", 'double', val)
    pretty_printers_dict[re.compile('^BTS::Image::Expected::Trilinear::Voxel$')] = lambda val: ImageVoxelPrinter("BTS::Image::Expected::Trilinear::Voxel", 'double', val)        
    pretty_printers_dict[re.compile('^BTS::Image::Expected::TopHat::Voxel$')] = lambda val: ImageVoxelPrinter("BTS::Image::Expected::TopHat::Voxel",'double',  val)        
    pretty_printers_dict[re.compile('^BTS::Image::Observed::Voxel$')] = lambda val: ImageVoxelPrinter("BTS::Image::Observed::Voxel", 'double', val)
    pretty_printers_dict[re.compile('^BTS::Image::Container::Voxel<BTS::Fibre::Strand>$')] = lambda val: ImageVoxelPrinter("BTS::Image::Container::Voxel<BTS::Fibre::Strand>", 'BTS::Fibre::Strand', val)
    pretty_printers_dict[re.compile('^BTS::Image::Container::Voxel<BTS::Fibre::Tractlet>$')] = lambda val: ImageVoxelPrinter("BTS::Image::Container::Voxel<BTS::Fibre::Tractlet>", 'BTS::Fibre::Tractlet', val)
    pretty_printers_dict[re.compile('^BTS::Image::Container::Voxel<BTS::Fibre::Strand::Tensor>$')] = lambda val: ImageVoxelPrinter("BTS::Image::Container::Voxel<BTS::Fibre::Strand::Tensor>", 'BTS::Fibre::Strand::Tensor', val)
    pretty_printers_dict[re.compile('^BTS::Image::Container::Voxel<BTS::Fibre::Tractlet::Tensor>$')] = lambda val: ImageVoxelPrinter("BTS::Image::Container::Voxel<BTS::Fibre::Tractlet::Tensor>", 'BTS::Fibre::Tractlet::Tensor', val)        
    pretty_printers_dict[re.compile('^BTS::Image::Container::Voxel< MR::Math::Vector<double> >$')] = lambda val: ImageVoxelPrinter("BTS::Image::Container::Voxel<MR::Math::Vector<double>>", ' MR::Math::Vector<double> ', val)                    
    pretty_printers_dict[re.compile('^BTS::Image::Container::Voxel< MR::Math::Matrix<double> >$')] = lambda val: ImageVoxelPrinter("BTS::Image::Container::Voxel<MR::Math::Matrix<double>>", ' MR::Math::Matrix<double> ', val)                    
    pretty_printers_dict[re.compile('^BTS::Image::Density::Voxel$')] = lambda val: ImageDensityVoxelPrinter(val)
    
pretty_printers_dict = {}

build_bts_image_dictionary ()