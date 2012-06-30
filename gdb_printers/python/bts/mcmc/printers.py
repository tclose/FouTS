# Pretty-printers for 'BTS::MCMC::*'.

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

class ProposalMomentumPrinter:
    "Print MCMC::Proposal::Momentum"

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
        
        if typename == 'MCMC::Proposal::Momentum':
            self.val = val
        else:
            self.val = val.cast(gdb.lookup_type('BTS::MCMC::Proposal::Momentum'))

    def children(self):
        return self._iterator(self.val['momen']['size'],
                              self.val['momen']['stride'],
                              self.val['momen']['data'])

    def to_string(self):
        return ("%s: size %d" % (self.typename, self.val['momen']['size']))


    def display_hint(self):
        return 'array'


def register_bts_mcmc_printers (obj):
    "Register Bayesian Tract Sampling (BTS) MCMC pretty-printers with objfile Obj."

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

def build_bts_mcmc_dictionary ():
    # bts objects requiring pretty-printing.
    # In order from:
    pretty_printers_dict[re.compile('^BTS::MCMC::Proposal::Momentum$')] = lambda val: ProposalMomentumPrinter("MCMC::Proposal::Momentum", val)
    pretty_printers_dict[re.compile('^BTS::MCMC::Proposal::Momentum::Weighted$')] = lambda val: ProposalMomentumPrinter("MCMC::Proposal::Momentum::Weighted", val)
    pretty_printers_dict[re.compile('^BTS::MCMC::Proposal::Momentum::Weighted::NonSeparable$')] = lambda val: ProposalMomentumPrinter("MCMC::Proposal::Momentum::Weighted::NonSeparable", val)
            
pretty_printers_dict = {}

build_bts_mcmc_dictionary ()
