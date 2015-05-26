
 
#=====================================================================================================================
# Test functions
#=====================================================================================================================
 
from fibre.base.object import Object
from fibre.base.set import Set, CantBeFrozenException, ElemSizeException, UnmatchedPropertiesException

import os.path
import random
import time

PROP_1 = "PROP_1"
PROP_2 = "PROP_2"
PROP_3 = "PROP_3"
PROP_4 = "PROP_4"
PROP_5 = "PROP_5"
PROP_6 = "PROP_6"
PROP_7 = "PROP_7"
PROP_8 = "PROP_8"
PROP_9 = "PROP_9"
PROP_10 = "PROP_10"

ELEM_PROP_1 = "ELEM_PROP_1"
ELEM_PROP_2 = "ELEM_PROP_2"
ELEM_PROP_3 = "ELEM_PROP_3"
ELEM_PROP_4 = "ELEM_PROP_4"
ELEM_PROP_5 = "ELEM_PROP_5"
ELEM_PROP_6 = "ELEM_PROP_6"
ELEM_PROP_7 = "ELEM_PROP_7"
ELEM_PROP_8 = "ELEM_PROP_8"
ELEM_PROP_9 = "ELEM_PROP_9"
ELEM_PROP_10 = "ELEM_PROP_10"

EXT_ELEM_PROP_1 = "EXT_ELEM_PROP_1"
EXT_ELEM_PROP_2 = "EXT_ELEM_PROP_2"
EXT_ELEM_PROP_3 = "EXT_ELEM_PROP_3"
EXT_ELEM_PROP_4 = "EXT_ELEM_PROP_4"
EXT_ELEM_PROP_5 = "EXT_ELEM_PROP_5"
EXT_ELEM_PROP_6 = "EXT_ELEM_PROP_6"
EXT_ELEM_PROP_7 = "EXT_ELEM_PROP_7"
EXT_ELEM_PROP_8 = "EXT_ELEM_PROP_8"
EXT_ELEM_PROP_9 = "EXT_ELEM_PROP_9"
EXT_ELEM_PROP_10 = "EXT_ELEM_PROP_10"


class EmptyListException(Exception):
  
  def __init__(self, value):
    self.value = value

  def __str__(self):
    return repr(self.value)
  
 


def random_size():
  
  return random.choice(range(10))

def random_value():
  
  return random.choice(range(1000))

 
 
class StressTestGenerator:
  
  SET_SIZE = 5
  OBJECT_SIZE = 5
  NUM_TESTS = 500
  
  def __init__(self, num_elems, elem_degree): 
    
    self.set = Set()
    self.set.first_set(num_elems, elem_degree)
    self.orig_size = num_elems
    
    
  def get_valid_elem_degree(self):
    
    if self.only_valid and self.set.elem_degree:
      elem_degree = self.set.elem_degree
    else:
      elem_degree = random_size()
    
    return elem_degree
    
    
  def get_valid_elem_index(self):

    if self.only_valid:
      if not self.set.size():
        raise EmptyListException('')
      elem_index = random.choice(xrange(self.set.size()))
    else:
      elem_index = random_size()
      
    return elem_index
      
      
  def get_valid_extend_elem_key(self): 
       
    if self.only_valid:
      if not self.set.elem_extend_prps.keys():
        raise EmptyListException('')  
      ext_elem_key = random.choice(self.set.elem_extend_prps.keys())   
    else:
      ext_elem_key = 'EXT_ELEM_PROP_%d' % random.choice(xrange(10)+1)
    
    return ext_elem_key
  
  
  def get_valid_elem_key(self): 
    
    if self.only_valid:
      if not self.set.elem_prp_keys:
        raise EmptyListException('')      
      elem_key = random.choice(self.set.elem_prp_keys)   
    else:
      elem_key = 'ELEM_PROP_%d' % random.choice(xrange(10)+1)
    
    return elem_key
  
  
  def get_valid_elem_key_index(self):
    
    if self.only_valid:
      if not self.set.elem_prp_keys:
        raise EmptyListException('')   
      elem_key_index = random.choice(xrange(len(self.set.elem_prp_keys)))   
    else:
      elem_key_index = random.choice(xrange(10))
    
    return elem_key_index
  
    
  def get_valid_key(self):
    
    if self.only_valid:
      if not self.set.prps.keys():
        raise EmptyListException('')         
      key = random.choice(self.set.prps.keys())   
    else:
      key = 'PROP_%d' % random.choice(xrange(10)+1)
    
    return key
    
    
  def get_invalid_key(self):
    
    if self.only_valid:
      invalid_keys = list()
      for i in range(1,11):
        key = 'PROP_%d' % i
        if key not in self.set.prps.keys():
          invalid_keys.append(key)
            
      if not invalid_keys:
        raise EmptyListException('')         
      key = random.choice(invalid_keys)   
    else:
      key = 'PROP_%d' % random.choice(xrange(10)+1)
    
    return key    
    
  def get_valid_key_index(self):
    
    if self.only_valid:
      if not self.set.prps:
        raise EmptyListException('')            
      key_index = random.choice(xrange(len(self.set.prps)))   
    else:
      key_index = random.choice(xrange(10))
    
    return key_index
      
    
  def clear(self):  
    self.size = 0
    
    
  def elem_resize_test(self):
    
    format_str = 'set.elem_resize(%d,%d)'
    
    self.elem_degree = random.choice(range(1,10))
    
    function_str = format_str % (self.elem_degree, random_value())
    
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,False)
    
  
  def push_back_test(self):
    
    format_str = 'set.push_back(objects[%d][%d])'
          
    function_str = format_str % (self.get_valid_elem_degree(), random_size())
    
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,False)


  def push_back_w_row_test(self):

    format_str = 'set.push_back(objects[%d][%d], ext_prop_row_maps[%d])'
            
    function_str = format_str % (self.get_valid_elem_degree(), random_size(), random_size())
        
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,False)
    
  
  def append_test(self):
    
    function_str = 'set.append(set2)'

    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,False)
  
  
  def erase_test(self):
    
    format_str = 'set.erase(%d)'
    
    function_str = format_str % self.get_valid_elem_index()
        
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,False)
  
    
  def insert_test(self):
    
    format_str = 'set.insert(objects[%d][%d],%d)'
    
    function_str = format_str % (self.get_valid_elem_degree(),random_size(),self.get_valid_elem_index())
        
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,False)
  
  
  def clear_test(self):
    
    format_str = 'set.clear()'
    
    function_str = format_str
    
    self.size = 0
    
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,False)
  
  
  def num_extend_elem_props_test(self):
    
    format_str = 'set.num_extend_elem_props()'
    
    function_str = format_str
    
    return (function_str,True,False)
    
  
  def get_extend_elem_prop_header_test(self):
    
    format_str = 'set.get_extend_elem_prop_header()'
        
    function_str = format_str
    
    return (function_str,True,False)
  
    
  def add_extend_elem_prop_test(self):
    
    format_str = 'set.add_extend_elem_prop(EXT_ELEM_PROP_%d,"-99")'
    
    function_str = format_str % (random_size()+1)
    
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,False)
  
  
  def add_extend_elem_props_test(self):
    
    format_str = 'set.add_extend_elem_props(set2)'
    
    function_str = format_str
    
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,False)
  
  
  def clear_extend_elem_props_test(self):
    
    format_str = 'set.clear_extend_elem_props()'
    
    function_str = format_str
    
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,False)
  
  
  def remove_extend_elem_prop_test(self):
    
    format_str = 'set.remove_extend_elem_prop(%s)'
    
    function_str = format_str % self.get_valid_extend_elem_key()
    
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,False)
  
  
  def set_extend_elem_prop_test(self):
    
    format_str = 'set.set_extend_elem_prop(%s,"%d",%d)'
    
    function_str = format_str % (self.get_valid_extend_elem_key(),random_value(),self.get_valid_elem_index())
    
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,False)
  
  
  def get_extend_elem_prop_test(self):
    
    format_str = 'set.get_extend_elem_prop(%s,%d)'
    
    function_str = format_str % (self.get_valid_extend_elem_key(),self.get_valid_elem_index())
    
    return (function_str,True,False)
    
  
  def has_extend_elem_prop_test(self):
    
    format_str = 'set.has_extend_elem_prop(EXT_ELEM_PROP_%d)'
    
    function_str = format_str % (random_size()+1)
    
    return (function_str,True,False)
  
  
  def get_extend_elem_prop_row_test(self):
    
    format_str = 'set.get_extend_elem_prop_row(%d)'
    
    function_str = format_str % self.get_valid_elem_index()
    
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,False)
  
  
  def set_extend_elem_prop_row_test(self):
    
    format_str = 'set.set_extend_elem_prop_row(ext_prop_row_maps[%d],%d)'
    
    function_str = format_str % (random_size(),self.get_valid_elem_index())
    
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,False)
  
  
  def copy_extend_elem_prop_row_test(self):
    
    format_str = 'set.copy_extend_elem_prop_row(set2,%d,%d)'
        
    function_str = format_str % (random.choice(range(self.orig_size)),self.get_valid_elem_index())
    
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,False)

  
  def copy_extend_elem_props_test(self):
    
    if self.set.size() != set2.size():
      raise ElemSizeException("'copy_extend_elem_props_test' could not be tested because set size has changed")
    
    format_str = 'set.copy_extend_elem_props(set2)'
    
    function_str = format_str
    
    if self.only_valid:
      eval('self.' + function_str)
  
      
    return (function_str,False,False)
  
  
  def append_extend_elem_props_test(self):
    
    if self.set.size() != set2.size():
      raise ElemSizeException("'append_extend_elem_props_test' could not be tested because set size has changed")
    
    format_str = 'set.append_extend_elem_props(set2)'
    
    function_str = format_str
    
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,False)

  
  def num_elem_props_test(self):
    
    format_str = 'set.num_elem_props()'
        
    function_str = format_str
    
    return (function_str,True,False)
  
  
  def add_elem_prop_test(self):
    
    format_str = 'set.add_elem_prop(%s,%d)'
    
    function_str = format_str % (self.get_valid_elem_key(),random_value())
    
    return (function_str,False,True)
  
  
  def has_elem_prop_test(self):
    
    format_str = 'set.has_elem_prop(ELEM_PROP_%d)'
        
    function_str = format_str % (random_size()+1)
    
    return (function_str,True,False)
  
  
  def remove_elem_prop_test(self):
    
    format_str = 'set.remove_elem_prop(%s)'
    
    function_str = format_str % self.get_valid_elem_key()
    
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,True)
  
  
  def clear_elem_props_test(self):
    
    format_str = 'set.clear_elem_props()'
    
    function_str = format_str
    
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,True)
  
  
  def elem_prop_index_test(self):
    
    format_str = 'set.elem_prop_index(%s)'
    
    function_str = format_str % self.get_valid_elem_key()
    
    return (function_str,True,False)
  
  
  def elem_prop_name_test(self):
    
    format_str = 'set.elem_prop_name(%d)'
        
    function_str = format_str % self.get_valid_elem_key_index()
    
    return (function_str,True,False)
  
  
  def elem_prop_keys_test(self):
    
    format_str = 'set.elem_prop_keys()'
        
    function_str = format_str
    
    return (function_str,True,False)
  
  
  def freeze_elem_degree_test(self):
    
    format_str = 'set.freeze_elem_degree()'
    
    function_str = format_str
    
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,True)
  
  
  def free_elem_degree_test(self):
    
    format_str = 'set.free_elem_degree()'
    
    function_str = format_str
    
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,False)
  
  
  def resize_test(self):
    
    format_str = 'set.resize(%d,%d,%d)'
    
    self.size = random_size()
    
    if self.set.elem_degree:
      elem_degree = 0
    else:
      elem_degree = self.get_valid_elem_degree()
      
    function_str = format_str % (random_size(),random_value(),elem_degree)
    
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,False)
    
    
  def size_test(self):

    format_str = 'set.size()'
    
    function_str = format_str
    
    return (function_str,True,False)
  
  
  def vsize_test(self):

    format_str = 'set.vsize()'
    
    function_str = format_str
    
    return (function_str,True,False)
  
  
  def bsize_test(self):

    format_str = 'set.bsize()'
    
    function_str = format_str
    
    return (function_str,True,False)
  
  
  def num_props_test(self):

    format_str = 'set.num_props()'
    
    function_str = format_str
    
    return (function_str,True,False)
  
  
  def add_prop_test(self):

    format_str = 'set.add_prop(%s,%d)'
    
    function_str = format_str % (self.get_invalid_key(),random_value())
    
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,False)
  
  
  def has_prop_test(self):

    format_str = 'set.has_prop(PROP_%d)'
    
    function_str = format_str % (random_size()+1)
    
    return (function_str,True,False)
  
  
  def remove_prop_test(self):

    format_str = 'set.remove_prop(%s)'
    
    function_str = format_str % self.get_valid_key()
    
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,False)
  
  
  def clear_props_test(self):

    format_str = 'set.clear_props()'
    
    function_str = format_str
    
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,False)
  
  
  def prop_index_test(self):

    format_str = 'set.prop_index(%s)'
    
    function_str = format_str % self.get_valid_key()
    
    return (function_str,True,False)
  
  
  def prop_test(self):

    format_str = 'set.prop(%s)'
    
    function_str = format_str % self.get_valid_key()
    
    return (function_str,True,False)
  
  
  def prop_name_test(self):

    format_str = 'set.prop_name(%d)'
    
    function_str = format_str % self.get_valid_key_index()
    
    if self.only_valid:
      eval('self.' + function_str)
    
    return (function_str,False,False)
  
  
  def prop_keys_test(self):

    format_str = 'set.prop_keys()'
    
    function_str = format_str
    
    return (function_str,True,False)


  def generate_stress_test(self,num_iterations, only_valid=True, include_clear=True, resize_on_freeze=True, print_every=10):
    
    #Store only valid in object so it can be accessed by '*_test()' functions
    self.only_valid = only_valid
    
    # Set up function list
    function_list = list()
    function_list.append('self.elem_resize_test()')
    function_list.append('self.push_back_test()')
    function_list.append('self.push_back_w_row_test()')
    function_list.append('self.append_test()')
    function_list.append('self.erase_test()')
    function_list.append('self.insert_test()')
    function_list.append('self.num_extend_elem_props_test()')
    function_list.append('self.get_extend_elem_prop_header_test()')
    function_list.append('self.add_extend_elem_prop_test()')
    function_list.append('self.add_extend_elem_props_test()')
    function_list.append('self.remove_extend_elem_prop_test()')
    function_list.append('self.set_extend_elem_prop_test()')
    function_list.append('self.get_extend_elem_prop_test()')
    function_list.append('self.has_extend_elem_prop_test()')
    function_list.append('self.get_extend_elem_prop_row_test()')
    function_list.append('self.set_extend_elem_prop_row_test()')
    function_list.append('self.copy_extend_elem_prop_row_test()')
    function_list.append('self.copy_extend_elem_props_test()')
    function_list.append('self.append_extend_elem_props_test()')
    function_list.append('self.num_elem_props_test()')
    function_list.append('self.add_elem_prop_test()')
    function_list.append('self.has_elem_prop_test()')
    function_list.append('self.remove_elem_prop_test()')
    function_list.append('self.elem_prop_index_test()')
    function_list.append('self.elem_prop_name_test()')
    function_list.append('self.elem_prop_keys_test()')
    function_list.append('self.freeze_elem_degree_test()')
    function_list.append('self.free_elem_degree_test()')
    function_list.append('self.resize_test()')    
    function_list.append('self.size_test()')
    function_list.append('self.vsize_test()')
    function_list.append('self.bsize_test()')
    function_list.append('self.num_props_test()')
    function_list.append('self.add_prop_test()')
    function_list.append('self.has_prop_test()')
    function_list.append('self.remove_prop_test()')
    function_list.append('self.prop_index_test()')
    function_list.append('self.prop_test()')
    function_list.append('self.prop_name_test()')
    function_list.append('self.prop_keys_test()')
    
    if include_clear:
      function_list.append('self.clear_test()')
      function_list.append('self.clear_extend_elem_props_test()')
      function_list.append('self.clear_elem_props_test()')
      function_list.append('self.clear_props_test()')

  
    cpp_code = open('/home/tclose/Code/Tractography/cmd/stress_tester.cpp', 'w')
    py_code = open('/home/tclose/Code/Tractography/python/fibre/base/test/stress_tester.py','w')
    
    cpp_code.write(cpp_header)
    py_code.write(py_header)

    cpp_code.write('  out << "\\n\\n-------------------------------------\\n            INITIAL SET             \\n-------------------------------------\\n\\n" << set << "\\n\\n";\n') # Write print statement in C++
    py_code.write  ('f.write("\\n\\n-------------------------------------\\n            INITIAL SET             \\n-------------------------------------\\n" + str(set) + "\\n\\n")\n')            
    cpp_code.write('  out << "\\n\\n-------------------------------------\\n            SECOND SET             \\n-------------------------------------\\n\\n" << set2 << "\\n\\n";\n') # Write print statement in C++
    py_code.write  ('f.write("\\n\\n-------------------------------------\\n            SECOND SET             \\n-------------------------------------\\n" + str(set2) + "\\n\\n")\n')            


    seed = time.time()

    seed = 1323911724.74
    random.seed(seed)

    print "random seed: " + str(seed)
  
    modify_count = 0
    for iter in range(num_iterations):
      test_function = random.choice(function_list)
      print test_function
      try:
        (function_str, printable, elem_props_changed)= eval(test_function)
        
        if printable:
          cpp_code.write('  out << "' + function_str + ': " << ' + function_str + ' << std::endl;\n') # Write the function to be evaluated in C++
          py_code.write("f.write('" + function_str + ": ' + str(" + function_str + ").replace(\"', '\",\" \").replace(\"['\",\"[ \").replace(\"']\",\" ]\") + '\\n')\n") # Write the function to be evaluated in Python
        else:
          modify_count = modify_count + 1
          cpp_code.write('  ' + function_str + ';\n') # Write the function to be evaluated in C++
          py_code.write(function_str + '\n') # Write the function to be evaluated in Python
          
          if not (modify_count % print_every):
            cpp_code.write('  out << "\\n\\n-------------------------------------\\n            %d MODFICATIONS             \\n-------------------------------------\\n\\n" << set << "\\n\\n";\n' % modify_count) # Write print statement in C++
            py_code.write ('f.write("\\n\\n-------------------------------------\\n            %d MODFICATIONS             \\n-------------------------------------\\n" + str(set) + "\\n\\n")\n' % modify_count)            

            
          
        if elem_props_changed:
          cpp_code.write('  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);\n')
          py_code.write('elem_props = set.get_elem_prop_header(); (objects,prop_row_maps,prop_rows) = recreate_elem_prop_objects(elem_props)\n')
          
      except EmptyListException:
        print 'Could not produce valid index as list was empty'
        
      except CantBeFrozenException:
        if resize_on_freeze:
          self.elem_resize_test()
          self.freeze_elem_degree_test()
        else:
          print 'Couldn''t freeze because rows were of different sizes'
          
      except ElemSizeException:
        print 'Couldn''t append set because element sizes had already changed'
        
      except UnmatchedPropertiesException:
        print 'Couldn''t append set because properties had already changed'
      
      
    cpp_code.write(cpp_footer)
    py_code.write(py_footer)  
            
    cpp_code.close()
    py_code.close()



#=====================================================================================================================
# Stress test functions
#=====================================================================================================================
 

def recreate_objects(elem_keys):

  object_list = list()

  for size_i in range(10):
    objects = list()
    
    for val_i in range(10):
      
      object = Object(size_i, val_i + 1001)
      
      for key in elem_keys:  
        object.prps[key] = val_i + 1001

      objects.append(object)
      
    object_list.append(objects)

  return object_list


def recreate_prop_row_maps(elem_keys):

  row_maps = list()

  key_i = 0
  for i in range(10):
    value_i = 0
    row_map = dict()
    for key in elem_keys:
      row_map[key] = value_i * 0.01 + key_i * 0.1 + 0.001
      value_i = value_i + 1
    row_maps.append(row_map)
    key_i = key_i + 1

  return row_maps


def recreate_prop_rows(prop_keys):

  prop_rows = list()

  key_i = 0
  for i in range(10):
    value_i = 0
    prop_row = list()
    for key in prop_keys:
      prop_row.append(value_i * 0.01 + key_i * 0.1 + 0.001)
      value_i = value_i + 1
    prop_rows.append(prop_row)

    key_i = key_i + 1

  return prop_rows


def recreate_elem_prop_objects(prop_keys):
  
  return (recreate_objects(prop_keys),recreate_prop_row_maps(prop_keys),recreate_prop_rows(prop_keys))
  
  
def create_ext_prop_row_maps():

  row_maps = list()

  for value_i in range(1,11):
    row_map = dict()
    for key_i in range(1,11):
      row_map["EXT_ELEM_PROP_" + str(key_i)] = value_i * 0.01 + key_i * 1.1 + 0.001
    row_maps.append(row_map)

  return row_maps;
  
  
props = list()

props.append(PROP_1)
props.append(PROP_2)
props.append(PROP_3)
props.append(PROP_4)
props.append(PROP_5)

elem_props = list()

elem_props.append(ELEM_PROP_1)
elem_props.append(ELEM_PROP_2)
elem_props.append(ELEM_PROP_3)
elem_props.append(ELEM_PROP_4)
elem_props.append(ELEM_PROP_5)

objects = recreate_objects(elem_props)
prop_rows = recreate_prop_rows(elem_props)
prop_row_maps = recreate_prop_row_maps(elem_props)
ext_prop_row_maps = create_ext_prop_row_maps()  
  
set2 = Set()
set2.second_set(StressTestGenerator.SET_SIZE,StressTestGenerator.OBJECT_SIZE)



cpp_header = '\n\
\n\
extern "C" {\n\
#include <gsl/gsl_rng.h>\n\
#include <gsl/gsl_randist.h>\n\
}\n\
\n\
#include "math/matrix.h"\n\
#include "exception.h"\n\
\n\
#include "bts/cmd.h"\n\
#include "bts/common.h"\n\
#include "bts/file.h"\n\
\n\
#include "bts/fibre/base/set.h"\n\
\n\
\n\
using namespace FTS;\n\
\n\
\n\
std::vector< std::vector<Fibre::Base::Object> >  recreate_objects(const std::vector<const char*>& elem_props) {\n\
\n\
  std::vector< std::vector<Fibre::Base::Object> > objects(10);\n\
\n\
  for (size_t size_i = 0; size_i < 10; ++size_i) {\n\
    objects[size_i];\n\
    for (size_t value_i = 0; value_i < 10; ++value_i) {\n\
      objects[size_i].push_back(Fibre::Base::Object(size_i, 3, elem_props));\n\
      objects[size_i][value_i].set(value_i+1001);\n\
    }\n\
  }\n\
\n\
  return objects;\n\
\n\
}\n\
\n\
std::vector< std::map<const char*, double> >  recreate_prop_row_maps(const std::vector<const char*>& elem_props) {\n\
\n\
  std::vector< std::map<const char*, double> > prop_row_maps(10);\n\
\n\
  for (size_t value_i = 0; value_i < 10; ++value_i) {\n\
    for (size_t key_i = 0; key_i < elem_props.size(); ++key_i)\n\
      prop_row_maps[value_i][elem_props[key_i]] = (double)(value_i+1) * 0.01 + (double)(key_i+1) * 0.1 + 0.001;\n\
  }\n\
\n\
  return prop_row_maps;\n\
\n\
}\n\
\n\
std::vector< std::vector<double> >  recreate_prop_rows(const std::vector<const char*>& elem_props) {\n\
\n\
  std::vector< std::vector<double> > prop_rows(10);\n\
\n\
  for (size_t value_i = 0; value_i < 10; ++value_i) {\n\
    for (size_t key_i = 0; key_i < elem_props.size(); ++key_i)\n\
      prop_rows[value_i].push_back((double)(value_i+1) * 0.01 + (double)(key_i+1) * 0.1 + 0.001);\n\
  }\n\
\n\
  return prop_rows;\n\
\n\
}\n\
\n\
\n\
std::vector< std::map<std::string,std::string> > create_ext_prop_row_maps() {\n\
\n\
  std::vector< std::map<std::string,std::string> > prop_rows(10);\n\
\n\
  for (size_t value_i = 0; value_i < 10; ++value_i) {\n\
    for (size_t key_i = 1; key_i < 11; ++key_i)\n\
      prop_rows[value_i]["EXT_ELEM_PROP_" +str(key_i)] = str((double)(value_i+1) * 0.01 + (double)key_i * 1.1 + 0.001);\n\
  }\n\
\n\
  return prop_rows;\n\
\n\
}\n\
\n\
void recreate_elem_prop_objects(const std::vector<const char*>& elem_props, std::vector<std::vector<Fibre::Base::Object> >& objects, std::vector< std::map<const char*, double> >& prop_row_maps, std::vector< std::vector<double> >& prop_rows) {\n\
\n\
  objects = recreate_objects(elem_props);\n\
  prop_row_maps = recreate_prop_row_maps(elem_props);\n\
  prop_rows = recreate_prop_rows(elem_props);\n\
\n\
}\n\
\n\
SET_VERSION_DEFAULT;\n\
SET_AUTHOR ("Thomas G. Close");\n\
SET_COPYRIGHT (NULL);\n\
\n\
DESCRIPTION = {\n\
  "Tests a random sequence of Fibre::Base::Set operations against those performed on a Python version.",\n\
  "",\n\
  NULL\n\
};\n\
\n\
ARGUMENTS = {\n\
\n\
\n\
  Argument()\n\
};\n\
\n\
\n\
OPTIONS = {\n\
\n\
Option() };\n\
\n\
\n\
EXECUTE {\n\
\n\
  size_t SET_SIZE = 5;\n\
  size_t OBJECT_SIZE = 5;\n\
\n\
  const char* PROP_1 = "PROP_1";\n\
  const char* PROP_2 = "PROP_2";\n\
  const char* PROP_3 = "PROP_3";\n\
  const char* PROP_4 = "PROP_4";\n\
  const char* PROP_5 = "PROP_5";\n\
  const char* PROP_6 = "PROP_6";\n\
  const char* PROP_7 = "PROP_7";\n\
  const char* PROP_8 = "PROP_8";\n\
  const char* PROP_9 = "PROP_9";\n\
  const char* PROP_10 = "PROP_10";\n\
\n\
  const char* ELEM_PROP_1 = "ELEM_PROP_1";\n\
  const char* ELEM_PROP_2 = "ELEM_PROP_2";\n\
  const char* ELEM_PROP_3 = "ELEM_PROP_3";\n\
  const char* ELEM_PROP_4 = "ELEM_PROP_4";\n\
  const char* ELEM_PROP_5 = "ELEM_PROP_5";\n\
  const char* ELEM_PROP_6 = "ELEM_PROP_6";\n\
  const char* ELEM_PROP_7 = "ELEM_PROP_7";\n\
  const char* ELEM_PROP_8 = "ELEM_PROP_8";\n\
  const char* ELEM_PROP_9 = "ELEM_PROP_9";\n\
  const char* ELEM_PROP_10 = "ELEM_PROP_10";\n\
\n\
  std::string EXT_ELEM_PROP_1 = "EXT_ELEM_PROP_1";\n\
  std::string EXT_ELEM_PROP_2 = "EXT_ELEM_PROP_2";\n\
  std::string EXT_ELEM_PROP_3 = "EXT_ELEM_PROP_3";\n\
  std::string EXT_ELEM_PROP_4 = "EXT_ELEM_PROP_4";\n\
  std::string EXT_ELEM_PROP_5 = "EXT_ELEM_PROP_5";\n\
  std::string EXT_ELEM_PROP_6 = "EXT_ELEM_PROP_6";\n\
  std::string EXT_ELEM_PROP_7 = "EXT_ELEM_PROP_7";\n\
  std::string EXT_ELEM_PROP_8 = "EXT_ELEM_PROP_8";\n\
  std::string EXT_ELEM_PROP_9 = "EXT_ELEM_PROP_9";\n\
  std::string EXT_ELEM_PROP_10 = "EXT_ELEM_PROP_10";\n\
\n\
  std::vector<const char*> props;\n\
\n\
  props.push_back(PROP_1);\n\
  props.push_back(PROP_2);\n\
  props.push_back(PROP_3);\n\
  props.push_back(PROP_4);\n\
  props.push_back(PROP_5);\n\
\n\
  std::vector<const char*> elem_props;\n\
\n\
  elem_props.push_back(ELEM_PROP_1);\n\
  elem_props.push_back(ELEM_PROP_2);\n\
  elem_props.push_back(ELEM_PROP_3);\n\
  elem_props.push_back(ELEM_PROP_4);\n\
  elem_props.push_back(ELEM_PROP_5);\n\
\n\
\n\
  std::vector< std::vector<Fibre::Base::Object> > objects = recreate_objects(elem_props);\n\
  std::vector<std::vector<double> > prop_rows = recreate_prop_rows(elem_props);\n\
  std::vector< std::map<const char*, double> > prop_row_maps = recreate_prop_row_maps(elem_props);\n\
\n\
  std::vector< std::map<std::string,std::string> > ext_prop_row_maps = create_ext_prop_row_maps();\n\
\n\
\n\
  //Iniate set\n\
  Fibre::Base::Set<Fibre::Base::Object> set(SET_SIZE, OBJECT_SIZE, 3, props, elem_props);\n\
\n\
  for (size_t elem_i = 0; elem_i < SET_SIZE; ++elem_i) {\n\
    Fibre::Base::Object ob = set[elem_i];\n\
    for (size_t elem_elem_i = 0; elem_elem_i < OBJECT_SIZE; ++elem_elem_i) {\n\
      Triple<double>& t = ob[elem_elem_i];\n\
      for (size_t elem_elem_elem_i = 0; elem_elem_elem_i < 3; ++elem_elem_elem_i)\n\
        t[elem_elem_elem_i] = elem_i * 100 + elem_elem_i * 10 + elem_elem_elem_i +1;\n\
    }\n\
\n\
    for (size_t prop_i = 0; prop_i < 5; ++prop_i)\n\
      ob.prop(prop_i) = (double)elem_i * 100.0 + (double)(prop_i+1) * .1 + 0.001;\n\
  }\n\
\n\
  for (size_t prop_i = 0; prop_i < 5; ++prop_i)\n\
    set.prop(prop_i) = (double)(prop_i+1) * 1.1 + 0.001;\n\
\n\
  for (size_t eprop_i = 1; eprop_i < 6; ++eprop_i) {\n\
    set.add_extend_elem_prop("EXT_ELEM_PROP_" + str(eprop_i),"NaN");\n\
    for (size_t row_i = 0; row_i < set.size(); ++row_i)\n\
      set.set_extend_elem_prop("EXT_ELEM_PROP_" + str(eprop_i), "Ooogle " + str(eprop_i) + "-" + str(row_i), row_i);\n\
  }\n\
\n\
\n\
  //Iniate set\n\
  Fibre::Base::Set<Fibre::Base::Object> set2(SET_SIZE, OBJECT_SIZE, 3, props, elem_props);\n\
\n\
  for (size_t elem_i = 0; elem_i < SET_SIZE; ++elem_i) {\n\
    Fibre::Base::Object ob = set2[elem_i];\n\
    for (size_t elem_elem_i = 0; elem_elem_i < OBJECT_SIZE; ++elem_elem_i) {\n\
      Triple<double>& t = ob[elem_elem_i];\n\
      for (size_t elem_elem_elem_i = 0; elem_elem_elem_i < 3; ++elem_elem_elem_i)\n\
        t[elem_elem_elem_i] = -(double)(elem_i * 100 + elem_elem_i * 10 + elem_elem_elem_i + 1);\n\
    }\n\
\n\
    for (size_t prop_i = 0; prop_i < 5; ++prop_i)\n\
      ob.prop(prop_i) = -(double)elem_i * 100.0 - (double)(prop_i+1) * .1 - 0.001;\n\
  }\n\
\n\
  for (size_t prop_i = 0; prop_i < 5; ++prop_i)\n\
    set2.prop(prop_i) = -(double)(prop_i+1) * 1.1 - 0.001;\n\
\n\
  for (size_t eprop_i = 6; eprop_i < 11; ++eprop_i) {\n\
    set2.add_extend_elem_prop("EXT_ELEM_PROP_" + str(eprop_i),"NaN");\n\
    for (size_t row_i = 0; row_i < set2.size(); ++row_i)\n\
      set2.set_extend_elem_prop("EXT_ELEM_PROP_" + str(eprop_i), "Boogle " + str(eprop_i) + "-" + str(row_i), row_i);\n\
  }\n\
\n\
  std::ofstream out;\n\
\n\
  out.open("/home/tclose/Code/Tractography/python/fibre/base/test/output/cpp_output.txt");\n\
\n\
\n\
\n\
  //====================================================================================================================\n\
  // Test functions start here\n\
  //====================================================================================================================\n\
\n'


cpp_footer='\
\n\
\n\
  out << std::endl << std::endl << "-------------------------------------" << std::endl;\n\
  out << "             FINAL PRINT             " << std::endl;\n\
  out << "-------------------------------------\\n" << std::endl;\n\
\n\
  out << set << std::endl;\n\
\n\
  std::cout << "Ran cpp test successfully" << std::endl;\n\
\n\
\n\
}\n'



py_header="\
from stress_test_generator import *\n\
\n\
set = Set()\n\
set.first_set(StressTestGenerator.SET_SIZE,StressTestGenerator.OBJECT_SIZE)\n\
global objects\n\
global prop_rows\n\
global prop_row_maps\n\
\n\
f = open('/home/tclose/Code/Tractography/python/fibre/base/test/output/python_output.txt','w')\n"



#=====================================================================================================================
# Generated functions start here
#=====================================================================================================================


py_footer="\
\n\
\n\
#=====================================================================================================================\n\
# Generated functions end here\n\
#=====================================================================================================================\n\
\n\
f.write('\\n\\n-------------------------------------\\n             FINAL PRINT             \\n-------------------------------------\\n')\n\
\n\
f.write(str(set))\n\
\n\
print 'Ran python test successfully!'\n\
\n\
f.write('\\n')\n\
\n\
f.close()\n"


if __name__ == '__main__':
  
  test_generator = StressTestGenerator(StressTestGenerator.SET_SIZE, StressTestGenerator.OBJECT_SIZE)
  test_generator.generate_stress_test(StressTestGenerator.NUM_TESTS, include_clear=False, print_every=5)

    
    