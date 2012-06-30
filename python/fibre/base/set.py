from copy import copy, deepcopy
from object import Object

class CantBeFrozenException(Exception):
  
  def __init__(self, value):
    self.value = value

  def __str__(self):
    return repr(self.value)


class ElemSizeException(Exception):
  
  def __init__(self, value):
    self.value = value

  def __str__(self):
    return repr(self.value)
  
  
class UnmatchedPropertiesException(Exception):
  
  def __init__(self, value):
    self.value = value

  def __str__(self):
    return repr(self.value)  


class Set(object):
  
  def __init__(self):
    
    self.elems = list()
    
    self.prps = dict()
    self.extend_prps = dict()

    self.elem_prp_keys = list()
    self.elem_extend_prps = dict()
      
    self.elem_degree = 0
    
    
  def first_set(self, num_elems, elem_degree):
    
    self.elems = list()
    
    self.prps = dict()
    self.extend_prps = dict()

    self.elem_prp_keys = list()
    self.elem_extend_prps = dict()
      
    for elem_i in range(num_elems):
      
      self.elems.append(Object(elem_degree))
      
      for elem_elem_i in range(elem_degree):
        base = elem_i * 100 + elem_elem_i * 10
        self.elems[elem_i].vector[elem_elem_i] = (base + 1, base + 2, base + 3)
        
        for prop_i in range(5):
          self.elems[elem_i].prps['ELEM_PROP_' + str(prop_i+1)] = elem_i * 100.0 + (prop_i+1) * .1 + 0.001;
      
    for prop_i in range(5):  
      self.elem_prp_keys.append('ELEM_PROP_' + str(prop_i+1))  
      self.prps['PROP_' + str(prop_i+1)] = (prop_i+1) * 1.1 + 0.001
      
      
    for i in range(1,6):
      self.add_extend_elem_prop("EXT_ELEM_PROP_" + str(i),"NaN")
      for j in range(len(self.elems)):
        self.set_extend_elem_prop("EXT_ELEM_PROP_" + str(i), "Ooogle " + str(i) + "-" + str(j), j)  
   
      
    self.elem_degree = elem_degree  
        

  def second_set(self, num_elems, elem_degree):
    
    self.elems = list()
    
    self.prps = dict()
    self.extend_prps = dict()

    self.elem_prp_keys = list()
    self.elem_extend_prps = dict()
      
    for elem_i in range(num_elems):
      
      self.elems.append(Object(elem_degree))
      
      for elem_elem_i in range(elem_degree):
        base = elem_i * 100 + elem_elem_i * 10
        self.elems[elem_i].vector[elem_elem_i] = (-(base + 1), -(base + 2), -(base + 3))
        
        for prop_i in range(5):
          self.elems[elem_i].prps['ELEM_PROP_' + str(prop_i+1)] = -(elem_i * 100.0 + (prop_i+1) * .1 + 0.001);
      
    for prop_i in range(5):  
      self.elem_prp_keys.append('ELEM_PROP_' + str(prop_i+1))  
      self.prps['PROP_' + str(prop_i+1)] = -((prop_i+1) * 1.1 + 0.001)
      
      
    for i in range(6,11):
      self.add_extend_elem_prop("EXT_ELEM_PROP_" + str(i),"NaN")
      for j in range(len(self.elems)):
        self.set_extend_elem_prop("EXT_ELEM_PROP_" + str(i), "Boogle " + str(i) + "-" + str(j), j)  
   
      
    self.elem_degree = elem_degree  


  def __str__(self):
      
    string = ''  
    
    if self.extend_prps:
      string = string + '('  
    for key_val in sorted(self.extend_prps.items()):  
      string = string + key_val[0] + ': ' + str(key_val[1]) + ', ' 
    if self.extend_prps:
      string = string + ')\n'       
      
    string = string + '\n'    
    
    for (key,val) in sorted(self.prps.items()):
      string = string + key + ': ' + str(val) + '\n' 
        
    elem_i = 0  
    
    for elem in self.elems:
      string = string + '\n--------------------\n'
      string = string + '        %d\n' % elem_i
      string = string + '--------------------\n'    
      
      if elem.extend_prps:
        string = string + '{'
      count = 0
      for (key,val) in sorted(elem.extend_prps.items()):
        string = string + key + ': ' + str(val) 
        count = count + 1
        if count != len(elem.extend_prps):
          string = string + ', ' 
      if elem.extend_prps:
        string = string + '}\n'        
      
      string = string + '\n'
      
      string = string + str(elem)
  
      string = string + '\n'    
  
      elem_i = elem_i + 1
  
    return string
      
  #=====================================================================================================================
  # Set functions
  #=====================================================================================================================


  def elem_resize(self, new_size, value):
    
    for elem_i in range(self.size()):
      
      elem_degree = len(self.elems[elem_i].vector)
      
      for i in range(elem_degree,new_size):
        self.elems[elem_i].vector.append((value,value,value)) 
      
      for i in range(elem_degree-1,new_size-1,-1):
        self.elems[elem_i].vector.pop(i)
        
    self.elem_degree = new_size
      
  
  def push_back(self, elem, ext_prop_row=dict()):
    
    # Make a copy of the element and reset its extended properties
    append_elem = deepcopy(elem)      
    append_elem.extend_prps = dict()
    
    # Loop through the keys present in the set
    for key in self.elem_extend_prps.keys():
      if ext_prop_row:    
        # Use the properties explicitly provided  
        if ext_prop_row.has_key(key):
          append_elem.extend_prps[key] = ext_prop_row[key]
        else:
          append_elem.extend_prps[key] = self.elem_extend_prps[key]      
      else:
        # Copy across only the extended properties present in the set, and populate the missing extended properties with defaults.
        if elem.extend_prps.has_key(key):
          append_elem.extend_prps[key] = elem.extend_prps[key]
        else:
          append_elem.extend_prps[key] = self.elem_extend_prps[key]      
                    
    self.elems.append(append_elem)
    
  
  def erase(self, index):
    
    self.elems.pop(index)
    
    
  def insert(self, elem, index):  
    
    # Make a copy of the element and reset its extended properties
    insert_elem = deepcopy(elem)      
    insert_elem.extend_prps = dict()
    
    # Loop through the keys present in the set
    for key in self.elem_extend_prps.keys():
      # Copy across only the extended properties present in the set, and populate the missing extended properties with defaults.
      if elem.extend_prps.has_key(key):
        insert_elem.extend_prps[key] = elem.extend_prps[key]
      else:
        insert_elem.extend_prps[key] = self.elem_extend_prps[key]   
    
    self.elems.insert(index, insert_elem)
  
  def append(self, set):
    
    if self.elem_prp_keys != set.elem_prp_keys or self.prps != set.prps:
      raise UnmatchedPropertiesException ('properties do not match')
    
    if self.elem_degree and self.elem_degree != set.elem_degree:
      raise ElemSizeException ('Element sizes do not match')
    
    for elem in set.elems:
      self.push_back(elem)    
  
  
  def get_extend_elem_prop_header(self):
    
    return sorted(self.elem_extend_prps.keys())
 
 
  def get_elem_prop_header(self):
    
    return sorted(self.elem_prp_keys)
  
  
  def add_extend_elem_prop(self, key, default_value):
    
    if not self.elem_extend_prps.has_key(key):
      for elem in self.elems:        
        elem.extend_prps[key] = default_value

    self.elem_extend_prps[key] = default_value
  
  def add_extend_elem_props(self,set):
    
    for (key,val) in set.elem_extend_prps.items():
      self.add_extend_elem_prop(key,val)
    
  
  def clear_extend_elem_props(self):
    
    self.elem_extend_prps = dict()
    
    for elem in self.elems:
      
      elem.extend_prps = dict()

  
  def remove_extend_elem_prop(self, key):
    
    self.elem_extend_prps.pop(key,None)
    
    for elem in self.elems:
      
      elem.extend_prps.pop(key,None)

  #######################
  def set_extend_elem_prop(self,key,value,index):
    
    self.elems[index].extend_prps[key] = value 

  
  def get_extend_elem_prop(self,key,index):
    
    return self.elems[index].extend_prps[key]

  
  def has_extend_elem_prop(self,key):
    
    return self.elem_extend_prps.keys().count(key)

  
  def get_extend_elem_prop_row(self,index):
    
    return self.elems[index].extend_prps

  
  def set_extend_elem_prop_row(self,prop_row,index):
    
    this_props = self.elems[index].extend_prps
    
    # Loop through the keys present in the set
    for key in self.elem_extend_prps.keys():
      # Copy across only the extended properties present in the set, and populate the missing extended properties with defaults.
      if prop_row.has_key(key):
        this_props[key] = prop_row[key]
      else:
        this_props[key] = self.elem_extend_prps[key]
    

  def copy_extend_elem_prop_row(self,set,their_index,this_index):

    self.set_extend_elem_prop_row(set.get_extend_elem_prop_row(their_index),this_index)
    

  def copy_extend_elem_props(self,set):
    
    self.elem_extend_prps = copy(set.elem_extend_prps)
    
    for (this_elem,their_elem) in zip(self.elems,set.elems):
      this_elem.extend_prps = copy(their_elem.extend_prps)
      
  
  def append_extend_elem_props(self,set):
    
    for (key,val) in set.elem_extend_prps.items():
#      if key not in self.elem_extend_prps:
      self.elem_extend_prps[key] = val
      for (this_elem,their_elem) in zip(self.elems,set.elems):
        this_elem.extend_prps[key] = their_elem.extend_prps[key]

  
  def add_elem_prop(self,key,value):
    
    if key not in self.elem_prp_keys:
      self.elem_prp_keys.append(key)

    for elem in self.elems:
      
      elem.prps[key] = value
  
  
  def has_elem_prop(self,key):
    
    return self.elem_prp_keys.count(key)
    
  
  def remove_elem_prop(self,key,ignore_missing=True):
    
    try:
      self.elem_prp_keys.remove(key)
    except ValueError:
      if not ignore_missing:
        raise Exception ("Could not remove element property '%s' as it was not present" % key)
            
    for elem in self.elems:
      
      elem.prps.pop(key,None)

  
  def clear_elem_props(self):
    
    self.elem_prp_keys = list()

    for elem in self.elems:
      
      elem.prps = dict()
  
  
  def elem_prop_index(self,key):
    
    keys = sorted(self.elem_prp_keys)
    
    for i in range(len(keys)):
      if keys[i] == key:
        return i
      
    raise Exception ("Key (" + key + ") was not found.")

  
  def elem_prop_name(self, index):
    
    keys = sorted(self.elem_prp_keys)
    
    return keys[index]


  def elem_prop_keys(self):
    
    keys = sorted(self.elem_prp_keys)
    
    return keys    

  
  def freeze_elem_degree(self):
    
    if self.size():
      elem_degree = self.elems[0].size()
      
      for elem in self.elems:
        if elem.size() != elem_degree:
          raise CantBeFrozenException ("Elements not all the same size")
        
      self.elem_degree = elem_degree
        
        
  def free_elem_degree(self):
    
    self.elem_degree = 0


  
  def resize(self, new_size, fill_value, new_elem_degree=0):
    
    if self.elem_degree:
      if new_elem_degree:
        raise Exception ("'new_elem_degree' parameter can only be supplied when element size is variable.");
      elem_degree = self.elem_degree
    else:
      elem_degree = new_elem_degree
      
    old_size = self.size()
    
    for i in range(old_size,new_size):
      object = Object(elem_degree,fill_value)
      for key in self.elem_prp_keys:
        object.prps[key] = fill_value
        
      for (key,val) in self.elem_extend_prps.items():
        object.extend_prps[key] = val
      
      self.elems.append(object)

    for i in range(old_size-1,new_size-1,-1):
      self.elems.pop(i)

    
  def clear(self):
       
    self.elems = list()
    self.extend_prps = dict()
    self.elem_extend_prps = dict()
       
  
  def size(self):
    
    return len(self.elems)
  
  
  def num_props(self):
    
    return len(self.prps)
  
  
  def num_elem_props(self):
    
    return len(self.elem_prp_keys)
  
  
  def num_extend_elem_props(self):
    
    return len(self.elem_extend_prps)


  def vsize(self):
 
    return self.bsize() + self.num_props()


  def bsize(self):
    
    base_size = 0
    
    for elem in self.elems:
      base_size = base_size + elem.size() * 3 + elem.num_props()
    
    return base_size


  def add_prop(self, prop_key,value):
    
    if self.prps.has_key(prop_key):
      raise Exception ("'" + prop_key + "' property already present.")
    
    self.prps[prop_key] = value
    

  def has_prop(self, prop_key):

    return int(self.prps.has_key(prop_key))


  def remove_prop(self, prop_key):

    if self.prps.pop(prop_key,-9999) == -9999:
      raise Exception ("key + (" + prop_key + ") was not found, use 'ignore_missing' to ignore.")
    

  def clear_props(self):

    self.prps.clear()


  def prop_index(self, prop_key):
    
    keys = sorted(self.prps.keys())

    for prop_i in range(len(keys)):
      if prop_key == keys[prop_i]:
        return prop_i
      
    raise Exception ("Key '" + str(prop_key) + "' not found.")


#  def prop(self, index):
#
#    keys = sorted(self.prps.keys())
#
#    return self.prps[keys[index]]
  

  def prop_name(self, index):
    
    return sorted(self.prps.keys())[index]


  def prop_keys(self):

    return sorted(self.prps.keys())


  def prop(self, key):

    return self.prps[key]

 
    
    