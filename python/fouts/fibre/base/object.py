


class Object(object):
  
  def __init__(self, size=0, value=0.0):
  
    self.vector = list()
    self.prps = dict()
    self.extend_prps = dict()
    
    for i in range(size):
      self.vector.append((value,value,value))
    
  def size(self):
    
    return len(self.vector)  
  
  def num_props(self):
    
    return len(self.prps)
  
    
  def __str__(self):
    
    string = ''
    
    for key_val in sorted(self.prps.items()):
      string = string + key_val[0] + ': ' + str(key_val[1]) + '\n' 
    
    for triple in self.vector:
      string = string + '[ %s, %s, %s ]\n' % (triple[0], triple[1], triple[2])

    return string