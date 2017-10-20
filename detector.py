class Detector:
  def __init__(self,l=[]):
    self.components = l
  def __repr__(self):
    return str(self.components)

  def add(self,c):
    if isinstance(c,list):
      for item in c:
        self.add(item)
    else:
      self.components.append(c)

  def throw(self,method):
    s = 0 
    for c in self.components:
      s += c.throw(method)
    return s #sum([c.throw() for c in self.components])

  def truth(self):
    s = 0 
    for c in self.components:
      s += c.truth()
    return s #sum([c.throw() for c in self.components])
    

class Component:
  def __init__(self,trueimp,eff,mass):
    self.trueimp = trueimp
    self.assay = None
    self.efficiency = eff
    self.mass = mass
  def __repr__(self):
    return str([self.mass,self.efficiency,self.trueimp,self.assay])
  def throw(self,method):
    return self.assay.throw(method) * self.efficiency * self.mass
  def truth(self):
    return self.trueimp * self.efficiency * self.mass

