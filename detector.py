from ROOT import * 
import math

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
      #print('throw c',c)
      s += c.throw(method)
    return s #sum([c.throw() for c in self.components])

  def truth(self):
    s = 0 
    for c in self.components:
      s += c.truth()
    return s #sum([c.throw() for c in self.components])

  def data(self):
    return [{'mass': comp.mass, 
             'eff': comp.efficiency,
             'trueimp': comp.trueimp, 
             'assay': comp.assay.params} for comp in self.components]

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

FeldmanCousins = TFeldmanCousins()
def feldman(counts,true_counts):
  FeldmanCousins.SetMuStep(0.005)
  FeldmanCousins.SetMuMax(max(5.*counts,10))
  fcll = -999
  while True:
    fcul = FeldmanCousins.CalculateUpperLimit(counts,true_counts)
    fcll = FeldmanCousins.GetLowerLimit()
    
    if fcul == 0:
      #print('mumax', FeldmanCousins.GetMuMax())
      FeldmanCousins.SetMuMax(FeldmanCousins.GetMuMax()*2)
      continue

    if fcll == -999: 
      #print('mustep', FeldmanCousins.GetMuStep())
      FeldmanCousins.SetMuStep(FeldmanCousins.GetMuStep()/2.)
      continue

    break

  if fcll == -999 or fcul == 0:
    print('!!!',counts,true_counts, fcll, fcul)

  return fcul

def define_detector(true_lambda, ncomp, spec_act, livetime):

  comps = []
  mass = 1

  if ncomp == 0: # Realistic detector type 0
    spec_acts = [5e-6, 20e-6, 100e-6]

    nregions = len(spec_acts)
    nparts = 10    # Number of parts 
    budgetPerPart = true_lambda/nregions/nparts
   
    effs = [1.*budgetPerPart/s/livetime/mass for s in spec_acts]

    #print('spec_acts',spec_acts)
    #print('effs',effs)

    for s,e in zip(spec_acts,effs):
      for i in range(nparts):
        comps.append(Component(s,e,mass))
  elif ncomp < 0:
    # Power low rate. Linear spec act distrib.
    # Type 1: small x with large eff
    # Type 2: large x with large eff

    nparts = 10   # Number of parts
    #spec_act = 100e-6

    b = 0.5 #math.log(1.7)
    
    allocations = [math.exp(-i*b) for i in range(nparts)]
    allocations = [a*true_lambda/sum(allocations) for a in allocations]
    
    if ncomp == -1: 
      spec_acts = [(i+1)*spec_act/nparts for i in range(nparts)]
    elif ncomp == -2: 
      spec_acts = [(nparts-i)*spec_act/nparts for i in range(nparts)]
    elif ncomp == -3: # flat
      spec_acts = [(nparts+1)/2*spec_act/nparts for i in range(nparts)]
    
    effs = [a/x/livetime/mass for a,x in zip(allocations,spec_acts)]
    
    print('spec_acts',spec_acts)
    print('effs',effs)
    print('c',[x*e*365*86400 for x,e in zip(spec_acts,effs)])

    #eff_0 = true_lambda/(spec_act/nparts)/(sum([math.exp(-i*b) for i in range(nparts)]))/livetime
    #print('x_0',spec_act/nparts)
    #print('eff_0',eff_0)

    for s,e in zip(spec_acts,effs):
      comps.append(Component(s,e,mass))

  else: # Detector with identical parts
    eff = 1.*true_lambda/ncomp/spec_act/livetime/mass
  
    #print('spec_act',spec_act)
    #print('eff',eff)
    # ===

    for ic in range(ncomp):
      comps.append(Component(spec_act,eff,mass))
  det = Detector(comps)
  #print('det',det)

  return det

if __name__ == '__main__':
  from assay import Assay
  import numpy as np

  assay = Assay({'limit':16.4})
  method = 'Plateau'
  #method = 'Delta'
  #method = 'Uniform'
  #method = 'Gaussian'

  comp = Component(1,1,1)
  comp.assay = assay
  true_counts = comp.throw(method)

  histo = TH1D('histo','histo',50,0,50)
  for i in range(100000):
    counts = np.random.poisson(true_counts)
    histo.Fill(counts)
    #fcul = feldman(counts,true_counts)
    #histo.Fill(fcul)

  canvas = TCanvas('canvas','canvas',1024,768)
  canvas.SetGridx()
  canvas.SetGridy()
  canvas.SetLogy()

  ''' 
  fGaus = TF1('fGaus','[0]/sqrt(2*pi)/[2] * exp(-((x-[1])**2)/[2]/[2]/2)',0,2) 
  #fGaus = TF1('fGaus','gaus',0,2) 
  fGaus.SetParameter(0,1)
  fGaus.SetParameter(1,0)
  fGaus.SetParameter(2,1)
  histo.Fit(fGaus,'r')
  '''
  histo.Draw()
  canvas.SaveAs('test.png')
