import numpy as np
#from TimeProfiler import *
import math,sys

from detector import Detector, Component
from assay import Assay
from facilities import HPGe
import toysens
#from toysens import calc_sens

from ROOT import *

# Units: 
#   impurity: Bq/kg
#   mass: kg
#   livetime: s
#   rate: Hz

#@profile
def main(ntoys, true_lambda, ncomp, spec_act, livetime, method, nsenstoys):

  # Define detector and its parts
  comps = []
  mass = 1  #per comp
  usetruth = (method == 'Truth')

  print('Arguments:')
  print('ntoys',ntoys)
  print('true_lambda',true_lambda)
  print('ncomp',ncomp)
  print('livetime',livetime)
  print('method',method)
  print('nsenstoys',nsenstoys)
  print('mass',mass)

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
  elif ncomp == -1 or ncomp == -2:  # Realistic detector types 1 and 2
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
    
    effs = [a/x/livetime/mass for a,x in zip(allocations,spec_acts)]
    
    #print('spec_acts',spec_acts)
    #print('effs',effs)

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
  print('det',det)

  # Assay settings
  ge = HPGe(10./86400.)
  
  if usetruth:
    # True sensitivity -- Perfect knowledge of impurities
    ult_truth = det.truth()*livetime
    print('Ultimate true counts', ult_truth)
    for i in range(ntoys):
      true_sens, true_uls = toysens.calc_sens(det, method, livetime, nsenstoys, True)
      print('Sensitivity', true_sens)
      #print('Upperlimits', true_uls)
      #print('RMS', np.std(true_uls))
      #print_histo(true_uls)
      print
  
  else:
  # Sensitivity in real life -- Knowledge of impurities from imperfect assay
    for i in range(ntoys):
  
      # Assay campaign
      nlimits = 0
      for comp in det.components:
        comp.assay = Assay(ge.count(comp.trueimp,1,livetime=14*86400))
        if 'limit' in comp.assay.params.keys():
          nlimits += 1
      #print(det)
  
      sens, uls = toysens.calc_sens(det, method, livetime, nsenstoys, False)

      print('Sensitivity', sens, nlimits)
      #print('Upperlimits', uls)
      print

if __name__ == '__main__':
  if len(sys.argv) < 8: 
    print('See code for usage!')
    exit(0)

  ntoys = int(sys.argv[1])
  true_lambda = float(sys.argv[2])
  ncomp = int(sys.argv[3])
  spec_act = float(sys.argv[4])
  livetime = float(sys.argv[5])*86400*365
  method = sys.argv[6]
  nsenstoys = int(sys.argv[7])
  #usetruth = int(sys.argv[8])
 
  main(ntoys, true_lambda, ncomp, spec_act, livetime, method, nsenstoys)
  #print_prof_data()
  print('cache_size',len(toysens.fc_cache))
