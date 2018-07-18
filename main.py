import numpy as np
import math,sys

from detector import Detector, Component, define_detector
from assay import Assay
from facilities import HPGe
import toysens

from ROOT import *

# Units: 
#   impurity: Bq/kg
#   mass: kg
#   livetime: s
#   rate: Hz

def main(ntoys, true_lambda, ncomp, spec_act, livetime, method, nsenstoys):

  # Define detector and its parts
  #comps = []
  #mass = 1  #per comp
  usetruth = (method == 'Truth')

  print('Arguments:')
  print('ntoys',ntoys)
  print('true_lambda',true_lambda)
  print('ncomp',ncomp)
  print('livetime',livetime)
  print('method',method)
  print('nsenstoys',nsenstoys)
  #print('mass',mass)

  det = define_detector(true_lambda, ncomp, spec_act, livetime)
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
 
  main(ntoys, true_lambda, ncomp, spec_act, livetime, method, nsenstoys)
  print('cache_size',len(toysens.fc_cache))
