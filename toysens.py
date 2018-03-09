import numpy as np
#from TimeProfiler import *
import math,sys

from detector import Detector, Component
from assay import Assay
from facilities import HPGe

from ROOT import *

# Units: 
#   impurity: Bq/kg
#   mass: kg
#   livetime: s
#   rate: Hz

gROOT.SetBatch(1)
gStyle.SetOptStat(0)
FeldmanCousins = TFeldmanCousins()
fc_cache = {}

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

  if ncomp == 0: # Realistic detector
    spec_acts = [5e-6, 20e-6, 100e-6]

    nregions = len(spec_acts)
    nparts = 10    # Number of parts 
    budgetPerPart = true_lambda/nregions/nparts
   
    effs = [1.*budgetPerPart/s/livetime/mass for s in spec_acts]

    print('spec_acts',spec_acts)
    print('effs',effs)

    for s,e in zip(spec_acts,effs):
      for i in range(nparts):
        comps.append(Component(s,e,mass))

  else: # Detector with identical parts
    eff = 1.*true_lambda/ncomp/spec_act/livetime/mass
  
    print('spec_act',spec_act)
    print('eff',eff)
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
      true_sens, true_uls = calc_sens(det, method, livetime, nsenstoys, True)
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
      print(det)
  
      sens, uls = calc_sens(det, method, livetime, nsenstoys, False)

      print('Sensitivity', sens, nlimits)
      #print('Upperlimits', uls)
      print

#@profile
def print_histo(l):
  histo = TH1D('histo','histo',50,0,50)
  for v in l:
    histo.Fill(v)
  canvas = TCanvas('canvas','canvas',1024,768)
  canvas.SetGridx()
  canvas.SetGridy()
  histo.Draw()
  canvas.SaveAs('histo.png')

#@profile
def feldman(counts,true_counts):

  cache_key = str([counts,"%.2f"%true_counts])
  if cache_key in fc_cache.keys():
    return fc_cache[cache_key]

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

  fc_cache[cache_key] = fcul
  return fcul

#@profile
def calc_sens(det, method, livetime, ntoys, truth):
  
  upperlimits = []
  meanul = 0

  ult_truth = det.truth()*livetime

  for i in range(ntoys):
    if i % (ntoys//10) == 0: print(i )

    true_counts = det.throw(method)*livetime if not truth else ult_truth
    counts = np.random.poisson(true_counts)

    fcul = feldman(counts,true_counts)
    print('fcul',fcul)

    #upperlimits.append(fcul)
    meanul += fcul

  meanul /= ntoys
  #return np.percentile(upperlimits,50), upperlimits
  #return np.mean(upperlimits), upperlimits
  return meanul, upperlimits
  
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
  print('cache_size',len(fc_cache))
