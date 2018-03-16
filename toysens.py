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

    #print('pre-fc',counts,true_counts)
    fcul = feldman(counts,true_counts)
    #print('fcul',fcul)

    #upperlimits.append(fcul)
    meanul += fcul

  meanul /= ntoys
  #return np.percentile(upperlimits,50), upperlimits
  #return np.mean(upperlimits), upperlimits
  return meanul, upperlimits
  
