import numpy as np
from math import sqrt
#import icpms
from ROOT import *

def main():
  gROOT.SetBatch(1)
  hMu = TH1F('hMu','Mu',1000,0,200e-6)

  ge = HPGe(10./86400.)  
  for i in range(100000):
    res = ge.count(100e-6,1,livetime=86400*365)
    if 'limit' not in res.keys():
      hMu.Fill(res['mu'])

  cMu = TCanvas('cMu','Mu',1024,768)
  hMu.Draw()
  cMu.SaveAs('Mu.png')

  #ic = ICPMS(200,10,1,10)
  #print(ic.count(1e04,1000))

class HPGe:
  def __init__(self,bkgrate):
    self.bkgrate = bkgrate

  def count(self, trueimp, mass, **options):

    livetime=options['livetime']

    srcrun_truecount = (trueimp+self.bkgrate)*mass*livetime
    bkgrun_truecount =  self.bkgrate*mass*livetime
    #print(srcrun_truecount, bkgrun_truecount)

    srcrun_count = np.random.poisson(srcrun_truecount) 
    bkgrun_count = np.random.poisson(bkgrun_truecount) 
    #print(srcrun_count, bkgrun_count)

    return self.report(srcrun_count, bkgrun_count, 1.*mass*livetime)

  def report(self,s,b,norm):
    sig = s-b
    sigerr = sqrt(s+b)

    # signal is above discovery threshold.
    orig = {'mu': sig/norm, 'sigma': sigerr/norm, 'src': s, 'bkg': b, 'norm': norm}
    if sig > 1.64 * sigerr:
      ans = {'mu': sig/norm, 'sigma': sigerr/norm, 'original': orig}
    elif sig > 0:
      ans = {'limit': (1.64 * sigerr + sig)/norm, 'original': orig}
    else:
      ans = {'limit': 1.64 * sigerr/norm, 'original': orig}
    return ans

'''
class ICPMS:
  def __init__(self,prob,blank,blanksigma,tracer):
    self.prob = prob
    self.blank = blank
    self.blanksigma = blanksigma
    self.tracer = tracer

  def count(self, trueimp, mass, **options):
    sample_mat = icpms.Material("sample",{"U238": [trueimp, 0]})
    blank_mat =  icpms.Material("blank", {"U238": [self.blank, self.blanksigma]})
    tracer_mat = icpms.Material("tracer",{"U233": [self.tracer, 0]})
    ans, fullresult, lod = icpms.assay(sample_mat, blank_mat, tracer_mat, mass=mass, icpms_prob=self.prob)
    return ans
'''

if __name__ == '__main__':
  main()
