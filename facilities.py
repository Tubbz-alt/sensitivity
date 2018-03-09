import numpy as np
from math import sqrt
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

if __name__ == '__main__':
  main()
