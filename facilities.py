import numpy as np
from math import sqrt
#import icpms

def main():
  ge = HPGe(0.01)  

  for i in range(100):
    print(ge.count(0.13e-3,1,livetime=86400*14))

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
    if sig > 1.64 * sigerr:
      ans = {'mu': sig/norm, 'sigma': sigerr/norm}
    else:
      ans = {'limit': 1.64 * sigerr/norm}
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
