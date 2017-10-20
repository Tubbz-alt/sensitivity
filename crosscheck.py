from ROOT import *
from facilities import HPGe
from detector import Component
from assay import Assay

niter = 10000
gROOT.SetBatch(1)

# Define component
spec_act=100e-6
eff=1.
mass=1.
comp = Component(spec_act,eff,mass)

# Define HPGe's
colors = [kRed, 8, kBlue, 6, 7, 9, 12, 28]
ge = [
  HPGe(1./86400.),
#  HPGe(2./86400.),
#  HPGe(5./86400.),
  HPGe(10./86400.),
#  HPGe(20./86400.),
#  HPGe(50./86400.),
  HPGe(100./86400.),
  HPGe(200./86400.),
]
livetime=14*86400

# Perform N HPGe countings
hMu = []
hSigma = []
hLimit = []
for ige in range(len(ge)):
  #hMu.append(   TH1D('hMu_%i'%ige,    'Mu_%i'%ige,    1000,0,200e-6))
  hMu.append(   TH1D('hMu_%i'%ige,    'Mu_%i'%ige,    int(500e-6*livetime),0,500e-6))
  hSigma.append(TH1D('hSigma_%i'%ige, 'Sigma_%i'%ige, 1000,0,100e-6))
  hLimit.append(TH1D('hLimit_%i'%ige, 'Limit_%i'%ige, 1000,0,200e-6))

for ige in range(len(ge)):
  hLimit[ige].SetLineColor(colors[ige])
  hMu[ige].SetLineColor(colors[ige])
  hSigma[ige].SetLineColor(colors[ige])
  hMu[ige].SetLineStyle(2)
  hSigma[ige].SetLineStyle(2)
  hLimit[ige].SetLineStyle(1)
  for i in range(niter):
    assay = Assay(ge[ige].count(comp.trueimp,1,livetime=livetime))
    if 'limit' in assay.params.keys():
      hLimit[ige].Fill(assay.params['limit'])
    else:
      hMu[ige].Fill(assay.params['mu'])
      hSigma[ige].Fill(assay.params['sigma'])
    #print assay

cCommon = TCanvas('cCommon','Common',1024,768)
cCommon.SetGridx()
cCommon.SetGridy()
cCommon.SetLogy()

fGaus = TF1('fGaus','gaus',0.08e-3,0.2e-3)
hMu[2].Fit(fGaus,'r')

for ige in range(len(ge)):
  hMu[ige].Draw('' if ige == 0 else 'same')
  #hSigma[ige].Draw('same') # if ige == 0 else 'same')
  hLimit[ige].Draw('same')
cCommon.SaveAs('Mu.png')

for ige in range(len(ge)):
  hSigma[ige].Draw('' if ige == 0 else 'same')
cCommon.SaveAs('Sigma.png')

for ige in range(len(ge)-1,-1,-1):
  hLimit[ige].Draw('' if ige == len(ge)-1 else 'same')
cCommon.SaveAs('Limit.png')

for ige in range(len(ge)):
  print ge[ige].bkgrate*86400, hMu[ige].GetEntries(), hLimit[ige].GetEntries()

for ige in range(len(ge)):
  print ge[ige].bkgrate*86400, hMu[ige].GetMean()/1e-6, hLimit[ige].GetMean()/1e-6

