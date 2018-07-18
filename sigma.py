import math,json,sys
import numpy as np
import pickle
import ROOT 

from detector import Detector, Component, define_detector
from facilities import HPGe
from assay import Assay
from toysens import calc_sens, calc_sens2

ROOT.gROOT.SetBatch(1)

def to_hl(x):
  return 1e6/136.*6.02e23*math.log(2)/(x/10.)/1e27

def optim_height(h):
  maxvalue = 0
  for i in range(1,h.GetNbinsX()):
    maxvalue = max(maxvalue,h.GetBinContent(i))
  return 0.9/maxvalue if maxvalue != 0 else 0
    

def find_chi2(det, livetime, nsenstoys, nrep):

  nparts = len(det.components)

  # Find chi2 curve
  results = {}
  methods = ['Truth']
  #methods += ['V:%.2f:%i' % (0.1*i,nparts) for i in range(30,-35,-5)]
  #methods += ['Central']
  #methods += ['TruncatedGaussian']
  methods += ['TruncatedGaussian']
  
  hULs = ROOT.TH1D('hULs','ULs',120,0,12)
  hSens = ROOT.TH1D('hSens','Sens',120,0,12)

  sens_lists = {}
  for method in methods:
    print(method)
    sens_list = []
    all_uls = []
    for j in range(nrep):
      if j % (nrep//10) == 0: print(j)
      #sens, uls = calc_sens2(det, method, livetime, nsenstoys * (10 if method[0] in ['T', 'V'] else 1), method == 'Truth')
      sens, uls = calc_sens2(det, method, livetime, 10000 if method == 'Truth' else nsenstoys, method == 'Truth')
      sens_list.append(to_hl(sens))
      #print 'iter',sens
      if method == 'TruncatedGaussian': all_uls += uls
    avgsens = float(np.mean(sens_list))
    #errsens = np.std(sens_list)/math.sqrt(nrep)
    errsens = float(np.std(sens_list))
    print('Sensitivity',method,'\t',avgsens,errsens)
    results[method] = (avgsens, errsens)
    sens_lists[method] = sens_list

    if method == 'TruncatedGaussian':
      for sens in sens_list:
        hSens.Fill(sens)
  
  lTruth = ROOT.TLine(results['Truth'][0],0,results['Truth'][0],10)
  lTruth.SetLineColor(ROOT.kBlack)
  lTruth.SetLineStyle(2)
  lTruth.SetLineWidth(2)
  lTruncatedGaussian = ROOT.TLine(results['TruncatedGaussian'][0],0,results['TruncatedGaussian'][0],10)
  lTruncatedGaussian.SetLineColor(ROOT.kBlue)
  lTruncatedGaussian.SetLineStyle(2)
  lTruncatedGaussian.SetLineWidth(2)
  
  cErrorBar = ROOT.TCanvas('cErrorBar','ErrorBar',1024,768)
  cErrorBar.SetGridx()
  cErrorBar.SetGridy()
  hSens.Draw()
  lTruth.Draw()
  lTruncatedGaussian.Draw()
  cErrorBar.SaveAs('ErrorBar.png')

  '''
  cChi2 = ROOT.TCanvas('cChi2','Chi2',1024,768)
  cChi2.SetGridx()
  cChi2.SetGridy()
  hULs.Scale(10.*optim_height(hULs))
  #hULs.GetXaxis().SetRangeUser(0,12)
  #hULs.GetYaxis().SetRangeUser(0,10)
  hULs.GetXaxis().SetTitle('Sensitivity [10^{27} y]')
  hULs.GetYaxis().SetTitle('#Delta#chi^{2}')
  hULs.Draw('hist')
  #gChi2.Draw('samel*')
  cChi2.Update()
  lTruth.Draw()
  lTruncatedGaussian.Draw()
  #lCentral.Draw()
  cChi2.SaveAs('Chi2.png')
  '''
  # == end temp

  #nsig_chi2 = math.sqrt(gChi2.Eval(results['Truth'][0]))
  #nsig_chi2 *= 1 if results['Truth'][0] > results['V:%.2f:%i' % (0,nparts)][0] else -1

  nsig_chi2 = (results['Truth'][0] - results['TruncatedGaussian'][0])/results['TruncatedGaussian'][1] if results['TruncatedGaussian'][1] != 0 else 0
  #print('nsig_chi2',nsig_chi2)

  #return {'det': det.data(), 'results': results}
  return {'det': det.data(), 'results': results, 'nsig_chi2': nsig_chi2, 'sens_lists': sens_lists}

if __name__ == '__main__':

  if len(sys.argv) < 8:
    print('See code for usage!')
    exit(0)

  # Parse arguments
  ntoys = int(sys.argv[1])
  true_lambda = float(sys.argv[2])
  ncomp = int(sys.argv[3])
  spec_act = float(sys.argv[4])
  livetime = float(sys.argv[5])*86400*365
  nsenstoys = int(sys.argv[6])
  nrep = int(sys.argv[7])
  #label = sys.argv[7]
  spy = 86400*365
  
  print(sys.argv)

  # Define det
  det = define_detector(true_lambda, ncomp, spec_act, livetime)

  # Find chi2 curve
  datasets = []
  for itoy in range(ntoys):
    # Perform assay
    ge = HPGe(10./86400.)
    for comp in det.components:
      comp.assay = Assay(ge.count(comp.trueimp,1,livetime=14*86400))
    nparts = len(det.components)

    #if itoy % (ntoys//10) == 0: print(itoy)
    dataset = find_chi2(det, livetime, nsenstoys, nrep)
    datasets.append(dataset)
    #print(dataset)
    results = dataset['results']
    nsig_chi2 = dataset['nsig_chi2']
    print('answer:',results['Truth'][0],results['Truth'][1],results['TruncatedGaussian'][0],results['TruncatedGaussian'][1], nsig_chi2)

  #pickle.dump(datasets,open('ds-%s.pickle' % label,'wb'))
  pickle.dump(datasets,open('ds.pickle','wb'),2)
 