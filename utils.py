def GetCL(h,cl=0.9):
  target = 1.*h.Integral()*cl
  area = 0 
  for ibin in range(h.GetNbinsX()):
    area += h.GetBinContent(ibin)
    if area >= target:
      x = (1 - (area - target)/h.GetBinContent(ibin)) if h.GetBinContent(ibin) > 0 else 0
      return h.GetBinLowEdge(ibin) + x*h.GetBinWidth(ibin)
      #return h.GetBinCenter(ibin)

