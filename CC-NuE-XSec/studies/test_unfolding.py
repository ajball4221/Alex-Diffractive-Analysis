import ROOT
import UnfoldUtils

def MySubTH1(h1,h2):
    h3 = h1.Clone()
    h3.Add(h2,-1)
    s = h3.GetSize()
    for i in range(s):
        h3.SetBinError(i,sqrt(h1.GetBinContent(i)+h2.GetBinContent(i)))
    return h3


MnvUnfold = UnfoldUtils.MnvUnfold()
