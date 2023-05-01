import ROOT,PlotUtils,UnfoldUtils
ROOT.TH1.AddDirectory(False)

def MnvUnfoldWrapper(ifile_name,ite):
  ifile = ROOT.TFile.Open(ifile_name)
  response=ifile.Get("migration")
  reco=ifile.Get("reco")
  truth=ifile.Get("reco_truth")
  hlin_measured = ifile.Get("data")
  unfold_old = ifile.Get("data_truth")
  mnvunfold=UnfoldUtils.MnvUnfold()
  h_unfolded =truth.Clone()
  covmx = ROOT.TMatrixD(500,500)
  mnvunfold.UnfoldHistoWithFakes(h_unfolded,covmx,response,response.ProjectionX("fuck",0,-1),ROOT.nullptr
,ROOT.nullptr
,ROOT.nullptr,ite,False,False)
  print((hlin_measured.Integral(),h_unfolded.Integral(0,-1),response.ProjectionY().Integral(0,-1)))
  for i in range(18):
    print((h_unfolded.GetBinContent(i),unfold_old.GetBinContent(i)))
  print ("done")
  ifile.Close()

MnvUnfoldWrapper("fin.root",1)
MnvUnfoldWrapper("fin.root",2)
#MnvUnfoldWrapper("fin.root",3)
#MnvUnfoldWrapper("fin.root",4)
#MnvUnfoldWrapper("fin.root",5)
