import ROOT
import PlotUtils
from config.AnalysisConfig import AnalysisConfig

ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.15)

playlist = AnalysisConfig.playlist
mc_file = "/minerva/data/users/ajball/nu_e/kin_dist_mc"+str(playlist)+"_collab1_fspline.root"
rootFile = ROOT.TFile.Open(mc_file,"READ")
t_numerator = rootFile.tDef1_diff.Clone()
t_numerator.Print("all")
t_denominator = rootFile.true_t.Clone()
t_denominator.Print("all")
pion_numerator = rootFile.PiZeroE_diffractive.Clone()
pion_denominator = rootFile.true_pizeroE.Clone()
t_numerator.Divide(t_numerator,t_denominator,1.0,1.0)
t_efficiency = t_numerator.Clone()
t_efficiency.GetXaxis().SetTitle("|t| (GeV/c)^{2}")
t_efficiency.GetYaxis().SetTitle("Efficiency")
t_efficiency.SetTitle("t Efficiency")
t_efficiency.SetName("t_efficiency")
t_efficiency.GetCVHistoWithError().Print("all")
pion_numerator.Divide(pion_numerator,pion_denominator,1.0,1.0)
pion_efficiency = pion_numerator.Clone()
pion_efficiency.GetYaxis().SetTitle("Efficiency")
pion_efficiency.SetTitle("Pion Energy Efficiency")
pion_efficiency.SetName("pion_efficiency")
pion_efficiency.GetCVHistoWithError().Print("all")
efficiencyFile = ROOT.TFile.Open("/minerva/data/users/ajball/nu_e/"+str(playlist)+"_calculatedEfficiency.root","RECREATE")
t_efficiency.GetXaxis().SetTitleFont(62)
t_efficiency.GetXaxis().SetTitleSize(0.045)
t_efficiency.GetYaxis().SetTitleFont(62)
t_efficiency.GetYaxis().SetTitleSize(0.045)
t_efficiency.SetTitleFont(62)
t_efficiency.SetLabelFont(62,"xyz")
t_efficiency.SetMarkerStyle(21)
t_efficiency.SetMarkerSize(2)
t_efficiency.SetLineWidth(4)
t_efficiency.SetStats(0)
t_efficiency.Write()
pion_efficiency.GetXaxis().SetTitleFont(62)
pion_efficiency.GetXaxis().SetTitleSize(0.045)
pion_efficiency.GetYaxis().SetTitleFont(62)
pion_efficiency.GetYaxis().SetTitleSize(0.045)
pion_efficiency.SetTitleFont(62)
pion_efficiency.SetLabelFont(62,"xyz")
pion_efficiency.SetMarkerStyle(21)
pion_efficiency.SetMarkerSize(2)
pion_efficiency.SetLineWidth(4)
pion_efficiency.SetStats(0)
pion_efficiency.Write()
t_canvas = ROOT.TCanvas("c1","c1",1600,1200)
t_eff_err = t_efficiency.GetCVHistoWithError()
t_eff_err.Draw("MIN0")
t_canvas.Print("/minerva/data/users/ajball/nu_e/plot/t_efficiency.png","png")
pion_canvas = ROOT.TCanvas("c2","c2",1600,1200)
pion_eff_err = pion_efficiency.GetCVHistoWithError()
pion_eff_err.Draw("MIN0")
pion_canvas.Print("/minerva/data/users/ajball/nu_e/plot/pion_efficiency.png","png")
efficiencyFile.Close()
rootFile.Close()