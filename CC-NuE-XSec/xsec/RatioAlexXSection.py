import ROOT
import PlotUtils
from config.AnalysisConfig import AnalysisConfig
from tools import Utilities, PlotTools
from tools.PlotLibrary import PLOT_SETTINGS
from config.SystematicsConfig import USE_NUE_CONSTRAINT,CONSOLIDATED_ERROR_GROUPS,DETAILED_ERROR_GROUPS,AnaNuPDG
from config import DrawingConfig
from config.CutConfig import NEUTRINO_ENERGY_RANGE,FIDUCIAL_Z_RANGE
from functools import partial

ROOT.TH1.AddDirectory(False)

USE_BIGNUE = True
threshold = 100 if USE_BIGNUE else 1
TARGET_UTILS = PlotUtils.TargetUtils.Get()
warping_errorband = ["fsi_weight","SuSA_Valencia_Weight","MK_model","LowQ2Pi_Joint","LowQ2Pi_NUPI0","LowQ2Pi_None"]
#warping_errorband = ["SuSA_Valencia_Weight"]
FLUX="minervame1d1m1nweightedave"
#FLUX="minervame1d"

def GetXSectionHistogram(unfolded,efficiency,is_mc):
    #divide by efficiency
    efficiency.AddMissingErrorBandsAndFillWithCV(unfolded)
    unfolded.Divide(unfolded,efficiency)
    #divide by flux
    DivideFlux(unfolded,is_mc)
    #divide by N nucleaon
    Nnucleon = TARGET_UTILS.GetTrackerNNucleons(FIDUCIAL_Z_RANGE[0],FIDUCIAL_Z_RANGE[1],is_mc)
    h_nucleon = GetNnucleonError(unfolded,Nnucleon)
    unfolded.Divide(unfolded,h_nucleon)
    unfolded = BinNormalize(unfolded)
    return unfolded

def GetNnucleonError(hist,ntargets):
    hist_target = hist.Clone("number_of_targets")
    hist_target.ClearAllErrorBands()
    hist_target.Reset()
    errband_name = "Target_Mass_CH"
    band_err = 0.014
    hist_target.AddVertErrorBand(errband_name,2)

    for i in range(hist_target.GetSize()):
        hist_target.SetBinContent(i,ntargets)
        hist_target.SetBinError(i,0)
        hist_target.GetVertErrorBand(errband_name).SetBinContent(i,ntargets)
        hist_target.GetVertErrorBand(errband_name).GetHist(0).SetBinContent(i,ntargets*(1-band_err))
        hist_target.GetVertErrorBand(errband_name).GetHist(1).SetBinContent(i,ntargets*(1+band_err))
    hist_target.AddMissingErrorBandsAndFillWithCV(hist)
    print("Target Normalization: {:.4E},{:.4E}".format(ntargets,ntargets*0.014))
    return hist_target

def DivideFlux(unfolded,is_mc):
    frw= PlotUtils.flux_reweighter(FLUX,AnaNuPDG,USE_NUE_CONSTRAINT) #playlist is dummy for now
    flux = frw.GetIntegratedFluxReweighted(AnaNuPDG,unfolded,NEUTRINO_ENERGY_RANGE[0],NEUTRINO_ENERGY_RANGE[1],False)
    #flux.PopVertErrorBand("Flux_BeamFocus")
    #flux.PopVertErrorBand("ppfx1_Total")
    flux.Scale(1e-4*(mc_pot if is_mc else data_pot)) #change unit to nu/cm^2
    print ("Flux Normalization: {:.4E},{:.4E}".format(flux.GetBinContent(1,1),flux.GetTotalError(False).GetBinContent(1,1)))
    unfolded.Divide(unfolded,flux)

def getPOTFromFile(filename):
    metatree = ROOT.TChain("Meta")
    if metatree.Add(filename,-1):
        return ROOT.PlotUtils.POTCounter().getPOTfromTChain(metatree)
    else:
        return None

def BinNormalize(hist):
    bin_sizes = hist.Clone()
    bin_sizes.ClearAllErrorBands()
    bin_sizes.Reset()
    
    xaxis = hist.GetXaxis()
    for i in range(hist.GetNbinsX()):
        bin_sizes.SetBinContent(i+1,xaxis.GetBinLowEdge(i+2)-xaxis.GetBinLowEdge(i+1))
    bin_sizes.AddMissingErrorBandsAndFillWithCV(hist)
    bin_sizes.Print("all")
    hist.Divide(hist,bin_sizes)
    return hist


playlist = "Alex_final_nx"
unfoldedFile = ROOT.TFile.Open("/minerva/data/users/ajball/nu_e/"+str(playlist)+"_unfolded.root","READ")
t_fitted_hist = unfoldedFile.t_fitted.Clone()
t_matrix_hist = unfoldedFile.t_recovered.Clone()
#pion_fitted_hist = unfoldedFile.pionE_fitted.Clone()
pion_matrix_hist = unfoldedFile.pionE_recovered.Clone()

mcTruthFile = ROOT.TFile.Open("/minerva/data/users/ajball/nu_e/kin_dist_mcAlex_tcut_nx_collab1_fspline.root","READ")
t_sim_hist = mcTruthFile.tDef1_diff
pion_sim_hist = mcTruthFile.PiZeroE_diffractive

efficiencyFile = ROOT.TFile.Open("/minerva/data/users/ajball/nu_e/"+str(playlist)+"_calculatedEfficiency.root","READ")
t_efficiency = efficiencyFile.t_efficiency.Clone()
pion_efficiency = efficiencyFile.pion_efficiency.Clone()
unfoldedFile.Close()
efficiencyFile.Close()

data_pot = getPOTFromFile("/minerva/data/users/ajball/nu_e/kin_dist_dataAlex_tcut_nx_collab1_fspline.root")
mc_pot = getPOTFromFile("/minerva/data/users/ajball/nu_e/kin_dist_mcAlex_tcutdfronly_nx_collab1_fspline.root")

t_fitted_xsec = GetXSectionHistogram(t_fitted_hist,t_efficiency,False)
t_matrix_xsec = GetXSectionHistogram(t_matrix_hist,t_efficiency,False)
t_sim_xsec = GetXSectionHistogram(t_sim_hist,t_efficiency,True)
#pion_fitted_xsec = GetXSectionHistogram(pion_fitted_hist,pion_efficiency,False)
pion_matrix_xsec = GetXSectionHistogram(pion_matrix_hist,pion_efficiency,False)
pion_sim_xsec = GetXSectionHistogram(pion_sim_hist,pion_efficiency,True)

t_fitted_xsec.Divide(t_fitted_xsec,t_sim_xsec)
t_matrix_xsec.Divide(t_matrix_xsec,t_sim_xsec)
t_sim_xsec.Divide(t_sim_xsec,t_sim_xsec)
pion_matrix_xsec.Divide(pion_matrix_xsec,pion_sim_xsec)
pion_sim_xsec.Divide(pion_sim_xsec,pion_sim_xsec)

outPath = "/minerva/data/users/ajball/nu_e/"+str(playlist)+"_xsec.root"
outFile = ROOT.TFile.Open(outPath,"RECREATE")

t_canvas = ROOT.TCanvas("c1","c1",1600,1200)
t_fitted_err = t_fitted_xsec.GetCVHistoWithError()
t_fitted_err.SetMarkerStyle(20)
t_fitted_err.SetMarkerColor(2)
t_fitted_err.SetLineColor(46)
t_fitted_err.SetMaximum(max([t_matrix_xsec.GetMaximum(),t_fitted_xsec.GetMaximum(),t_sim_xsec.GetMaximum()])*1.25)
t_fitted_err.SetMinimum(min([t_matrix_xsec.GetMinimum(),t_fitted_xsec.GetMinimum(),t_sim_xsec.GetMinimum()])/1.25)
t_fitted_err.GetXaxis().SetTitle("|t| (GeV/c)^2")
t_fitted_err.GetYaxis().SetTitle("dSigma/d|t| (cm^2/(GeV/c)^2/nucleon)")
t_fitted_err.SetTitle("|t| Cross Section from Exponential Fit")

t_recovered_err = t_matrix_xsec.GetCVHistoWithError()
t_recovered_err.SetMarkerStyle(21)
t_recovered_err.SetMarkerColor(4)
t_recovered_err.SetLineColor(9)
t_recovered_err.GetXaxis().SetTitle("|t| (GeV/c)^2")
t_recovered_err.GetYaxis().SetTitle("dSigma/d|t| (cm^2/(GeV/c)^2/nucleon)")
t_recovered_err.SetTitle("|t| Cross Section from Inverse Matrix")

t_sim_err = t_sim_xsec.GetCVHistoWithError()
t_sim_err.SetMarkerStyle(21)
t_sim_err.SetMarkerColor(8)
t_sim_err.SetLineColor(8)
t_sim_err.GetXaxis().SetTitle("|t| (GeV/c)^2")
t_sim_err.GetYaxis().SetTitle("dSigma/d|t| (cm^2/(GeV/c)^2/nucleon)")
t_sim_err.SetTitle("|t| Cross Section from Simulation")

t_fitted_err.Draw("MIN0 E1")
t_recovered_err.Draw("SAME MIN0 E1")
#t_sim_err.Draw("SAME MIN0 E1")
ROOT.gPad.BuildLegend()
t_canvas.Print("datafolder/xsecplots/t_xsec_ratio.png","png")
t_canvas.SetLogy(1)
t_canvas.Print("datafolder/xsecplots/t_xsec_ratio_log.png","png")

pion_canvas = ROOT.TCanvas("c2","c2",1600,1200)
'''pion_fitted_err = pion_fitted_xsec.GetCVHistoWithError()
pion_fitted_err.SetMarkerStyle(20)
pion_fitted_err.SetMarkerColor(2)
pion_fitted_err.SetLineColor(46)
pion_fitted_err.SetMaximum(max([pion_matrix_xsec.GetMaximum(),pion_fitted_xsec.GetMaximum()])*1.25)
pion_fitted_err.GetXaxis().SetTitle("Pion Energy (GeV)")
pion_fitted_err.GetYaxis().SetTitle("dSigma/dpionE (cm^2/GeV/nucleon)")
pion_fitted_err.SetTitle("PionE Cross Section from Exponential Fit")'''

pion_recovered_err = pion_matrix_xsec.GetCVHistoWithError()
pion_recovered_err.SetMarkerStyle(21)
pion_recovered_err.SetMarkerColor(4)
pion_recovered_err.SetLineColor(9)
pion_recovered_err.SetMaximum(max([pion_matrix_xsec.GetMaximum(),pion_sim_xsec.GetMaximum()])*1.7)
pion_recovered_err.SetMinimum(min([pion_sim_xsec.GetMinimum()])/1.25)
pion_recovered_err.GetXaxis().SetTitle("Pion Energy (GeV)")
pion_recovered_err.GetYaxis().SetTitle("dSigma/dpionE (cm^2/GeV/nucleon)")
pion_recovered_err.SetTitle("PionE Cross Section from Inverse Matrix")

pion_sim_err = pion_sim_xsec.GetCVHistoWithError()
pion_sim_err.SetMarkerStyle(21)
pion_sim_err.SetMarkerColor(8)
pion_sim_err.SetLineColor(8)
pion_sim_err.GetXaxis().SetTitle("Pion Energy (GeV)")
pion_sim_err.GetYaxis().SetTitle("dSigma/dpionE (cm^2/GeV/nucleon)")
pion_sim_err.SetTitle("PionE Cross Section from Simulation")

#pion_fitted_err.Draw("MIN0 E1")
pion_recovered_err.Draw("MIN0 E1")
#pion_sim_err.Draw("SAME MIN0 E1")
#ROOT.gPad.BuildLegend()
ROOT.gPad.BuildLegend(0.45,0.75,0.85,0.85)
pion_canvas.Print("datafolder/xsecplots/pion_xsec_ratio.png","png")
pion_canvas.SetLogy(1)
pion_canvas.Print("datafolder/xsecplots/pion_xsec_ratio_log.png","png")

t_fitted_xsec.Write()
t_matrix_xsec.Write()
#pion_fitted_xsec.Write()
pion_matrix_xsec.Write()
t_sim_xsec.Write()
pion_sim_xsec.Write()
outFile.Close()
