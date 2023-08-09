import ROOT
import PlotUtils
from config.AnalysisConfig import AnalysisConfig
from tools import Utilities, PlotTools
from tools.PlotLibrary import PLOT_SETTINGS
from config.SystematicsConfig import USE_NUE_CONSTRAINT,CONSOLIDATED_ERROR_GROUPS,DETAILED_ERROR_GROUPS,AnaNuPDG
from config import DrawingConfig
from config.CutConfig import NEUTRINO_ENERGY_RANGE,FIDUCIAL_Z_RANGE,FIDUCIAL_APOTHEM
from functools import partial
import numpy as np

ROOT.TH1.AddDirectory(False)

USE_BIGNUE = True
threshold = 100 if USE_BIGNUE else 1
TARGET_UTILS = PlotUtils.TargetUtils.Get()
warping_errorband = ["fsi_weight","SuSA_Valencia_Weight","MK_model","LowQ2Pi_Joint","LowQ2Pi_NUPI0","LowQ2Pi_None"]
#FLUX="minervame1d1m1nweightedave"
FLUX="minervame1d" # Switched to this one for closure

def GetXSectionHistogram(unfolded,efficiency,is_mc):
    # Got rid of efficiency because this should be all truth not reco
    #divide by flux
    DivideFlux(unfolded,is_mc)
    #divide by N Hydrogen
    nplanes = 2*(80-27+1) # Fiducial Volume -> Modules 27-80
    mass = TARGET_UTILS.GetTrackerMass(nplanes,is_mc,FIDUCIAL_APOTHEM)
    massFrac = TARGET_UTILS.GetTrackerElementMassFraction(1,is_mc)
    targets = mass*massFrac*TARGET_UTILS.GetTrackerElementAtomsPerGram(1)*1 # That last *1 is because there's 1 proton/atom for H
    h_nucleon = GetNnucleonError(unfolded,targets)
    unfolded.Divide(unfolded,h_nucleon)
    #Adding Bin Normalization
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
    #bin_sizes.Print("all")
    hist.Divide(hist,bin_sizes)
    return hist

def Integral(hist,name):
    integral = 0
    errors = np.zeros(hist.GetNbinsX())
    for i in range(hist.GetNbinsX()):
        integral += hist.GetBinContent(i+1)*(hist.GetXaxis().GetBinLowEdge(i+2)-hist.GetXaxis().GetBinLowEdge(i+1))
        errors[i] = hist.GetBinError(i+1)*(hist.GetXaxis().GetBinLowEdge(i+2)-hist.GetXaxis().GetBinLowEdge(i+1))
    error = np.sqrt(np.sum(np.square(errors)))
    print("The integral of the",name,"cross section is",integral,"and the error is",error)

playlist = AnalysisConfig.playlist
mcTruthFile = ROOT.TFile.Open("/minerva/data/users/ajball/nu_e/kin_dist_mc"+str(playlist)+"_collab1_fspline.root","READ")
t_sim_hist = mcTruthFile.true_t
pion_sim_hist = mcTruthFile.true_pizeroE

#mc_pot = getPOTFromFile("/minerva/data/users/ajball/nu_e/kin_dist_mcAlex_dfronlyclosure_nx_collab1_fspline.root")
mc_pot = getPOTFromFile("/minerva/data/users/ajball/nu_e/kin_dist_mc"+str(playlist)+"_collab1_fspline.root")

t_evtRate = t_sim_hist.Clone()
pion_evtRate = pion_sim_hist.Clone()
mytXSec = GetXSectionHistogram(t_sim_hist,None,True) # the None is where efficiency would go, I don't feel like rewriting it
myPionXSec = GetXSectionHistogram(pion_sim_hist,None,True)

GENIEXSecFile = ROOT.TFile.Open("/minerva/app/users/ajball/ncdiff/Alex_GENIEClosure_xsec.root")

print("Event Rate Ratios (For Debugging):")
GENIEPionEvtRate = GENIEXSecFile.cc_pionE_xsec_evRate.Clone()
pion_evtRate.Divide(pion_evtRate,GENIEPionEvtRate)
pion_evtRate.Print("all")

GENIEtEvtRate = GENIEXSecFile.cc_TDef1_xsec_evRate.Clone()
tEventRatio = t_evtRate/GENIEtEvtRate
tEventRatio.Print("all")

print("Cross Section Ratios (The Closure Test):")
GENIEPionXSec = GENIEXSecFile.cc_pionE_xsec.Clone()
pionClosure = myPionXSec/GENIEPionXSec
pionClosure.Print("all")

GENIEtXSec = GENIEXSecFile.cc_TDef1_xsec.Clone()
tClosure = mytXSec/GENIEtXSec
tClosure.Print("all")

mcTruthFile.Close()
GENIEXSecFile.Close()
