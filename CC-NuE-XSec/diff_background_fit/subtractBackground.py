import os
import sys
import ROOT
import PlotUtils
import math
import copy
from array import array
from collections import OrderedDict

from tools.PlotLibrary import HistHolder
from config.AnalysisConfig import AnalysisConfig
from tools import Utilities,PlotTools
from config.SignalDef import SIGNAL_DEFINATION

def BackgroundSubtraction(data_hists, mc_hists):
    data_hists.POTScale(False)
    mc_hists.POTScale(False)
    out_data = data_hists.GetHist().Clone()
    out_mc = mc_hists.hists["Total"].Clone()
    out_mc.GetVertErrorBand("Flux").GetHist(0).Print("all")
    out_data.AddMissingErrorBandsAndFillWithCV(out_mc)
    out_data.Print("all")
    out_data.GetVertErrorBand("Flux").GetHist(0).Print("all")
    first_subtraction = True
    for group in mc_hists.hists:
        if group == "Total":
                continue
        elif group not in SIGNAL_DEFINATION:
            out_data = SubtractPoissonHistograms(out_data,mc_hists.hists[group])
            out_mc = SubtractPoissonHistograms(out_mc,mc_hists.hists[group],first_subtraction)
            first_subtraction = False
    for i in range(out_data.GetSize()): 
        if out_data.GetBinContent(i) < 0: 
            out_data.SetBinContent(i,0.01)
    return out_data,out_mc

def SubtractPoissonHistograms(signal,background_band,first_subtraction = False): # Sets error of signal where signal - background_band is being calculat
    errors = []
    for i in range(signal.GetSize()): # Adds errors in quadrature
        if first_subtraction: errors.append(math.sqrt(math.sqrt(pot_scale)*signal.GetBinError(i)**2 + math.sqrt(pot_scale)*background_band.GetBinError(i)**2)) # Scales the mc total error by the POT scale just once
        else: errors.append(math.sqrt(signal.GetBinError(i)**2 + math.sqrt(pot_scale)*background_band.GetBinError(i)**2))
    signal.Add(background_band,-1) # Subtracts the band from signal
    for i in range(signal.GetSize()): # Sets the error of signal to be what was calculated before
        signal.SetBinError(i,errors[i])
    return signal

def GetBins(hist,print_results = False):
    bins = []
    for i in range(hist.GetSize()):
        bins.append(hist.GetBinContent(i))
    if print_results: print(hist,bins)
    return bins

def getPOTFromFile(filename):
    metatree = ROOT.TChain("Meta")
    if metatree.Add(filename,-1):
        return ROOT.PlotUtils.POTCounter().getPOTfromTChain(metatree)
    else:
        return None

def GetPOTScale(data_path,mc_path):
    pots = [None,None]
    for i,t in enumerate(["data","mc"]):
        path = [data_path,mc_path][i]
        try:
            pots[i]= getPOTFromFile(path) or getPOT(playlist,t,ntuple_tag)
        except KeyError:
            pots[i]=None
    pot_scale = pots[0]/pots[1] if pots.count(None) == 0 else 1.0
    print(pot_scale)
    return pot_scale

playlist = AnalysisConfig.playlist
mc_path = "/minerva/data/users/ajball/nu_e/kin_dist_mc"+str(playlist)+"_collab1_fspline.root"
mc_file = ROOT.TFile.Open(mc_path)
data_path = "/minerva/data/users/ajball/nu_e/kin_dist_data"+str(playlist)+"_collab1_fspline.root"
data_file = ROOT.TFile.Open(data_path)

pot_scale = GetPOTScale(data_path,mc_path)

mcEel_hists = HistHolder("Lepton Energy",mc_file,"Signal",is_mc = True,pot = pot_scale) # Gets the two electron energy histograms
dataEel_hist = HistHolder("Lepton Energy",data_file,"Signal",is_mc = False) # I'm using HistHolders so I can access the backgrounds

mcIUE_hists = HistHolder("Inline Upstream Energy",mc_file,"Signal",is_mc = True,pot = pot_scale) # Gets the two Inline Upstream Energy Histograms
dataIUE_hist = HistHolder("Inline Upstream Energy",data_file,"Signal",is_mc = False)

new_dataEel_hist,new_mcEel_hists = BackgroundSubtraction(dataEel_hist,mcEel_hists)
new_dataIUE_hist,new_mcIUE_hists = BackgroundSubtraction(dataIUE_hist,mcIUE_hists)
GetBins(new_dataEel_hist,True)
GetBins(new_dataIUE_hist,True)

mc_file.Close()
data_file.Close()

output_datafile = ROOT.TFile.Open("/minerva/data/users/ajball/nu_e/kin_dist_data"+str(playlist)+"_subtracted_collab1_fspline.root","RECREATE")
new_dataEel_hist.GetXaxis().SetTitle("Electron Cone Energy (GeV)")
new_dataEel_hist.GetXaxis().SetTitleFont(62)
new_dataEel_hist.GetXaxis().SetTitleSize(0.045)
new_dataEel_hist.GetYaxis().SetTitle("Background Subtracted Bin Counts")
new_dataEel_hist.GetYaxis().SetTitleFont(62)
new_dataEel_hist.GetYaxis().SetTitleSize(0.045)
new_dataEel_hist.SetTitle("Subtracted Electron Cone Energy Data")
new_dataEel_hist.SetTitleFont(62)
new_dataEel_hist.SetLabelFont(62,"xyz")
new_dataEel_hist.SetMarkerStyle(21)
new_dataEel_hist.SetMarkerSize(2)
new_dataEel_hist.SetLineWidth(4)
new_dataEel_hist.SetStats(0)
new_dataEel_hist.Write()
Eel_canvas = ROOT.TCanvas("c1","c1",1600,1200)
new_dataEel_hist.Draw()
Eel_canvas.Print("datafolder/nu_e/plot/Eel_subtracted.png","png")


new_dataIUE_hist.GetXaxis().SetTitle("Inline Upstream Energy (MeV)")
new_dataIUE_hist.GetXaxis().SetTitleFont(62)
new_dataIUE_hist.GetXaxis().SetTitleSize(0.045)
new_dataIUE_hist.GetYaxis().SetTitle("Background Subtracted Bin Counts")
new_dataIUE_hist.GetYaxis().SetTitleFont(62)
new_dataIUE_hist.GetYaxis().SetTitleSize(0.045)
new_dataIUE_hist.SetTitle("Subtracted Inline Upstream Energy Data")
new_dataEel_hist.SetTitleFont(62)
new_dataIUE_hist.SetLabelFont(62,"xyz")
new_dataIUE_hist.SetMarkerStyle(21)
new_dataIUE_hist.SetMarkerSize(2)
new_dataIUE_hist.SetLineWidth(4)
new_dataIUE_hist.SetStats(0)
new_dataIUE_hist.Write()
output_datafile.Close()
IUE_canvas = ROOT.TCanvas("c2","c2",1600,1200)
new_dataIUE_hist.Draw()
IUE_canvas.Print("datafolder/nu_e/plot/IUE_subtracted.png","png")

output_mcfile = ROOT.TFile.Open("/minerva/data/users/ajball/nu_e/kin_dist_mc"+str(playlist)+"_subtracted_collab1_fspline.root","RECREATE")
new_mcEel_hists.Write()
new_mcIUE_hists.Write()
output_mcfile.Close()

