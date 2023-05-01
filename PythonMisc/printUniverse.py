import os
import sys
import ROOT
import PlotUtils
import math
import copy
from array import array
from collections import OrderedDict
import numpy as np

from tools.PlotLibrary import HistHolder
from config.AnalysisConfig import AnalysisConfig
from tools import Utilities,PlotTools
from config.SignalDef import SIGNAL_DEFINATION

def getPOTFromFile(filename):
    metatree = ROOT.TChain("Meta")
    if metatree.Add(filename,-1):
        return ROOT.PlotUtils.POTCounter().getPOTfromTChain(metatree)
    else:
        return None

def GetPOTScale(data_path,mc_path): # If you're a future REU student reading this, I have no idea what this does but it works now
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


weighting_playlist = "Alex_systematics_nx"
mc_path = "/minerva/data/users/ajball/nu_e/kin_dist_mc"+str(weighting_playlist)+"_collab1_fspline.root"
mc_file = ROOT.TFile.Open(mc_path)
data_path = "/minerva/data/users/ajball/nu_e/kin_dist_data"+str(weighting_playlist)+"_collab1_fspline.root"
data_file = ROOT.TFile.Open(data_path)
pot_scale = GetPOTScale(data_path,mc_path)

for i in range(5):
    outfile = ROOT.TFile.Open("/minerva/data/users/ajball/nu_e/kin_dist_mcFlux"+str(i)+"_collab1_fspline.root","RECREATE")
    sidebands = ["Coherent","Electron_Neutrino","Pi0","Signal"]
    for sideband in sidebands:
        mcEel_hists = HistHolder("Lepton Energy",mc_file,sideband,is_mc = True,pot = pot_scale) # Gets the two electron energy histograms
        dataEel_hist = HistHolder("Lepton Energy",data_file,sideband,is_mc = False) # I'm using HistHolders to access the backgrounds
        mcEel_hists.POTScale(False)
        dataEel_hist.POTScale(False)
        outhistEel = mcEel_hists.GetVertErrorBand("Flux").GetHist(i)
        outfile.Write(outhistEel)
        outhistIUE = HistHolder("Inline Upstream Energy",mc_file,sideband,is_mc = True,pot = pot_scale)
        outfile.Write(outhistIUE)
    outfile.Close()


mc_file.Close()
data_file.Close()
