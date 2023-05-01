import os
import sys
import ROOT
from ROOT import *
import PlotUtils
from functools import partial
from config.AnalysisConfig import AnalysisConfig
from tools import Utilities
from tools.PlotLibrary import HistHolder
from config.PlotConfig import ELECTRON_ENERGY_BINNING, INLINE_UPSTREAM_ENERGY_BINS
from PlotUtils.HistWrapper import HistWrapper
import tools.TruthTools as TruthTools

def GetDfrWeight(event):
    if event.mc_intType == 10:
        dfrFile = ROOT.TFile.Open('/minerva/data/users/ajball/nu_e/kin_dist_mcdfrAlex_rebinned_nx_collab1_fspline.root')
        weightHist = dfrFile.Get('weight')
        elEnergy = event.kin_cal.reco_E_lep
        pionE = TruthTools.PiZeroE(event)       
#        if(pionE < 2):
#            eBin = 0
#        elif(pionE < 2.5):
#            eBin = 1
#        elif(pionE < 3):
#            eBin = 2
#        elif(pionE < 3.5):
#            eBin = 3
#        elif(pionE < 4):
#            eBin = 4
#        elif(pionE < 5):
#            eBin = 5
#        elif(pionE < 6):
#            eBin = 6
#        elif(pionE < 7):
#            eBin = 7
#        elif(pionE < 10):
#            eBin = 8
#        elif(pionE < 14):
#            eBin = 9
#        else:
#            eBin = 10'''
        if(pionE < 2):
            eBin = 0
        elif(pionE < 2.5):
            eBin = 1
        elif(pionE < 3.25):
            eBin = 2
        elif(pionE < 4):
            eBin = 3
        elif(pionE < 5):
            eBin = 4
        elif(pionE < 6):
            eBin = 5
        elif(pionE < 7):
            eBin = 6
        elif(pionE < 10):
            eBin = 7
        elif(pionE < 14):
            eBin = 8
        else:
            eBin = 9 
        #eBin = weightHist.FindBin(elEnergy)
        #print(elEnergy, eBin, weightHist.GetBinContent(eBin - 1))
        #print elEnergy, eBin
        #print weightHist.GetBinContent(eBin)
	
        return weightHist.GetBinContent(eBin)
    else:
        return 1.0    


