import os
import sys
import ROOT
from ROOT import *
import PlotUtils
from functools import partial
from config.AnalysisConfig import AnalysisConfig
from tools import Utilities
from tools.PlotLibrary import HistHolder
from PlotUtils.HistWrapper import HistWrapper
import tools.TruthTools as TruthTools

IsCC = lambda event: event.mc_current == 1
IsNC = lambda event: event.mc_current == 2
IsCoherent = lambda event: event.mc_intType == 4
IsNuE = lambda event: abs(event.mc_incoming) == 12
IsPi0InFinalState = lambda event: 111 in event.mc_FSPartPDG

IsHeavyBaryon = lambda event: 3112 in event.mc_FSPartPDG or 3122 in event.mc_FSPartPDG or 3212 in event.mc_FSPartPDG or 3222 in event.mc_FSPartPDG or 4112 in event.mc_FSPartPDG or 4122 in event.mc_FSPartPDG or 4212 in event.mc_FSPartPDG or 4222 in event.mc_FSPartPDG
IsMeson = lambda event: 211 in event.mc_FSPartPDG or -211 in event.mc_FSPartPDG or 321 in event.mc_FSPartPDG or -321 in event.mc_FSPartPDG or 323 in event.mc_FSPartPDG or -323 in event.mc_FSPartPDG  or 111 in event.mc_FSPartPDG or 130 in event.mc_FSPartPDG or 310 in event.mc_FSPartPDG or 311 in event.mc_FSPartPDG
IsPhoton = lambda event: 22 in event.mc_FSPartPDG and event.mc_FSPartPDG[0] != 22

IsCCQEnu_e = lambda event: IsCC(event) and IsNuE(event) and not IsPi0InFinalState(event) and not IsMeson(event) and not IsHeavyBaryon(event) and not IsPhoton(event)
IsNCPi0 = lambda event: IsNC(event) and IsPi0InFinalState(event)

def GetAlexWeights(event):
    if IsCoherent(event): # Coherent event
        Eel = event.kin_cal.reco_E_lep # Gets electron cone energy
        if Eel < 2:
            weight = 1.7
        elif Eel < 2.5:
            weight = 1.77
        elif Eel < 3.25:
            weight = 1.81
        elif Eel < 4:
            weight = 1.97
        elif Eel < 5:
            weight = 2.45
        elif Eel < 6:
            weight = 1.55
        elif Eel < 7:
            weight = 2.13
        elif Eel <= 14:
            weight = 3
        else: weight = 1

    elif IsCCQEnu_e(event): # CCQE nu_e event
        Eel = event.kin_cal.reco_E_lep # Gets electron cone energy
        if Eel < 2:
            weight = 1.25
        elif Eel < 2.5:
            weight = 1.4
        elif Eel < 3.25:
            weight = 1.5
        elif Eel < 4:
            weight = 1.4
        elif Eel < 5:
            weight = 1.1
        elif Eel < 6:
            weight = 1.5
        elif Eel < 7:
            weight = 1.35
        elif Eel < 10:
            weight = 1.25
        elif Eel < 14:
            weight = 1.4
        elif Eel <= 24:
            weight = 1.75
        else: weight = 1

        if 10 <= event.UpstreamInlineEnergy < 15:
            weight *= 1.15
        elif 15 <= event.UpstreamInlineEnergy < 25:
            weight *= 1.2
    elif IsNCPi0(event):
        Eel = event.kin_cal.reco_E_lep # Gets electron cone energy
        if Eel < 2:
            weight = 0.9
        elif Eel < 2.5:
            weight = 0.8
        elif 4 < Eel < 7:
            weight = 1.25
        elif 7 < Eel < 14:
            weight = 4.5 
        else: weight = 1

    else: weight = 1
    return weight
