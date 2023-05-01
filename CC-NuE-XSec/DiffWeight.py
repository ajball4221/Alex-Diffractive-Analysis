import os
import sys
import ROOT
import PlotUtils
from functools import partial
from config.AnalysisConfig import AnalysisConfig
from tools import Utilities
from tools.PlotLibrary import HistHolder
from config.PlotConfig import ELECTRON_ENERGY_BINNING, INLINE_UPSTREAM_ENERGY_BINS
from PlotUtils.HistWrapper import HistWrapper
import math
import array

channel = "Inline Upstream vs Electron Energy"

REGION = "Signal"

PLAYLIST_TO_WRITE = ["simAlex_rebinned_nx"]
PLAYLIST_TO_WRITEDFR = ["dfrAlex_rebinned_nx"]

fitPath = '/minerva/data/users/ajball/nu_e/'

def CalcDfrExcess(mc_hist,data_hist,channel):
    data_hist = data_hist[channel]
    excess_hist = data_hist.Clone()
    excess_hist.AddMissingErrorBandsAndFillWithCV(mc_hist[channel])
    excess_hist.Add(mc_hist[channel],-1.0)
    return excess_hist

def GetMCHist(channel,data_file,mc_file):
    print("Channel: ", channel)
    print("data_file: ", data_file)
    print("mc_file: ", mc_file)
    hist = HistHolder(channel,mc_file,REGION,True)
    mc_hist = hist.GetHist().Clone()
    return mc_hist

def GetDataHist(channel,data_file,mc_file):
    hist = HistHolder(channel,data_file,REGION,True)
    data_hist = hist.GetHist().Clone()
    return data_hist

def GetShapeWeight(num,den):
    local_num = num.Clone()
    local_den = den.Clone()
    norm_num = local_num.Integral(0,-1)
    norm_den = local_den.Integral(0,-1)

    local_num.Divide(local_num,local_den)
    return local_num, norm_num,norm_den

def GetDfrWeight(num,den):
    local_num = num.Clone()
    local_den = den.Clone()
    norm_num = local_num.Integral(0,-1)
    norm_den = local_den.Integral(0,-1)

    local_num.Divide(local_num,local_den)
    return local_num, norm_num, norm_den

def GetDfrWeightHist():
    fdfr = ROOT.TFile.Open("{}/{}".format(fitPath,weight_dfr))
    weighthist = fdfr.Get("weight_Eel_HighInline_Excess")
    fdfr.Close()
    return weighthist


def DiffWeight(excess, dfr):
    weights = []
    excess_tot = []
    dfr_tot = []

    for i in range(0, 12):
        d_count = 0
        e_count = 0
        for j in range(0, 16):
            ex = excess.GetBinContent(j, i)
            diff = dfr.GetBinContent(j, i)
            d_count += diff
            e_count += ex
            #print "ex: ", ex, "diff: ", diff, "el bin: ", i, "inline bin: ", j
        if d_count != 0:
            weight = abs(e_count/d_count)
            print("Excess: ", ex, " Dfr: ", diff, " Weight: ", weight)
        else:
            weight = 0

        weights.append(weight)
        excess_tot.append(e_count)
        dfr_tot.append(d_count)
    print("Diffractive: ", dfr_tot)
    return weights, excess_tot, dfr_tot 


def CalcUncert(data, dfr):
    u = []
    for i in range(0, 19):
        da = 0
        df = 0
        for j in range(1, 16):
            data_count = data.GetBinContent(j, i)
            dfr_count = dfr.GetBinContent(j, i)
            da += data_count
            df += dfr_count
        if(df != 0):
            nextVal = math.sqrt(da)/df
            u.append(nextVal)
    return u


def eventCount(data, mc):
    d = []
    m = []
    print("Event Count!!!")
    for i in range(0, 19):
        d_count = 0
        m_count = 0
        for j in range(1, 16):
            dct = data.GetBinContent(j, i)
            mct = mc.GetBinContent(j, i)
            d_count += dct
            m_count += mct
        d.append(d_count)
        m.append(m_count)
    print("Data count: ", d[1:12])
    print("Monte carlo count: ", m[1:12])
    print("\n")
    return d, m


if __name__ == "__main__":

    #data_playlist = 'Alex_nx' #not used don't worry about it
    dfr_playlist = 'dfrAlex_rebinned_nx'
    
    data_file = ROOT.TFile.Open('/minerva/data/users/ajball/nu_e/kin_dist_dataAlex_rebinned_nx_collab1_fspline.root')
    mc_file = ROOT.TFile.Open('/minerva/data/users/ajball/nu_e/kin_dist_mcsimAlex_rebinned_nx_collab1_fspline.root')
    dfr_file = ROOT.TFile.Open('/minerva/data/users/ajball/nu_e/kin_dist_mcdfrAlex_rebinned_nx_collab1_fspline.root')

    data_pot = Utilities.getPOTFromFile('/minerva/data/users/ajball/nu_e/kin_dist_dataAlex_rebinned_nx_collab1_fspline.root')
    mc_pot = Utilities.getPOTFromFile('/minerva/data/users/ajball/nu_e/kin_dist_mcsimAlex_rebinned_nx_collab1_fspline.root')
    dfr_pot = Utilities.getPOTFromFile('/minerva/data/users/ajball/nu_e/kin_dist_mcdfrAlex_rebinned_nx_collab1_fspline.root')
 
    mc_hists = dict(zip([channel],map(partial(GetMCHist,data_file=data_file, mc_file=mc_file), [channel]))) 
    data_hists = dict(zip([channel],map(partial(GetDataHist,data_file=data_file, mc_file=mc_file), [channel])))
    dfr_hists = dict(zip([channel],map(lambda channel: HistHolder(channel,dfr_file,REGION,False).GetHist(),[channel])))
    
    mc_hists[channel].Scale(data_pot/mc_pot)
    dfr_hists[channel].Scale(mc_pot/dfr_pot)
    dfr_hists[channel].Scale(data_pot/mc_pot)
    #print" "
    #print(channel)
    #print(mc_hists[channel])
    #print(dfr_hists[channel])
    #print(data_hists[channel])
    #print " "

    for (coh_playlist, dfr_playlist) in zip(PLAYLIST_TO_WRITE, PLAYLIST_TO_WRITEDFR):
        TruthBin = []
        TruthBin.append(False)

        while all(TruthBin) == False:
            TruthBin = []
            #coh_weight = GetWeightHist()
            #dfr_weight = GetDfrWeightHist()
            dfr_weight = False
            #for channel in CHANNELS:
            myarray = ELECTRON_ENERGY_BINNING    
            
            dfr_hists[channel].AddMissingErrorBandsAndFillWithCV(mc_hists[channel])

            if dfr_weight:
                mc_hists[channel].Add(dfr_hists[channel],-1.0)
            excess_hist = CalcDfrExcess(mc_hists,data_hists,channel)
            #print(mc_hists[channel])
            data_mc_counts = eventCount(data_hists[channel], mc_hists[channel])
            diff_weights = DiffWeight(excess_hist, dfr_hists[channel])[0]
            uncertainties = CalcUncert(data_hists[channel], dfr_hists[channel])
            print('Diff_weights: ', diff_weights)
            print('Uncertainties: ', uncertainties) 
            h = ROOT.TH1F("h", "stuff", len(myarray), array.array('d', myarray))
            
            diff_weights = diff_weights[1:]
            for i in range(0, 11):
                print(i, "Weight: ", diff_weights[i])
                h.SetBinContent(i, diff_weights[i])
            
            dfr_file = ROOT.TFile.Open('/minerva/data/users/ajball/nu_e/kin_dist_mcdfrAlex_rebinned_nx_collab1_fspline.root', "recreate")
            h.Write("weight")
            dfr_file.Close()
            print("data: ", data_mc_counts[0])
            print("mc: ", data_mc_counts[1])
            print("weights: ", diff_weights[0])
            print("excess: ", diff_weights[1])
            print("diff: ", diff_weights[2])
            #print(weights) 
