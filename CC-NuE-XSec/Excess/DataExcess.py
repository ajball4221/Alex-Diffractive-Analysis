
"""
Characterize the dE/dX excess in data

"""

import os
del os.environ['DISPLAY']# don't turn on x11. this slow down processing.
import sys
import ROOT
import PlotUtils
#insert path for modules of this package.
from config import PlotConfig
from config.AnalysisConfig import AnalysisConfig
from tools import Utilities,PlotTools,PlotLibrary
from config import SystematicsConfig


MNVPLOTTER = PlotUtils.MnvPlotter()
#config MNVPLOTTER:
MNVPLOTTER.draw_normalized_to_bin_width=False
MNVPLOTTER.legend_text_size = 0.02
MNVPLOTTER.extra_top_margin = -.035# go slightly closer to top of pad
MNVPLOTTER.mc_bkgd_color = 46
MNVPLOTTER.mc_bkgd_line_color = 46

MNVPLOTTER.data_bkgd_color = 12 #gray
MNVPLOTTER.data_bkgd_style = 24 #circle  

#legend entries are closer
MNVPLOTTER.height_nspaces_per_hist = 1.2
MNVPLOTTER.width_xspace_per_letter = .4
MNVPLOTTER.legend_text_size        = .03

CANVAS = ROOT.TCanvas("c2","c2",1600,1000)

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)

EXCESS_PLOTS = {
    "Electron Theta":
    {   
        "plot_type" : ["difference"]
    },
    "Electron Energy":
    {   
        "plot_type" : ["difference"],
    },
    "EeTheta":
    {
        "plot_type" : ["difference"]
    },
    "Median Plane Shower Width":
    { 
        "plot_type" : ["difference"]
    },
    "Transverse Shower Asymmetry":
    {
        "plot_type" : ["difference"]
    },     
}


ADDITIONAL_DIR_STUBS = {
    "PC": {
        "PCPhoton": "me6a_PCPhoton_FlatFinal_nx",
        "PCPi0": "me6a_PCPi0_FlatFinal_nx",
        "PCElectron": "me6a_PCElectron_FlatFinal_nx",
    },
}


SELECTED_SIDEBANDS = AnalysisConfig.sidebands

class HistHolder:
    def __init__(self, name, f,is_mc = False):
        self.plot_name = PlotLibrary.PLOT_SETTINGS[name]["name"]
        self.f = f
        self.is_mc = is_mc
        self.hists= {}
        if f is not None:
            self.ReadHistograms()

    def ReadHistograms(self):
        for side_band in ["Signal"]+SELECTED_SIDEBANDS:
            namestring = "{}{}".format(self.plot_name,"" if side_band =="Signal" else "_"+side_band)
            self.hists[side_band] = Utilities.GetHistogram(self.f, namestring)
            if self.is_mc:
                for cate in PlotConfig.INT_COLORS:
                    cate_namestring = namestring+"_"+cate
                    self.hists[side_band+cate]=Utilities.GetHistogram(self.f, cate_namestring)

    def Scale(self,pot_scale, bin_width_normalize = True):
        if bin_width_normalize:
            for k, v in self.hists.items():
                if v is not None:
                    v.Scale(pot_scale,"width")
        else:
            for k, v in self.hists.items():
                if v is not None:
                    v.Scale(pot_scale)

    def GetListHistInput(self,sideband):
        if self.f is None or self.hists[sideband] is None: 
            return None
        if isinstance(self.hists[sideband],ROOT.TH2D):
            return PlotTools.Make2DSlice(self.hists[sideband],False,False)
        else:  
            return [self.hists[sideband]]

    def GetCateTOArray(self,sideband):
        if self.f is None or self.hists[sideband] is None or not self.is_mc:
            return None,None
        _mc_ints = []
        _colors =[]
    
        for cate in list(PlotConfig.INT_COLORS.keys())[::-1]:
            temp = PlotTools.Make2DSlice(self.hists[sideband+cate],False,False) if isinstance(self.hists[sideband+cate],ROOT.TH2D) else [self.hists[sideband+cate]]
            if len(_mc_ints) ==0:
                _mc_ints = [ROOT.TObjArray() for i in range(len(temp))]
            for i in range(len(temp)):
                hist = temp[i]
                hist.SetTitle(PlotConfig.INT_NAMES[cate])
                _mc_ints[i].Add(hist)
            _colors.append(PlotConfig.INT_COLORS[cate])
        return _mc_ints,_colors


def CalcExcess(data_hist, mc_ints , mc_hist, output,plot):
    
    have_mc = mc_ints is not None
    have_data = data_hist is not None
    
    if not have_mc:
        raise KeyError("Cannot make difference plot w/o MC")
    
    slices = data_hist if data_hist is not None else mc_ints
    N_plots = len(slices)

    for k in range(N_plots):
        data_h = data_hist[k]
        mc_h = mc_hist[k]
        data_h.Add(mc_h,-1) 
        
        return data_h


def FitExcessKinematics(excess_plot,plot):
    c = ROOT.TCanvas()
    f_fits = ROOT.TFile(AnalysisConfig.SelectionHistoPath(AnalysisConfig.playlist+"_ExcessSplineFits", True, True), "update")

    h_plot = excess_plot.GetCVHistoWithStatError()  
    spline = ROOT.TSpline3(h_plot)
    
    spline.Draw("P")  
    MNVPLOTTER.MultiPrint(c,  AnalysisConfig.PlotPath(plot,"fit")) 

    excess_plot.Write("%s" % (plot))
    espline_h = spline.GetHistogram()
 
    return espline_h.Clone()
    
   
def FitPCToExcess(pc_hists,output,excess,plot,playlist,excess_spline):
    c = ROOT.TCanvas("c2","c2",1600,1000)
    t_fits = ROOT.TFile(AnalysisConfig.SelectionHistoPath(AnalysisConfig.playlist+"_PCPolynomialSF", True, True), "update")
   
    have_pc = pc_hists is not None
    N_plots = len(pc_hists)
 
    for k in range(N_plots):
        pc = pc_hists[k] 

        pc_h = pc.GetCVHistoWithStatError()
        #pc_h.GetXaxis().SetLimits(1.5,10) 
        pcspline = ROOT.TSpline3(pc_h)

        pcspline.Draw("P")
        MNVPLOTTER.MultiPrint(c,  AnalysisConfig.PlotPath(plot,"PCFit"))

        pcspline_h = pcspline.GetHistogram()   
        pcspline_sf = pcspline_h.Clone()

        excessspline_sf = excess_spline.Clone() 
        #pcspline_sf.Divide(excessspline_sf) #Only here for the option to produce TSpline poly sf 
        
        pcspline_sf.GetXaxis().SetTitle("Excess / %s sf" % (playlist)) 
        pcspline_sf.Draw()
        
        pc_h.Write("%s_%s" % (plot,playlist))

    return excess, pc, pcspline_sf

def MakeShapeComparisons(plotname,excess,pc_h,tag,c,pcspline_sf):
   
    leg = ROOT.TLegend(0.55, 0.65, 0.85, 0.9)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)


    if pc_h.Integral() == 0:
        return
    
    pc_h.Scale(1./pc_h.Integral())
    excess.Scale(1./excess.Integral())
    pc_h.GetYaxis().SetRangeUser(0,0.3)
    pc_h.GetYaxis().SetTitle("Fraction of Events") 


    for i in range(0,pc_h.GetNbinsX()): #going through root hoops. Spline is th1f, mnvh1dCV is th1d
        bins = pc_h.GetBinContent(i)
        center = pc_h.GetBinCenter(i)
        binPC = pcspline_sf.GetXaxis().FindBin(center) 
        sf = pcspline_sf.GetBinContent(binPC) 
        
    
    if tag == "PCPi0":  
        pc_h.SetFillStyle(0)
        pc_h.SetLineColor(ROOT.kGreen)
        pc_h.SetLineWidth(3)
        pc_h.SetTitle("PCPi0")         
        pc_h.Draw("same hist")
        
    if tag == "PCPhoton":
        pc_h.SetFillStyle(0)
        pc_h.SetLineColor(ROOT.kBlue)
        pc_h.SetLineWidth(3)
        pc_h.SetTitle("PCPhoton")
        pc_h.Draw("HIST same")

    if tag == "PCElectron":
        pc_h.Scale(1./pc_h.Integral())
        pc_h.SetFillStyle(0)
        pc_h.SetLineColor(ROOT.kRed)
        pc_h.SetLineWidth(3)
        pc_h.SetTitle("PCElectron")
        pc_h.Draw("HIST same")
      
        excess.Draw("same")
   
    ROOT.gPad.BuildLegend() 
    #MNVPLOTTER.WriteNorm("Unit Normalized", "TL", 0.035)
    MNVPLOTTER.MultiPrint(c,  AnalysisConfig.PlotPath(plot,"overlay"))
   

if __name__ == "__main__":

    #input knobs
    playlist=AnalysisConfig.playlist
    
    if len(AnalysisConfig.data_types) == 2:
        try:
            data_pot=Utilities.getPOT(playlist,"data",AnalysisConfig.ntuple_tag)[0]
            mc_pot=Utilities.getPOT(playlist,"mc",AnalysisConfig.ntuple_tag)[0]
            pot_scale = data_pot/mc_pot

        except TypeError:
            # the playlist is not in dict, which means it is a merge playlist
            # merging multiple playlist already scaled pot
            pot_scale = 1.0

    else:
        # if only data or mc, there is no pot scale
        pot_scale =1.0
    
    data_file = ROOT.TFile.Open(AnalysisConfig.SelectionHistoPath(playlist,True,False)) if "data" in AnalysisConfig.data_types else None
    mc_file = ROOT.TFile.Open(AnalysisConfig.SelectionHistoPath(playlist,False,False)) if "mc" in AnalysisConfig.data_types else None

 
    for plot,config in EXCESS_PLOTS.items():
        data_hists = HistHolder(plot,data_file)
        mc_hists = HistHolder(plot,mc_file,True)
        
        #scale by pot and bin width
        data_hists.Scale(1.0)
        mc_hists.Scale(pot_scale)
        
        for sideband in ["Signal"]+SELECTED_SIDEBANDS:
            for i in config["plot_type"]:
                (mc_ints, colors) = mc_hists.GetCateTOArray(sideband)
            
        excess_plot = CalcExcess(data_hists.GetListHistInput(sideband),mc_ints,mc_hists.GetListHistInput(sideband),AnalysisConfig.PlotPath(data_hists.plot_name,sideband,i),plot) #this placement gets only sideband shape plots
        
        spline = FitExcessKinematics(excess_plot,plot) 
        
        c = ROOT.TCanvas()   
        for groupname, group in ADDITIONAL_DIR_STUBS.items():
            for tag, playlist in group.items():
                
                pc_file = ROOT.TFile.Open(AnalysisConfig.SelectionHistoPath(playlist,False,False)) #if "data" in AnalysisConfig.data_types else None
                pc_hists = HistHolder(plot,pc_file)
                
                excess_h, pc_h, pcspline_sf = FitPCToExcess(pc_hists.GetListHistInput(sideband),AnalysisConfig.PlotPath(data_hists.plot_name,sideband,i),excess_plot,plot,playlist,excess_plot) #Added this
               
                MakeShapeComparisons(plot, excess_h, pc_h,tag,c,pcspline_sf) 


