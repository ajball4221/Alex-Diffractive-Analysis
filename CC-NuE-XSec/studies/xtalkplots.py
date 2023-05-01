import ROOT
import math
from functools import partial
import PlotUtils
from config import PlotConfig
from config.AnalysisConfig import AnalysisConfig
from config.DrawingConfig import Default_Plot_Type,Default_Scale,DefaultSlicer
from tools import Utilities,PlotTools
from tools.PlotLibrary import HistHolder
ROOT.TH1.AddDirectory(False)


PLOTS_TO_MAKE = [
    {"name":
     "Visible Energy"},
     {"name":
     "Lepton Pt"},
     {"variables":["Visible Energy","Lepton Pt"],}
]

def xtalkCompPlotter(mnvplotter,data_hist,mc_hist):
    tmpMC = PlotUtils.MnvH1D(mc_hist.GetCVHistoWithStatError())
    tmpXtalk = mc_hist.GetVertErrorBand("response_xtalk")
    tmpMC.PushErrorBand("xtalk",tmpXtalk)
    mnvplotter.DrawDataMCWithErrorBand(data_hist.GetCVHistoWithStatError(),tmpMC.GetCVHistoWithError(False),1.0,"TR")

def xtalkRatioPlotter1(mnvplotter,data_hist,mc_hist):
    tmpMC = PlotUtils.MnvH1D(mc_hist.GetCVHistoWithStatError())
    tmpXtalk = mc_hist.GetVertErrorBand("response_xtalk")
    tmpMC.PushErrorBand("xtalk",tmpXtalk)
    tmpDe = PlotUtils.MnvH1D(mc_hist.GetCVHistoWithStatError())
    tmpDe.AddMissingErrorBandsAndFillWithCV(tmpMC)
    tmpData=PlotUtils.MnvH1D(data_hist.GetCVHistoWithStatError())
    tmpData.AddMissingErrorBandsAndFillWithCV(tmpMC)
    tmpData.Divide(tmpData,tmpDe)
    tmpMC.Divide(tmpMC,tmpDe)
    mnvplotter.DrawDataMCWithErrorBand(tmpData.GetCVHistoWithStatError(),tmpMC.GetCVHistoWithError(False),1.0,"TR")

def xtalkRatioPlotter(mnvplotter,data_hist,mc_hist):
    tmpMC = PlotUtils.MnvH1D(mc_hist.GetCVHistoWithStatError())
    tmpXtalk = mc_hist.GetVertErrorBand("response_xtalk")
    tmpMC.PushErrorBand("xtalk",tmpXtalk)
    tmpData=data_hist.GetCVHistoWithStatError()
    for i in range(tmpData.GetSize()):
        tmpData.SetBinError(i,0)
    mnvplotter.DrawDataMCRatio(tmpData,tmpMC.GetCVHistoWithError(False),1.0,True,0.8,1.2)

def MakePlot(data_hists,mc_hists,config):
    #scale
    if "scale" in config:
        config["scale"](data_hists)
        config["scale"](mc_hists)
    else:
        Default_Scale(data_hists)
        Default_Scale(mc_hists)
    CanvasConfig = config.setdefault("canvasconfig",lambda x:True)
    #PlotType = config.setdefault("plot_type",Default_Plot_Type)
    slicer = config.setdefault("slicer", DefaultSlicer(data_hists))
    draw_seperate_legend = config.setdefault("draw_seperate_legend",data_hists.dimension!=1)
    try:
        custom_tag = "xtalk_comp"
        plotfunction = xtalkCompPlotter
        hists = data_hists.GetHist(),mc_hists.GetHist()
        PlotTools.MakeGridPlot(slicer,plotfunction,hists,CanvasConfig,draw_seperate_legend)
        PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name,sideband,custom_tag))
        custom_tag = "xtalk_comp_ratio"
        plotfunction = xtalkRatioPlotter
        PlotTools.MakeGridPlot(slicer,plotfunction,hists,CanvasConfig,draw_seperate_legend)
        PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name,sideband,custom_tag))
        print("plot {} made.".format(data_hists.plot_name))
    except KeyError as e:
        print("plot {} not made.".format(data_hists.plot_name))
        raise e
        return False
    return True

if __name__ == "__main__":
    #input knobs
    playlist=AnalysisConfig.playlist
    type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag,True)
    standPOT = data_pot

    sideband_map = {}
    for config in PLOTS_TO_MAKE:
        sideband_group =  ["Signal"]
        if isinstance(sideband_group,list):
            for sideband in sideband_group:
                data_hists = HistHolder(config["name"] if "name" in config else config,data_file,sideband,False,data_pot,standPOT)
                mc_hists = HistHolder(config["name"] if "name" in config else config,mc_file,sideband,True,mc_pot,standPOT)
                MakePlot(data_hists,mc_hists,config)
        else:
            #assuing sideband_group is a tuple of name, and list of sidebands
            sideband = sideband_group[0]
            sidebands = sideband_group[1]
            data_hists = HistHolder(config["name"] if "name" in config else config,data_file,sidebands[0],False,data_pot,standPOT)
            mc_hists = HistHolder(config["name"] if "name" in config else config,mc_file,sidebands[0],True,mc_pot,standPOT)
            for _ in range(1,len(sidebands)):
                data_hists.Add(HistHolder(config["name"] if "name" in config else config,data_file,sidebands[_],False,data_pot,standPOT))
                mc_hists.Add(HistHolder(config["name"] if "name" in config else config,mc_file,sidebands[_],True,mc_pot,standPOT))
            MakePlot(data_hists,mc_hists,config)
