import ROOT
import PlotUtils

ROOT.TH1.AddDirectory(False)
ROOT.gROOT.SetBatch(True)
MNVPLOTTER = PlotUtils.MnvPlotter()
CANVAS = ROOT.TCanvas("canvas","canvas",1600,1000)

def MakeErrPlot(hist,mnvplotter=MNVPLOTTER,canvas=CANVAS):
    mnvplotter.axis_maximum = 0.25 #round((hist.GetTotalError().GetMaximum()/hist.GetBinContent(hist.GetTotalError().GetMaximumBin()))+0.05,1)
    mnvplotter.DrawErrorSummary(hist)

def updatePlotterErrorGroup(mnvplotter=MNVPLOTTER):
    group = {'Detector model': ['eltheta', 'beam_angle', 'elE_ECAL', 'elE_HCAL', 'birks', 'response_p', 'response_meson', 'response_em', 'response_other', 'response_xtalk', 'Leakage_Uncertainty', 'Target_Mass_CH', 'GEANT_Neutron', 'GEANT_Pion', 'GEANT_Proton'], 'Interaction model': ['GENIE_AGKYxF1pi', 'GENIE_AhtBY', 'GENIE_BhtBY', 'GENIE_CCQEPauliSupViaKF', 'GENIE_CV1uBY', 'GENIE_CV2uBY', 'GENIE_EtaNCEL', 'GENIE_MaNCEL', 'GENIE_MaRES', 'GENIE_MvRES', 'GENIE_NormCCRES', 'GENIE_NormDISCC', 'GENIE_NormNCRES', 'GENIE_RDecBR1gamma', 'GENIE_Rvn2pi', 'GENIE_Rvp2pi', 'GENIE_Theta_Delta2Npi', 'GENIE_VecFFCCQEshape', 'GENIE_MaCCQE', 'GENIE_Rvn1pi', 'GENIE_Rvp1pi', 'GENIE_FrAbs_N', 'GENIE_FrAbs_pi', 'GENIE_FrCEx_N', 'GENIE_FrCEx_pi', 'GENIE_FrElas_N', 'GENIE_FrElas_pi', 'GENIE_FrInel_N', 'GENIE_FrPiProd_N', 'GENIE_FrPiProd_pi', 'GENIE_MFP_N', 'GENIE_MFP_pi'], 'MnvTunes': ['RPA_HighQ2', 'RPA_LowQ2', 'Low_Recoil_2p2h_Tune', 'LowQ2Pi'], 'Alternative Tunning methods': ['bkg_tune']}
    mnvplotter.error_summary_group_map.clear();
    for k,v in group.items():
        vec = ROOT.vector("std::string")()
        for vs in v :
            vec.push_back(vs)
        mnvplotter.error_summary_group_map[k]= vec

#configuring MNVPLOTTER
MNVPLOTTER.draw_normalized_to_bin_width=False
MNVPLOTTER.legend_text_size = 0.04
MNVPLOTTER.extra_top_margin = -.035 # go slightly closer to top of pad
MNVPLOTTER.mc_bkgd_color = 46
MNVPLOTTER.mc_bkgd_line_color = 46

MNVPLOTTER.data_bkgd_color = 12 #gray
MNVPLOTTER.data_bkgd_style = 24 #circle

#legend entries are closer
MNVPLOTTER.height_nspaces_per_hist = 1.2
MNVPLOTTER.width_xspace_per_letter = .4

MNVPLOTTER.good_colors.pop_back()
updatePlotterErrorGroup()


dataFile = ROOT.TFile.Open("datafolder/nu_e/Alex_final_nx_calculatedEfficiency.root","READ")
#hists_to_print = ["t_fitted","t_recovered","pionE_recovered","IUE_fit_prediction","IUE_matrixpredicted","Eel_matrixpredicted"]
#hists_to_print = ["Eel","InlineUpstream_Energy"]
hists_to_print = ["t_efficiency","pion_efficiency"]
for i in range(2):
    hist = eval("dataFile."+hists_to_print[i])
    CANVAS = ROOT.TCanvas("canvas"+str(i),"canvas"+str(i),1600,1000)
    MakeErrPlot(hist)
    CANVAS.Print("datafolder/errplots/"+hists_to_print[i]+".png","png")

