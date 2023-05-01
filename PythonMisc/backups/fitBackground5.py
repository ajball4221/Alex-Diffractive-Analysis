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
    print("POT Scale:",pot_scale)
    return pot_scale

def AddTwoHistograms(band1,band2): # Adds band2 to band 1 and adds the pot scaled errors in quadrature
    errors = []
    for i in range(band1.GetSize()): # Adds errors in quadrature
        errors.append(math.sqrt(math.sqrt(pot_scale)*band1.GetBinError(i)**2 + math.sqrt(pot_scale)*band2.GetBinError(i)**2))
    band1.Add(band2,1) # Adds band2 to band1
    for i in range(band1.GetSize()): # Sets the error of signal to be what was calculated before
        band1.SetBinError(i,errors[i])
    return band1

def SubtractPoissonHistograms(signal,background_band,first_subtraction = False): # Sets error of signal where signal - background_band is being calculat
    errors = []
    for i in range(signal.GetSize()): # Adds errors in quadrature
        if first_subtraction: errors.append(math.sqrt(math.sqrt(pot_scale)*signal.GetBinError(i)**2 + math.sqrt(pot_scale)*background_band.GetBinError(i)**2)) # Scales the mc total error by the POT scale just once
        else: errors.append(math.sqrt(signal.GetBinError(i)**2 + math.sqrt(pot_scale)*background_band.GetBinError(i)**2))
    signal.Add(background_band,-1) # Subtracts the band from signal
    for i in range(signal.GetSize()): # Sets the error of signal to be what was calculated before
        signal.SetBinError(i,errors[i])
    return signal

def PrimeExtractHists(data_hists, mc_hists):
    data = data_hists.GetHist().Clone()
    bands = ["NCCohPi0","NCDiff","CCQElike","notCCQElike","NCPi0","CCPi0"]
    coherent = mc_hists.hists["NCCohPi0"].Clone()
    ncdiff = mc_hists.hists["NCDiff"].Clone()
    nu_e = AddTwoHistograms(mc_hists.hists["CCQElike"].Clone(),mc_hists.hists["notCCQElike"].Clone())
    piZero = AddTwoHistograms(mc_hists.hists["NCPi0"].Clone(),mc_hists.hists["CCPi0"].Clone())
    total_mc = mc_hists.hists["Total"].Clone()
    for i,band in enumerate(bands): # Here, each of the subbands that are being calculated are removed from the total
        total_mc = SubtractPoissonHistograms(total_mc,mc_hists.hists[band],i) 
    augment = SubtractPoissonHistograms(data,total_mc) # This subtracts the LHS const from the data to give the RHS of the equation
    return [coherent,nu_e,piZero,ncdiff,augment]

def ErrExtractHists(data_hists, mc_hists, err_type, universe_num):
    #data_hists.POTScale(False)
    #mc_hists.POTScale(False)
    data = data_hists.GetHist().Clone()
    coherent = mc_hists.hists["NCCohPi0"].GetVertErrorBand(err_type).GetHist(universe_num).Clone()
    ncdiff = mc_hists.hists["NCDiff"].GetVertErrorBand(err_type).GetHist(universe_num).Clone()
    nu_e = AddTwoHistograms(mc_hists.hists["CCQElike"].GetVertErrorBand(err_type).GetHist(universe_num).Clone(),mc_hists.hists["notCCQElike"].GetVertErrorBand(err_type).GetHist(universe_num).Clone())
    piZero = AddTwoHistograms(mc_hists.hists["NCPi0"].GetVertErrorBand(err_type).GetHist(universe_num).Clone(),mc_hists.hists["CCPi0"].GetVertErrorBand(err_type).GetHist(universe_num).Clone())
    total_mc = mc_hists.hists["Total"].GetVertErrorBand(err_type).GetHist(universe_num).Clone()
    bands = ["NCCohPi0","NCDiff","CCQElike","notCCQElike","NCPi0","CCPi0"]
    for i,band in enumerate(bands): # Here, each of the subbands that are being calculated are removed from the total
        #print(band,universe_num)
        #mc_hists.hists[band].GetVertErrorBand(err_type).GetHist(universe_num).Clone().Print("all")
        total_mc = SubtractPoissonHistograms(total_mc,mc_hists.hists[band].GetVertErrorBand(err_type).GetHist(universe_num),not bool(i))
    other = RejectSmallBins(total_mc) # There were some bins that were like 1e-16 so I get rid of them here
    other = PlotUtils.MnvH1D(other) # I have no idea why the add command can't do this itself but whatever
    augment = SubtractPoissonHistograms(data,other) # This subtracts the LHS const from the data to give the RHS of the equation
    return [coherent,nu_e,piZero,ncdiff,augment]

def RejectSmallBins(hist):
    for i in range(hist.GetSize()):
        if hist.GetBinContent(i) < 1e-3:
            hist.SetBinContent(i,0)
    return hist

def HistListtoArray(histlist):
    histarray = np.zeros((len(histlist),histlist[0].GetSize()-1))
    for i,hist in enumerate(histlist):
        for j in range(hist.GetSize()-1):
            histarray[i,j] = hist.GetBinContent(j+1)
    if len(histlist) == 1:
        histarray = histarray.flatten()
    return histarray

def SetupMatrices(dataarray):
    numofmatrices = len(dataarray[0,0])
    matrices = np.zeros((numofmatrices,4,4))
    vectors = np.zeros((numofmatrices,4))
    for i in range(numofmatrices):
        for j in range(4):
            matrices[i,j] = dataarray[j,0:4,i]
            vectors[i,j] = dataarray[j,4,i]

    for i in range(numofmatrices):
        for j in range(len(matrices[0])):
            vectors[i,j] -= np.sum(matrices[i,j,0:4]) # This makes it so that the matrices are finding (weight-1) 
    return matrices,vectors

def DoLinearAlgebra(inversesArray,vectorArray):
    solutions = np.zeros(vectorArray.shape)
    for i in range(len(inversesArray)):
        #invMatrix = np.linalg.inv(matrixArray[i])
        solutions[i] = np.dot(inversesArray[i],vectorArray[i])
    final_weights = np.add(np.ones(solutions.shape),solutions) # Corrects for the fact that the matrices are solving for weight-1
    return final_weights

def MPInverse(matrices):
    svds = []
    for matrix in matrices: # Calculates all of the SVDs and puts them in a list
        decomposition = np.linalg.svd(matrix)
        svds.append(decomposition)
    discarded = np.array([[False,False,False,False],
                          [False,False,False,False],
                          [False,False,False,False],
                          [False,False,False,False],
                          [False,False,False,False],
                          [False,False,False,False],
                          [False,False,False,False],
                          [False,False,False, True],
                          [False,False,False, True],
                          [False, True, True, True],
                          [False,False, False, True]],dtype = bool) #Picks singular values to be discarded
    for i,svd in enumerate(svds):
        for j in range(len(svd[1])):
            if discarded[i,j]:
                svd[1][j] = 0
    inverses = []
    for svd in svds: # Constructs the MP Pseudoinverse using the SVD
        Ustar = np.transpose(svd[0])
        V = np.transpose(svd[2])
        sigma = svd[1]
        sigmaPlus = np.zeros((len(sigma),len(sigma)))
        for i in range(len(sigma)):
            if sigma[i] == 0:
                sigmaPlus[i,i] = 0
            else: sigmaPlus[i,i] = 1/sigma[i]
        mpinverse = np.matmul(np.matmul(V,sigmaPlus),Ustar)
        inverses.append(mpinverse)
    return inverses

def GetBins(hist,print_results = False): # If you want to print out the binds in a hist, use this
    bins = []
    for i in range(hist.GetSize()):
        bins.append(hist.GetBinContent(i))
    if print_results: print(hist,bins)
    return bins

np.set_printoptions(precision=6,suppress = True,threshold = np.inf)

playlist = AnalysisConfig.playlist
mc_path = "/minerva/data/users/ajball/nu_e/kin_dist_mc"+str(playlist)+"_collab1_fspline.root"
mc_file = ROOT.TFile.Open(mc_path)
data_path = "/minerva/data/users/ajball/nu_e/kin_dist_data"+str(playlist)+"_collab1_fspline.root"
data_file = ROOT.TFile.Open(data_path)
pot_scale = GetPOTScale(data_path,mc_path)

sidebands = ["Coherent","Electron_Neutrino","Pi0","Signal"]

# First going through the weights for the prime universe because I want to have this be its own thing
sideband_list = []
for sideband in sidebands:
    mcEel_hists = HistHolder("Lepton Energy",mc_file,sideband,is_mc = True,pot = pot_scale) # Gets the two electron energy histograms
    dataEel_hist = HistHolder("Lepton Energy",data_file,sideband,is_mc = False) # I'm using HistHolders to access the backgrounds
    mcEel_hists.POTScale(False)
    dataEel_hist.POTScale(False)
    listofhists = PrimeExtractHists(dataEel_hist, mcEel_hists)
    arrayofhists = HistListtoArray(listofhists)
    sideband_list.append(arrayofhists)
bigArray = np.array(sideband_list)
matrices,vectors = SetupMatrices(bigArray)
primeweights = DoLinearAlgebra(MPInverse(matrices),vectors)
del mcEel_hists, dataEel_hist

# Now for the Systematic Universes
err_weights_dict = {}
universe_list = [["Flux",200],["GENIE_AGKYxF1pi",2],["GENIE_AhtBY",2],["GENIE_BhtBY",2],["GENIE_CCQEPauliSupViaKF",2],["GENIE_CV1uBY",2],["GENIE_CV2uBY",2],["GENIE_EtaNCEL",2],["GENIE_FrAbs_N",2],["GENIE_FrAbs_pi",2],["GENIE_FrCEx_N",2],["GENIE_FrCEx_pi",2],["GENIE_FrElas_N",2],["GENIE_FrElas_pi",2],["GENIE_FrInel_N",2],["GENIE_FrPiProd_N",2],["GENIE_FrPiProd_pi",2],["GENIE_MFP_N",2],["GENIE_MFP_pi",2],["GENIE_MaNCEL",2],["GENIE_MaRES",2],["GENIE_MvRES",2],["GENIE_NormCCRES",2],["GENIE_NormDISCC",2],["GENIE_NormNCRES",2],["GENIE_RDecBR1gamma",2],["GENIE_Rvn2pi",2],["GENIE_Rvp2pi",2],["GENIE_Theta_Delta2Npi",2],["GENIE_VecFFCCQEshape",2],["GENIE_Rvn1pi",2],["GENIE_Rvp1pi",2],["GENIE_MaCCQE",2],["Low_Recoil_2p2h_Tune",3],["RPA_HighQ2",2],["RPA_LowQ2",2],["LowQ2Pi",2],["LowQ2Pi_None",1],["MK_model",1],["fsi_weight",3],["SuSA_Valencia_Weight",1],["GEANT_Proton",2],["GEANT_Neutron",2],["GEANT_Pion",2],["bkg_tune",2],["Target_Mass_CH",2],["Muon_Energy_MINERvA",2],["Muon_Energy_MINOS",2],["beam_angle",4],["response_p",2],["response_meson",2],["response_em",2],["response_other",2],["response_xtalk",2],["Leakage_Uncertainty",2]]
for error_type,error_universes in universe_list:
    print(error_type,error_universes)
    err_type_weights = []
    for universe_num in range(error_universes):
        sideband_list = []
        for sideband in sidebands:
            mcEel_hists = HistHolder("Lepton Energy",mc_file,sideband,is_mc = True,pot = pot_scale) 
            dataEel_hist = HistHolder("Lepton Energy",data_file,sideband,is_mc = False)
            listofhists = ErrExtractHists(dataEel_hist, mcEel_hists, error_type, universe_num)
            arrayofhists = HistListtoArray(listofhists)
            sideband_list.append(arrayofhists)
        bigArray = np.array(sideband_list)
        matrices,vectors = SetupMatrices(bigArray)
        err_type_weights.append(DoLinearAlgebra(MPInverse(matrices),vectors))
        del mcEel_hists, dataEel_hist
    err_weights_dict[error_type] = np.array(err_type_weights)

print(err_weights_dict)
print(err_weights_dict["Flux"][44])

'''np.savetxt("/minerva/app/users/ajball/ncdiff/CC-NuE-XSec/tools/weights.txt",weights,delimiter=",")
np.set_printoptions(precision=6,suppress = True)
print(np.array2string(weights,separator = ","))''' # Some stuff for exporting

mc_file.Close()
data_file.Close()
