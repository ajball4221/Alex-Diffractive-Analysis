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

def ExtractHists(data_hists, mc_hists):
    data_hists.POTScale(False)
    mc_hists.POTScale(False)
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
    return matrices,vectors

def DoLinearAlgebra(matrixArray,vectorArray):
    solutions = np.zeros(vectorArray.shape)
    for i in range(len(matrixArray)):
        invMatrix = np.linalg.inv(matrixArray[i])
        solutions[i] = np.dot(invMatrix,vectorArray[i])
    return solutions

def MPInverse(matrices):
    svds = []
    discard_ratio = 500 # Ratio at which an singular value will be removed
    for matrix in matrices: # Calculates all of the SVDs and puts them in a list
        decomposition = np.linalg.svd(matrix)
        svds.append(decomposition)
    for svd in svds:
        largestSV = svd[1][0] # If the ratio between the largest and current SV is too large, the SV gets set to 0
        for i in range(len(svd[1])):
            if largestSV/svd[1][i] >= discard_ratio: 
                svd[1][i] = 0
    inverses = []
    for svd in svds: # Constructs the MP Pseudoinverse using the SVD
        Ustar = np.transpose(svd[0])
        V = np.transpose(svd[2])
        sigma = svd[1]
        sigmaPlus = np.zeros((len(sigma),len(sigma)))
        for i in range(len(sigma)):
            sigmaPlus[i,i] = 1/sigma[i]
        mpinverse = np.matmul(np.matmul(V,sigmaPlus),Ustar)
        inverses.append(mpinverse)
    return inverses

def GetBins(hist,print_results = False): # If you want to print out the binds in a hist, use this
    bins = []
    for i in range(hist.GetSize()):
        bins.append(hist.GetBinContent(i))
    if print_results: print(hist,bins)
    return bins

np.set_printoptions(precision=2,suppress = True)

playlist = AnalysisConfig.playlist
mc_path = "/minerva/data/users/ajball/nu_e/kin_dist_mc"+str(playlist)+"_collab1_fspline.root"
mc_file = ROOT.TFile.Open(mc_path)
data_path = "/minerva/data/users/ajball/nu_e/kin_dist_data"+str(playlist)+"_collab1_fspline.root"
data_file = ROOT.TFile.Open(data_path)
pot_scale = GetPOTScale(data_path,mc_path)

sidebands = ["Coherent","Electron_Neutrino","Pi0","Signal"]
sideband_list = []
for sideband in sidebands:
    mcEel_hists = HistHolder("Lepton Energy",mc_file,sideband,is_mc = True,pot = pot_scale) # Gets the two electron energy histograms
    dataEel_hist = HistHolder("Lepton Energy",data_file,sideband,is_mc = False) # I'm using HistHolders to access the backgrounds
    listofhists = ExtractHists(dataEel_hist, mcEel_hists)
    arrayofhists = HistListtoArray(listofhists)
    sideband_list.append(arrayofhists)
bigArray = np.array(sideband_list)
matrices,vectors = SetupMatrices(bigArray)
for svd in SVD(matrices):
    print(svd[0],".\n",svd[1],".\n",svd[2],"\n\n")
weights = DoLinearAlgebra(matrices,vectors)
np.savetxt("/minerva/app/users/ajball/ncdiff/CC-NuE-XSec/tools/weights.txt",weights,delimiter=",")
print(weights)

np.set_printoptions(precision=0,suppress = True)
print(matrices,vectors)

mc_file.Close()
data_file.Close()

