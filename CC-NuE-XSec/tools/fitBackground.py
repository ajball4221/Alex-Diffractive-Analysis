import ROOT
import PlotUtils
import math
from array import array
import numpy as np
import os

from tools.PlotLibrary import HistHolder
from config.AnalysisConfig import AnalysisConfig
from config.PlotConfig import ELECTRON_ENERGY_BINNING

ROOT.TH1.AddDirectory(True) # If you make this false, the error messages go away but the fitter just stops working because it depends on weird ROOT behavior with the POT scaling applying in ways I don't understand but it works as is
# If you really want to make this false, you gotta go back and rethink how the fitter grabs the root files and turn on the POT scaling. If you don't fix how it grabs ROOT files after making it false, this code will take up all the RAM on your system. 

def getPOTFromFile(filename):
    metatree = ROOT.TChain("Meta")
    if metatree.Add(filename,-1):
        return ROOT.PlotUtils.POTCounter().getPOTfromTChain(metatree)
    else:
        return None

def GetPOTScale(data_path,mc_path): # This is to correct for the fact that the mc and data have different simulated numbers of target protons in the neutrino maker
    pots = [None,None]
    for i,t in enumerate(["data","mc"]):
        path = [data_path,mc_path][i]
        try:
            pots[i]= getPOTFromFile(path) or getPOT(playlist,t,ntuple_tag)
        except KeyError:
            pots[i]=None
    pot_scale = pots[0]/pots[1] if pots.count(None) == 0 else 1.0
    print("POT Scale:",pot_scale, "Data POT:",pots[0],"Sim POT:",pots[1])
    return pot_scale

def AddTwoHistograms(band1,band2,pot_scale): # Adds band2 to band 1 and adds the pot scaled errors in quadrature
    errors = []
    for i in range(band1.GetSize()): # Adds errors in quadrature
        errors.append(math.sqrt(math.sqrt(pot_scale)*band1.GetBinError(i)**2 + math.sqrt(pot_scale)*band2.GetBinError(i)**2))
    band1.Add(band2,1) # Adds band2 to band1
    for i in range(band1.GetSize()): # Sets the error of signal to be what was calculated before
        band1.SetBinError(i,errors[i])
    return band1

def SubtractPoissonHistograms(signal,background_band,pot_scale,first_subtraction = False): # Sets error of signal where signal - background_band is being calculat
    errors = []
    for i in range(signal.GetSize()): # Adds errors in quadrature
        if first_subtraction: errors.append(math.sqrt(math.sqrt(pot_scale)*signal.GetBinError(i)**2 + math.sqrt(pot_scale)*background_band.GetBinError(i)**2)) # Scales the mc total error by the POT scale just once
        else: errors.append(math.sqrt(signal.GetBinError(i)**2 + math.sqrt(pot_scale)*background_band.GetBinError(i)**2))
    signal.Add(background_band,-1) # Subtracts the band from signal
    for i in range(signal.GetSize()): # Sets the error of signal to be what was calculated before
        signal.SetBinError(i,errors[i])
    return signal

def PrimeExtractHists(data_hists, mc_hists):
    pot_scale = mc_hists.pot_scale
    data = data_hists.GetHist().Clone()
    bands = ["NCCohPi0","NCDiff","CCQElike","notCCQElike","NCPi0","CCPi0"]
    coherent = mc_hists.hists["NCCohPi0"].Clone()
    ncdiff = mc_hists.hists["NCDiff"].Clone()
    nu_e = AddTwoHistograms(mc_hists.hists["CCQElike"].Clone(),mc_hists.hists["notCCQElike"].Clone(),pot_scale)
    piZero = AddTwoHistograms(mc_hists.hists["NCPi0"].Clone(),mc_hists.hists["CCPi0"].Clone(),pot_scale)
    total_mc = mc_hists.hists["Total"].Clone()
    for i,band in enumerate(bands): # Here, each of the subbands that are being calculated are removed from the total
        total_mc = SubtractPoissonHistograms(total_mc,mc_hists.hists[band],pot_scale,not bool(i)) 
    augment = SubtractPoissonHistograms(data,total_mc,pot_scale) # This subtracts the LHS const from the data to give the RHS of the equation
    return [coherent,nu_e,piZero,ncdiff,augment]

def ErrExtractHists(data_hists, mc_hists, err_type, universe_num):
    pot_scale = mc_hists.pot_scale
    data = data_hists.GetHist().Clone()
    coherent = mc_hists.hists["NCCohPi0"].GetVertErrorBand(err_type).GetHist(universe_num).Clone()
    ncdiff = mc_hists.hists["NCDiff"].GetVertErrorBand(err_type).GetHist(universe_num).Clone()
    nu_e = AddTwoHistograms(mc_hists.hists["CCQElike"].GetVertErrorBand(err_type).GetHist(universe_num).Clone(),mc_hists.hists["notCCQElike"].GetVertErrorBand(err_type).GetHist(universe_num).Clone(),pot_scale)
    piZero = AddTwoHistograms(mc_hists.hists["NCPi0"].GetVertErrorBand(err_type).GetHist(universe_num).Clone(),mc_hists.hists["CCPi0"].GetVertErrorBand(err_type).GetHist(universe_num).Clone(),pot_scale)
    total_mc = mc_hists.hists["Total"].GetVertErrorBand(err_type).GetHist(universe_num).Clone()
    bands = ["NCCohPi0","NCDiff","CCQElike","notCCQElike","NCPi0","CCPi0"]
    for i,band in enumerate(bands): # Here, each of the subbands that are being calculated are removed from the total
        #print(band,universe_num)
        #mc_hists.hists[band].GetVertErrorBand(err_type).GetHist(universe_num).Clone().Print("all")
        total_mc = SubtractPoissonHistograms(total_mc,mc_hists.hists[band].GetVertErrorBand(err_type).GetHist(universe_num),pot_scale,not bool(i))
    other = PlotUtils.MnvH1D(total_mc) # I have no idea why the add command can't do this itself but whatever
    augment = SubtractPoissonHistograms(data,other,pot_scale) # This subtracts the LHS const from the data to give the RHS of the equation
    return [coherent,nu_e,piZero,ncdiff,augment]

#def RejectSmallBins(hist): # This shouldnt be needed
#    for i in range(hist.GetSize()):
#        if hist.GetBinContent(i) < 1e-3:
#            hist.SetBinContent(i,0)
#    return hist

def HistListtoArray(histlist): # This converts the root histograms being passed to it into an array
    histarray = np.zeros((len(histlist),histlist[0].GetSize()-1))
    for i,hist in enumerate(histlist):
        for j in range(hist.GetSize()-1):
            histarray[i,j] = hist.GetBinContent(j+1)
    if len(histlist) == 1:
        histarray = histarray.flatten()
    return histarray

def SetupMatrices(dataarray): # This converts the weird data array that has a weird data structure into matrices that get the weights for each bin
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

def NormalInverse(matrices): #This shouldn't be called anywhere it's useful for testing how much the SVD made things weird
    inverses = []
    for matrix in matrices:
        inverse = np.linalg.inv(matrix)
        inverses.append(inverse)
    return inverses

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
    # These may look different for you and you may need to rethink these
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

def GetPreweightedROOTFiles(weighting_playlist):
    if weighting_playlist == None: #If not specified, use the default weighting files
        mc_path = str(os.environ["CCNUEROOT"]+"/tools/kin_dist_mcAlex_prenoq3_nx_collab1_fspline.root")
        mc_file = ROOT.TFile.Open(mc_path)
        data_path = str(os.environ["CCNUEROOT"]+"/tools/kin_dist_dataAlex_prenoq3_nx_collab1_fspline.root")
        data_file = ROOT.TFile.Open(data_path)
        pot_scale = GetPOTScale(data_path,mc_path)
    else:
        mc_path = "/minerva/data/users/ajball/nu_e/kin_dist_mc"+str(weighting_playlist)+"_collab1_fspline.root"
        mc_file = ROOT.TFile.Open(mc_path)
        data_path = "/minerva/data/users/ajball/nu_e/kin_dist_data"+str(weighting_playlist)+"_collab1_fspline.root"
        data_file = ROOT.TFile.Open(data_path)
        pot_scale = GetPOTScale(data_path,mc_path)
    return mc_file,data_file,pot_scale

def GetCVWeights(mc_file,data_file,pot_scale):
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
    WasIRight(matrices,primeweights,vectors)
    return primeweights

def WasIRight(matrices,weights,vectors):
    for i in range(len(matrices)):
        print("MC: ",np.dot(matrices[i],weights[i]-np.ones(4)))
        print("Data: ",vectors[i])

def GetSystematicUniverseWeights(mc_file,data_file,pot_scale):
    # Now for the Systematic Universes
    sidebands = ["Coherent","Electron_Neutrino","Pi0","Signal"]
    err_weights_dict = {}
    universe_list = [["Flux",200],["GENIE_AGKYxF1pi",2],["GENIE_AhtBY",2],["GENIE_BhtBY",2],["GENIE_CCQEPauliSupViaKF",2],["GENIE_CV1uBY",2],["GENIE_CV2uBY",2],["GENIE_EtaNCEL",2],["GENIE_FrAbs_N",2],["GENIE_FrAbs_pi",2],["GENIE_FrCEx_N",2],["GENIE_FrCEx_pi",2],["GENIE_FrElas_N",2],["GENIE_FrElas_pi",2],["GENIE_FrInel_N",2],["GENIE_FrPiProd_N",2],["GENIE_FrPiProd_pi",2],["GENIE_MFP_N",2],["GENIE_MFP_pi",2],["GENIE_MaNCEL",2],["GENIE_MaRES",2],["GENIE_MvRES",2],["GENIE_NormCCRES",2],["GENIE_NormDISCC",2],["GENIE_NormNCRES",2],["GENIE_RDecBR1gamma",2],["GENIE_Rvn2pi",2],["GENIE_Rvp2pi",2],["GENIE_Theta_Delta2Npi",2],["GENIE_VecFFCCQEshape",2],["GENIE_Rvn1pi",2],["GENIE_Rvp1pi",2],["GENIE_MaCCQE",2],["Low_Recoil_2p2h_Tune",3],["RPA_HighQ2",2],["RPA_LowQ2",2],["LowQ2Pi",2],["LowQ2Pi_None",1],["MK_model",1],["fsi_weight",3],["SuSA_Valencia_Weight",1],["GEANT_Proton",2],["GEANT_Neutron",2],["GEANT_Pion",2],["bkg_tune",2],["Target_Mass_CH",2],["Muon_Energy_MINERvA",2],["Muon_Energy_MINOS",2],["beam_angle",4],["response_p",2],["response_meson",2],["response_em",2],["response_other",2],["response_xtalk",2],["Leakage_Uncertainty",2]]
    for error_type,error_universes in universe_list:
        err_type_weights = []
        for universe_num in range(error_universes):
            sideband_list = []
            for sideband in sidebands:
                mcEel_hists = HistHolder("Lepton Energy",mc_file,sideband,is_mc = True,pot = pot_scale) 
                dataEel_hist = HistHolder("Lepton Energy",data_file,sideband,is_mc = False)
                #mcEel_hists.POTScale(False)
                #dataEel_hist.POTScale(False)
                listofhists = ErrExtractHists(dataEel_hist, mcEel_hists, error_type, universe_num)
                arrayofhists = HistListtoArray(listofhists)
                sideband_list.append(arrayofhists)
            bigArray = np.array(sideband_list)
            matrices,vectors = SetupMatrices(bigArray)
            err_type_weights.append(DoLinearAlgebra(MPInverse(matrices),vectors))
            del mcEel_hists, dataEel_hist
        err_weights_dict[error_type] = np.array(err_type_weights)
    return err_weights_dict

def RunAlexWeighter():
    print("Starting Alex's Reweighting Weighter")
    mc_file,data_file,pot_scale = GetPreweightedROOTFiles(None)
    cvweights = GetCVWeights(mc_file,data_file,pot_scale)
    universe_weights_dict = GetSystematicUniverseWeights(mc_file,data_file,pot_scale)
    #mc_file.Close()
    #data_file.Close()
    print("Weighter Ran Successfully")
    return cvweights,universe_weights_dict

if __name__ == "__main__":
    preweighted_playlist = AnalysisConfig.playlist
    mc_file,data_file,pot_scale = GetPreweightedROOTFiles(preweighted_playlist)
    cvweights = GetCVWeights(mc_file,data_file,pot_scale)
    universe_weights_dict = GetSystematicUniverseWeights(mc_file,data_file,pot_scale)
    np.set_printoptions(precision=6,suppress = True,threshold = np.inf)
    print("Weights for Main Universe:\n",cvweights,"\nSystematics Universes Error Weights:\n",universe_weights_dict)
    Eel_hist = mc_file.Eel.Clone()
    #mc_file.Close()
    #data_file.Close()
    weightFile = ROOT.TFile.Open("datafolder/nu_e/"+str(preweighted_playlist)+"_weights.root","RECREATE")
    coherent = PlotUtils.MnvH1D("coherent","Weights on Coherent Events",len(ELECTRON_ENERGY_BINNING)-1,array("d",ELECTRON_ENERGY_BINNING))
    electron_neutrino = PlotUtils.MnvH1D("electron_neutrino","Weights on Electron Neutrino Events",len(ELECTRON_ENERGY_BINNING)-1,array("d",ELECTRON_ENERGY_BINNING))
    pizero = PlotUtils.MnvH1D("pizero","Weights on PiZero Events",len(ELECTRON_ENERGY_BINNING)-1,array("d",ELECTRON_ENERGY_BINNING))
    ncdiff = PlotUtils.MnvH1D("ncdiff","Weights on NCDiff Events",len(ELECTRON_ENERGY_BINNING)-1,array("d",ELECTRON_ENERGY_BINNING))
    for i in range(len(ELECTRON_ENERGY_BINNING)-1):
        coherent.SetBinContent(i+1,cvweights[i,0])
        electron_neutrino.SetBinContent(i+1,cvweights[i,1])
        pizero.SetBinContent(i+1,cvweights[i,2])
        ncdiff.SetBinContent(i+1,cvweights[i,3])
    coherent.AddMissingErrorBandsAndFillWithCV(Eel_hist)
    electron_neutrino.AddMissingErrorBandsAndFillWithCV(Eel_hist)
    pizero.AddMissingErrorBandsAndFillWithCV(Eel_hist)
    ncdiff.AddMissingErrorBandsAndFillWithCV(Eel_hist)
    for error_type in coherent.GetErrorBandNames():
        error_type = str(error_type) # I really hate that this line has to be here
        for k in range(coherent.GetVertErrorBand(error_type).GetNHists()):
            for i in range(len(ELECTRON_ENERGY_BINNING)-1):
                coherent.GetVertErrorBand(error_type).GetHist(k).SetBinContent(i+1,universe_weights_dict[error_type][k][i,0])
                electron_neutrino.GetVertErrorBand(error_type).GetHist(k).SetBinContent(i+1,universe_weights_dict[error_type][k][i,1])
                pizero.GetVertErrorBand(error_type).GetHist(k).SetBinContent(i+1,universe_weights_dict[error_type][k][i,2])
                ncdiff.GetVertErrorBand(error_type).GetHist(k).SetBinContent(i+1,universe_weights_dict[error_type][k][i,3])
    coherent.GetXaxis().SetTitle("Lepton Energy")
    coherent.SetMarkerStyle(21)
    electron_neutrino.GetXaxis().SetTitle("Lepton Energy")
    electron_neutrino.SetMarkerStyle(21)
    pizero.GetXaxis().SetTitle("Lepton Energy")
    pizero.SetMarkerStyle(21)
    ncdiff.GetXaxis().SetTitle("Lepton Energy")
    ncdiff.SetMarkerStyle(21)
    coherent.Write()
    electron_neutrino.Write()
    pizero.Write()
    ncdiff.Write()
    weightFile.Close()
