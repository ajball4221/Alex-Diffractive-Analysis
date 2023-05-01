# Alex Ball
# Matrix Unfolding with Exponential Fit for Diffractive Pion Production

import ROOT
import numpy as np
import math
import PlotUtils
from collections import OrderedDict
from array import array

#from config.AnalysisConfig import AnalysisConfig

ROOT.TH1.AddDirectory(False)
# Classes to convert ROOT histograms into numpy arrays with values, error, and binning
class Variable1D:
    def __init__(self,name,filePath,skipFirst = False,isEfficiency = False,univ_type = None,univ_num = None):
        self.name = name
        self.isCV = not bool(univ_type)
        self.isEfficiency = isEfficiency
        if self.isCV:
            self.values,self.error,self.bins = self.GetCVValsErrorandBins(name,filePath,skipFirst)
            self.universes = self.UniverseListAssembler(name,filePath,skipFirst)
        else: 
            self.univ_type = univ_type
            self.univ_num = univ_num
            self.values,self.error,self.bins = self.GetUnivValsErrorandBins(name,filePath,skipFirst,univ_type,univ_num)
        #if self.isEfficiency: self.values = np.ones(len(self.values)) # For testing with no efficiency

    def GetCVValsErrorandBins(self,name,filePath,skipFirst):
        dataFile = ROOT.TFile.Open(filePath,"READ")
        print(dataFile)
        dataHist = eval("dataFile."+str(name)+".Clone()")
        self.hist = dataHist.Clone()
        dataArray = np.zeros(dataHist.GetNbinsX()-skipFirst)
        errorArray = np.zeros(dataHist.GetNbinsX()-skipFirst)
        dataBins = np.zeros(len(dataArray)+1)
        for i in range(len(dataArray)):
            dataArray[i] = dataHist.GetBinContent(i+1+skipFirst)
            errorArray[i] = dataHist.GetBinError(i+1+skipFirst)
            dataBins[i] = dataHist.GetBinLowEdge(i+1+skipFirst)
        dataBins[i+1] = dataHist.GetBinLowEdge(i+2+skipFirst) # Gets upper edge of last bin
        return dataArray,errorArray,dataBins

    def UniverseListAssembler(self,name,filePath,skip1st):
        dataFile = ROOT.TFile.Open(filePath,"READ")
        dataHist = eval("dataFile."+str(name)+".Clone()")
        universedict = {}
        for error_type in dataHist.GetErrorBandNames():
            error_type_univ_list = []
            for k in range(dataHist.GetVertErrorBand(error_type).GetNHists()):
                error_type_univ_list.append(Variable1D(name,filePath,skipFirst = skip1st, isEfficiency = self.isEfficiency,univ_type = error_type,univ_num = k))
            universedict[str(error_type)] = error_type_univ_list
        return universedict

    def GetUnivValsErrorandBins(self,name,filePath,skipFirst,error_type,k):
        dataFile = ROOT.TFile.Open(filePath,"READ")
        dataHist = eval("dataFile."+str(name)+".Clone()")
        dataArray = np.zeros(dataHist.GetVertErrorBand(error_type).GetHist(k).GetNbinsX()-skipFirst)
        errorArray = np.zeros(dataHist.GetVertErrorBand(error_type).GetHist(k).GetNbinsX()-skipFirst)
        dataBins = np.zeros(len(dataArray)+1)
        for i in range(len(dataArray)):
            dataArray[i] = dataHist.GetVertErrorBand(error_type).GetHist(k).GetBinContent(i+1+skipFirst)
            errorArray[i] = dataHist.GetVertErrorBand(error_type).GetHist(k).GetBinError(i+1+skipFirst)
            dataBins[i] = dataHist.GetVertErrorBand(error_type).GetHist(k).GetBinLowEdge(i+1+skipFirst)
        dataBins[i+1] = dataHist.GetVertErrorBand(error_type).GetHist(k).GetBinLowEdge(i+2+skipFirst) # Gets upper edge of last bin
        return dataArray,errorArray,dataBins

class MigrationMatrix:
    def __init__(self,name,filePath,skipFirstX = False, univ_type = None,univ_num = None):
        self.name = name
        self.isCV = not bool(univ_type)
        if self.isCV:
            self.values,self.error,self.xbins,self.ybins = self.GetCVValsErrorandBins(name,filePath,skipFirstX)
            self.universes = self.UniverseListAssembler(name,filePath,skipFirstX)
        else:
            self.univ_type = univ_type
            self.univ_num = univ_num
            self.values,self.error,self.xbins,self.ybins = self.GetUnivValsErrorandBins(name,filePath,skipFirstX,univ_type,univ_num)

    def GetCVValsErrorandBins(self,name,filePath,skipFirstX):
        dataFile = ROOT.TFile.Open(filePath,"READ")
        dataHist = eval("dataFile."+str(name)+".Clone()")
        self.hist = dataHist.Clone()
        dataArray = np.zeros((dataHist.GetNbinsY(),dataHist.GetNbinsX()-skipFirstX))
        errorArray = np.zeros((dataHist.GetNbinsY(),dataHist.GetNbinsX()-skipFirstX))
        for i in range(len(dataArray)):
            for j in range(len(dataArray[0])):
                dataArray[i,j] = dataHist.GetBinContent(j+1+skipFirstX,i+1)
                errorArray[i,j] = dataHist.GetBinError(j+1+skipFirstX,i+1)
            errorArray[i] /= np.sum(dataArray[i])
            dataArray[i] /= np.sum(dataArray[i])
        ybins = np.zeros(len(dataArray)+1)
        xbins = np.zeros(len(dataArray[0])+1)
        yaxis = dataHist.GetYaxis()
        xaxis = dataHist.GetXaxis()
        for i in range(len(ybins)):
            ybins[i] = yaxis.GetBinLowEdge(i+1)
        for j in range(len(xbins)):
            xbins[j] = xaxis.GetBinLowEdge(j+1+skipFirstX)
        #for i in range(len(xbins)-1): # this chunk is for xbin normalized migration
        #    dataArray[:,i] /= (xbins[i+1]-xbins[i])
        #for i in range(len(dataArray)): dataArray[i] /= np.sum(dataArray[i])
        return dataArray,errorArray,xbins,ybins

    def UniverseListAssembler(self,name,filePath,skip1stX):
        dataFile = ROOT.TFile.Open(filePath,"READ")
        dataHist = eval("dataFile."+str(name)+".Clone()")
        universedict = {}
        for error_type in dataHist.GetErrorBandNames():
            error_type_univ_list = []
            for k in range(dataHist.GetVertErrorBand(error_type).GetNHists()):
                error_type_univ_list.append(MigrationMatrix(name,filePath,skipFirstX = skip1stX,univ_type = error_type,univ_num = k))
            universedict[str(error_type)] = error_type_univ_list
        return universedict

    def GetUnivValsErrorandBins(self,name,filePath,skipFirstX,error_type,k):
        dataFile = ROOT.TFile.Open(filePath,"READ")
        dataHist = eval("dataFile."+str(name)+".Clone()")
        dataArray = np.zeros((dataHist.GetNbinsY(),dataHist.GetNbinsX()-skipFirstX))
        errorArray = np.zeros((dataHist.GetNbinsY(),dataHist.GetNbinsX()-skipFirstX))
        for i in range(len(dataArray)):
            for j in range(len(dataArray[0])):
                dataArray[i,j] = dataHist.GetVertErrorBand(error_type).GetHist(k).GetBinContent(j+1+skipFirstX,i+1)
                errorArray[i,j] = dataHist.GetVertErrorBand(error_type).GetHist(k).GetBinError(j+1+skipFirstX,i+1)
            errorArray[i] /= np.sum(dataArray[i])
            dataArray[i] /= np.sum(dataArray[i])
        ybins = np.zeros(len(dataArray)+1)
        xbins = np.zeros(len(dataArray[0])+1)
        yaxis = dataHist.GetYaxis()
        xaxis = dataHist.GetXaxis()
        for i in range(len(ybins)):
            ybins[i] = yaxis.GetBinLowEdge(i+1)
        for j in range(len(xbins)):
            xbins[j] = xaxis.GetBinLowEdge(j+1+skipFirstX)
        #for i in range(len(xbins)-1): # this chunk is for xbin normalized migration
        #    dataArray[:,i] /= (xbins[i+1]-xbins[i])
        #for i in range(len(dataArray)): dataArray[i] /= np.sum(dataArray[i])
        return dataArray,errorArray,xbins,ybins

def dN_dt(N,b,t):
    return N/b*np.exp(-t/b)

def IntegralofdN_dt(N,b,ti,tf):
    return N*(np.exp(-ti/b)-np.exp(-tf/b))

def pres_format(value,uncertainty): # Only works for floats
    place = int(math.floor(math.log10(uncertainty))) - 1
    final_value = round(value,-place)
    final_uncertainty = round(uncertainty,-place)
    if final_value == int(final_value) and final_uncertainty == int(final_uncertainty): 
        final_value = int(final_value)
        final_uncertainty = int(final_uncertainty)
        return str(final_value) +" ± "+ str(final_uncertainty)
    if str(final_uncertainty)[-2] == "0":
        extrauncert0 = "0"
    else: extrauncert0 = ""
    if len((str(final_value).split("."))[1]) != len((str(final_uncertainty).split("."))[1])+len(extrauncert0):
        zerostoadd = str("0"*(len((str(final_uncertainty).split("."))[1])-(len((str(final_value).split("."))[1])+len(extrauncert0))))
        return str(final_value) + zerostoadd + " ± "+ str(final_uncertainty) + extrauncert0
    return str(final_value) +" ± "+ str(final_uncertainty)

class Fitter(ROOT.Math.Functor):
    def Populate(self,data,migration,efficiency):
        self.data = data
        self.migration = migration
        self.efficiency = efficiency
        print(self.efficiency.values)
        self.populated = True

    def DoEval(self,r_vec):
        return self.CalcChi2(r_vec[0],r_vec[1])

    def CalcChi2(self,N,b):
        fit = self.x_bin_maker(N,b)
        total_variance = 0
        for i in range(len(fit)):
            total_variance += (self.data.values[i]-fit[i])**2/((self.data.error[i])**2)
        return total_variance

    def x_bin_maker(self,N,b): # r_vec is [N,b]
        y_bin_vals = self.y_bin_maker(N,b)
        x_bin_vals = np.dot(self.migration.values,y_bin_vals)
        return x_bin_vals

    def y_bin_maker(self,N,b,reco = True): # r_vec is [N,b]
        true_y_bin_vals = np.zeros(len(self.migration.values))
        for i in range(len(self.migration.values)):
            true_y_bin_vals[i] = IntegralofdN_dt(N,b,self.migration.ybins[i],self.migration.ybins[i+1])
        if reco:
            reco_y_bin_vals = true_y_bin_vals * self.efficiency.values
            return reco_y_bin_vals
        else: return true_y_bin_vals
    
    def GetFittedBinErrorX(self,N,b,dN,db,rho): # Finds error in truth space based on exponential distribution
        jacobian = np.zeros((len(self.migration.values),2))
        for i in range(len(self.migration.values)):
            ti = self.migration.ybins[i] # ti and tf are the bin edges
            tf = self.migration.ybins[i+1]
            dF_dN = (np.exp(-ti/b)-np.exp(-tf/b))*self.efficiency.values[i]
            dF_db = (N*(ti*np.exp(-ti/b)-tf*np.exp(-tf/b))/(b**2))*self.efficiency.values[i]
            jacobian[i] = [dF_dN,dF_db]
        cov_matrix_param = np.array([[dN**2,rho*dN*db],[rho*dN*db,db**2]])
        cov_matrix_y = jacobian @ cov_matrix_param @ np.transpose(jacobian)
        cov_matrix_x = self.migration.values @ cov_matrix_y @ np.transpose(self.migration.values)
        return np.sqrt(np.diagonal(cov_matrix_x))

    def GetFittedBinErrorY(self,N,b,dN,db,rho): # Finds error in truth space based on exponential distribution 
        binerrors = np.zeros(len(self.migration.values))
        for i in range(len(self.migration.values)):
            ti = self.migration.ybins[i] # ti and tf are the bin edges
            tf = self.migration.ybins[i+1]
            dF_db = N*(ti*np.exp(-ti/b)-tf*np.exp(-tf/b))/(b**2)
            dF_dN = (np.exp(-ti/b)-np.exp(-tf/b))
            binerrors[i] = np.sqrt((dF_dN*dN)**2+(dF_db*db)**2+2*rho*(dF_dN*dN)*(dF_db*db))
        return binerrors

    def NDim(self): return 2

    def Reset(self):
        self.data = None
        self.migration = None 
        self.efficiency = None
        self.populated = False


def RunMinimizer(fitter):
    minimizer = ROOT.Math.Factory.CreateMinimizer("Minuit2")
    minimizer.SetTolerance(0.0001)
    Nstart = 1000
    bstart = 1
    startingStep = 0.01
    minimizer.SetFunction(fitter)
    minimizer.SetVariable(0,"N",Nstart,startingStep)
    minimizer.SetVariable(1,"b",bstart,startingStep)
    minimizer.SetVariableLowerLimit(0,0.0)
    minimizer.SetVariableLowerLimit(1,0.0)
    minimizer.Minimize()
    correlation = minimizer.Correlation(0,1)
    fit_result = [minimizer.X()[0],minimizer.X()[1]]
    fit_error = [minimizer.Errors()[0],minimizer.Errors()[1]]
    return fit_result,fit_error,correlation

def ExponentialTruthFit(data,migration,efficiency,xname,yname): # These have to be the main data classes
    fitter = Fitter()
    fitter.Populate(data,migration,efficiency) # Until reset, the fitter holds the relevant distributions and migration matrix
    result,error,rho = RunMinimizer(fitter)
    y_fitted_dist = fitter.y_bin_maker(result[0],result[1],reco = False)
    y_fitted_error = fitter.GetFittedBinErrorY(result[0],result[1],error[0],error[1],rho)
    x_predicted_dist = fitter.x_bin_maker(result[0],result[1])
    x_predicted_error = fitter.GetFittedBinErrorX(result[0],result[1],error[0],error[1],rho)

    y_fitted_hist = PlotUtils.MnvH1D(str(yname)+"_fitted","True "+str(yname)+" Fitted from Exponential",len(migration.values),array("d",migration.ybins))
    x_predicted_hist = PlotUtils.MnvH1D(str(xname)+"_fit_prediction",str(xname)+" Predicted from "+str(yname)+" Fit",len(migration.values[0]),array("d",migration.xbins))
    y_fitted_hist.AddMissingErrorBandsAndFillWithCV(efficiency.hist)
    x_predicted_hist.AddMissingErrorBandsAndFillWithCV(data.hist)
    for i in range(len(y_fitted_dist)):
        y_fitted_hist.SetBinContent(i+1,y_fitted_dist[i])
        y_fitted_hist.SetBinError(i+1,y_fitted_error[i])
    for j in range(len(x_predicted_dist)):
        x_predicted_hist.SetBinContent(j+1,x_predicted_dist[j])
        x_predicted_hist.SetBinError(j+1,x_predicted_error[j])
    print("N = "+pres_format(result[0],error[0])+", b = "+pres_format(result[1],error[1])+" cv")
    fitter.Reset()
    
    y_fitted_distcv = y_fitted_dist
    y_fitted_errorcv = y_fitted_error
    x_predicted_distcv = x_predicted_dist
    x_predicted_errorcv = x_predicted_error
    for error_type in data.hist.GetErrorBandNames():
        error_type = str(error_type)
        for i in range(len(y_fitted_distcv)): # This is dumb that it makes me do this JUST LOOK AT THE MAIN CV HIST TO CALC ERROR
            y_fitted_hist.GetVertErrorBand(error_type).SetBinContent(i+1,y_fitted_distcv[i])
            y_fitted_hist.GetVertErrorBand(error_type).SetBinError(i+1,y_fitted_errorcv[i])
        for j in range(len(x_predicted_distcv)):
            x_predicted_hist.GetVertErrorBand(error_type).SetBinContent(j+1,x_predicted_distcv[j])
            x_predicted_hist.GetVertErrorBand(error_type).SetBinError(j+1,x_predicted_errorcv[j])
        for k in range(data.hist.GetVertErrorBand(error_type).GetNHists()):
            fitter.Populate(data.universes[error_type][k],migration.universes[error_type][k],efficiency.universes[error_type][k])
            result,error,rho = RunMinimizer(fitter)
            y_fitted_dist = fitter.y_bin_maker(result[0],result[1],reco = False)
            y_fitted_error = fitter.GetFittedBinErrorY(result[0],result[1],error[0],error[1],rho)
            x_predicted_dist = fitter.x_bin_maker(result[0],result[1])
            x_predicted_error = fitter.GetFittedBinErrorX(result[0],result[1],error[0],error[1],rho)

            for i in range(len(y_fitted_dist)):
                y_fitted_hist.GetVertErrorBand(error_type).GetHist(k).SetBinContent(i+1,y_fitted_dist[i])
                y_fitted_hist.GetVertErrorBand(error_type).GetHist(k).SetBinError(i+1,y_fitted_error[i])
            for j in range(len(x_predicted_dist)):
                x_predicted_hist.GetVertErrorBand(error_type).GetHist(k).SetBinContent(j+1,x_predicted_dist[j])
                x_predicted_hist.GetVertErrorBand(error_type).GetHist(k).SetBinError(j+1,x_predicted_error[j])
            print("N = "+pres_format(result[0],error[0])+", b = "+pres_format(result[1],error[1]),error_type,k)
            fitter.Reset()
    return y_fitted_hist,x_predicted_hist

def InverseMigration(data,migration,efficiency):
    invmatrix = np.linalg.inv(migration.values)
    reco_y_dist = np.dot(invmatrix,data.values)
    true_y_dist = reco_y_dist / efficiency.values
    reco_y_error = np.zeros(len(migration.values[0]))
    for j in range(len(migration.values[0])):
       error_sum = 0
       for i in range(len(migration.values)):
           error_sum += (migration.values[j,i]**2)*data.values[i]
       reco_y_error[j] = np.sqrt(error_sum)
    true_y_error = reco_y_error / efficiency.values
    return true_y_dist, true_y_error
    
def MatrixTruthRecovery(data,migration,efficiency,xname,yname): # These have to be the CV classes that include the universes
    y_recovered_dist, y_recovered_error = InverseMigration(data,migration,efficiency)
    print(y_recovered_dist, yname)
    y_recovered_hist = PlotUtils.MnvH1D(str(yname)+"_recovered","True "+str(yname)+" From Inverse Migration Matrix",len(migration.values),array("d",migration.ybins))
    y_recovered_hist.AddMissingErrorBandsAndFillWithCV(efficiency.hist)
    for i in range(len(y_recovered_dist)):
        y_recovered_hist.SetBinContent(i+1,y_recovered_dist[i])
        y_recovered_hist.SetBinError(i+1,y_recovered_error[i])
    y_recovered_distcv = y_recovered_dist
    y_recovered_errorcv = y_recovered_error
    for error_type in data.hist.GetErrorBandNames():
        error_type = str(error_type) # I really hate that this line has to be here
        for i in range(len(y_recovered_distcv)): # This is dumb that it makes me do this JUST LOOK AT THE MAIN CV HIST TO CALC ERROR
            y_recovered_hist.GetVertErrorBand(error_type).SetBinContent(i+1,y_recovered_distcv[i])
            y_recovered_hist.GetVertErrorBand(error_type).SetBinError(i+1,y_recovered_errorcv[i])
        for k in range(data.hist.GetVertErrorBand(error_type).GetNHists()):
            y_recovered_dist, y_recovered_error = InverseMigration(data.universes[error_type][k],migration.universes[error_type][k],efficiency.universes[error_type][k])
            for i in range(len(y_recovered_dist)):
                y_recovered_hist.GetVertErrorBand(error_type).GetHist(k).SetBinContent(i+1,y_recovered_dist[i])
                y_recovered_hist.GetVertErrorBand(error_type).GetHist(k).SetBinError(i+1,y_recovered_error[i])
    return y_recovered_hist

playlist = "Alex_tcut_nx"
efficiencyPath = "/minerva/app/users/ajball/ncdiff/CC-NuE-XSec/efficiency/"+str(playlist)+"_calculatedEfficiency.root"
subtractedDataPath = "/minerva/data/users/ajball/nu_e/kin_dist_data"+str(playlist)+"_subtracted_collab1_fspline.root"
migrationPath = "/minerva/data/users/ajball/nu_e/kin_dist_mc"+str(playlist)+"_collab1_fspline.root"

t_efficiency = Variable1D("t_efficiency",efficiencyPath,isEfficiency = True)
IUE = Variable1D("InlineUpstream_Energy",subtractedDataPath,skipFirst = True)
t_migration = MigrationMatrix("inline_vs_t",migrationPath,skipFirstX = True)

pion_efficiency = Variable1D("pion_efficiency",efficiencyPath,isEfficiency = True)
Eel = Variable1D("Eel",subtractedDataPath)
pion_migration = MigrationMatrix("ee_vs_pionE",migrationPath)

t_fitted_hist, IUE_predicted_hist = ExponentialTruthFit(IUE,t_migration,t_efficiency,"IUE","t")
t_recovered_hist = MatrixTruthRecovery(IUE,t_migration,t_efficiency,"IUE","t")
pion_fitted_hist, Eel_predicted_hist = ExponentialTruthFit(Eel,pion_migration,pion_efficiency,"Eel","pionE")
pion_recovered_hist = MatrixTruthRecovery(Eel,pion_migration,pion_efficiency,"Eel","pionE")

outPath = "/minerva/data/users/ajball/nu_e/"+str(playlist)+"_unfolded.root"
outFile = ROOT.TFile.Open(outPath,"RECREATE")
t_fitted_hist.Write() 
IUE_predicted_hist.Write()
t_recovered_hist.Write()
pion_fitted_hist.Write()
Eel_predicted_hist.Write()
pion_recovered_hist.Write()
outFile.Close()

t_fitted_err = t_fitted_hist.GetCVHistoWithError()
t_fitted_err.SetMarkerStyle(20)
t_fitted_err.SetMarkerColor(2)
#t_fitted_err.SetMaximum(36000)
t_fitted_err.SetMaximum(600)
t_fitted_err.GetXaxis().SetTitle("|t| (GeV/c)^2")

t_recovered_err = t_recovered_hist.GetCVHistoWithError()
t_recovered_err.SetMarkerStyle(21)
t_recovered_err.SetMarkerColor(4)
t_fitted_err.Draw("MIN0 E1")
t_recovered_err.Draw("SAME MIN0 E1")
ROOT.gPad.BuildLegend()

pion_canvas = ROOT.TCanvas()
pion_fitted_err = pion_fitted_hist.GetCVHistoWithError()
pion_fitted_err.SetMarkerStyle(20)
pion_fitted_err.SetMarkerColor(2)
#pion_fitted_err.SetMaximum(7500)
pion_fitted_err.GetXaxis().SetTitle("Pion Energy (GeV)")

pion_recovered_err = pion_recovered_hist.GetCVHistoWithError()
pion_recovered_err.SetMarkerStyle(21)
pion_recovered_err.SetMarkerColor(4)
pion_fitted_err.Draw("MIN0 E1")
pion_recovered_err.Draw("SAME MIN0 E1")
ROOT.gPad.BuildLegend()

t_fit_accuracy_canvas = ROOT.TCanvas()
IUE_predicted_err = IUE_predicted_hist.GetCVHistoWithError()
IUE_predicted_err.SetMarkerStyle(20)
IUE_predicted_err.SetMarkerColor(2)
IUE_predicted_err.GetXaxis().SetTitle("Inline Upstream Energy (MeV)")

IUE_data_err = IUE.hist.GetCVHistoWithError()
IUE_data_err.SetMarkerStyle(21)
IUE_data_err.SetMarkerColor(4)
IUE_predicted_err.Draw("MIN0 E1")
IUE_data_err.Draw("SAME MIN0 E1")
ROOT.gPad.BuildLegend()

pion_fit_accuracy_canvas = ROOT.TCanvas()
Eel_predicted_err = Eel_predicted_hist.GetCVHistoWithError()
Eel_predicted_err.SetMarkerStyle(20)
Eel_predicted_err.SetMarkerColor(2)
Eel_predicted_err.GetXaxis().SetTitle("Electron Cone Energy (GeV)")

Eel_data_err = Eel.hist.GetCVHistoWithError()
Eel_data_err.SetMarkerStyle(21)
Eel_data_err.SetMarkerColor(4)
Eel_predicted_err.Draw("MIN0 E1")
Eel_data_err.Draw("SAME MIN0 E1")
ROOT.gPad.BuildLegend()

input("Press any key to exit")
