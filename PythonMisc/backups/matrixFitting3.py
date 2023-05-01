# Alex Ball
# Matrix Unfolding with Exponential Fit for Diffractive Pion Production

import ROOT
import numpy as np
import math
import PlotUtils
from collections import OrderedDict
#from config.AnalysisConfig import AnalysisConfig

ROOT.TH1.AddDirectory(False)
# Classes to convert ROOT histograms into numpy arrays with values, error, and binning
class Variable1D:
    def __init__(self,name,filePath,skipFirst = False,isEfficiency = False):
        self.name = name
        self.values,self.error,self.bins = self.GetCVValsErrorandBins(name,filePath,skipFirst)
        self.isEfficiency = isEfficiency
        self.univvalues,self.univerror = self.GetUnivValsandError(name,filePath,skipFirst)

    def GetCVValsErrorandBins(self,name,filePath,skipFirst):
        dataFile = ROOT.TFile.Open(filePath,"READ")
        dataHist = eval("dataFile."+str(name)+".Clone()")
        self.cvhist = dataHist.Clone()
        dataArray = np.zeros(dataHist.GetNbinsX()-skipFirst)
        errorArray = np.zeros(dataHist.GetNbinsX()-skipFirst)
        dataBins = np.zeros(len(dataArray)+1)
        for i in range(len(dataArray)):
            dataArray[i] = dataHist.GetBinContent(i+1+skipFirst)
            errorArray[i] = dataHist.GetBinError(i+1+skipFirst)
            dataBins[i] = dataHist.GetBinLowEdge(i+1+skipFirst)
        dataBins[i+1] = dataHist.GetBinLowEdge(i+2+skipFirst) # Gets upper edge of last bin
        return dataArray,errorArray,dataBins

    def GetUnivValsandError(self,name,filePath,skipFirst):
        dataFile = ROOT.TFile.Open(filePath,"READ")
        dataHist = eval("dataFile."+str(name)+".Clone()")
        univvalues = {}
        univerrors = {}
        for error_type in dataHist.GetErrorBandNames():
            error_type_values_list = []
            error_type_errors_list = []
            for k in range(dataHist.GetVertErrorBand(error_type).GetNHists()):
                dataArray = np.zeros(dataHist.GetNbinsX()-skipFirst)
                errorArray = np.zeros(dataHist.GetNbinsX()-skipFirst)
                for i in range(len(dataArray)):
                    dataArray[i] = dataHist.GetVertErrorBand(error_type).GetHist(k).GetBinContent(i+1+skipFirst)
                    errorArray[i] = dataHist.GetVertErrorBand(error_type).GetHist(k).GetBinError(i+1+skipFirst)
                error_type_values_list.append(dataArray)
                error_type_errors_list.append(errorArray)
            univvalues[str(error_type)] = np.array(error_type_values_list)
            univerrors[str(error_type)] = np.array(error_type_errors_list)
        return univvalues,univerrors

class MigrationMatrix:
    def __init__(self,name,filePath,skipFirstX = False):
        self.name = name
        self.values,self.error,self.xbins,self.ybins = self.GetCVValsErrorandBins(name,filePath,skipFirstX)

    def GetCVValsErrorandBins(self,name,filePath,skipFirstX):
        dataFile = ROOT.TFile.Open(filePath,"READ")
        dataHist = eval("dataFile."+str(name)+".Clone()")
        self.cvhist = dataHist.Clone()
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

#Math for the matrix stuff
def simpson_integration(x0,xf,N,integrand): # Defines a generalized Simpson's rule algorithm that takes an integrand function
    xvals = np.linspace(x0,xf,N+1,endpoint=True) # Finds evenly spaced x values
    yvals = []
    for xval in xvals: # Creates a list of y values for each x value
        yvals.append(integrand(xval))
    h = (xf-x0)/N # Finds the distance between each x value
    integral = h/3*(yvals[0]+yvals[-1]+4*sum(yvals[1:N:2])+2*sum(yvals[2:N:2])) # Finds the integral with Simpson's rule
    return integral

def dN_dt(N,b,t):
    return N/b*np.exp(-t/b)

def y_bin_maker(N,b,migration,efficiency): # r_vec is [N,b]
    dist = lambda t: dN_dt(N,b,t) # The t distribution function
    true_y_bin_vals = np.zeros(len(migration.values))
    for i in range(len(migration.values)):
        true_y_bin_vals[i] = simpson_integration(migration.ybins[i],migration.ybins[i+1],100,dist)
    reco_y_bin_vals = true_y_bin_vals * efficiency.values
    return reco_y_bin_vals

def x_bin_maker(N,b,migration,efficiency): # r_vec is [N,b]
    y_bin_vals = y_bin_maker(N,b,migration,efficiency)
    x_bin_vals = np.dot(migration.values,y_bin_vals)
    return x_bin_vals

class Fitter(ROOT.Math.Functor):
    def Populate(self,data,migration,efficiency):
        self.data = data
        self.migration = migration
        self.efficiency = efficiency
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

    def y_bin_maker(self,N,b): # r_vec is [N,b]
        dist = lambda t: dN_dt(N,b,t) # The t distribution function
        true_y_bin_vals = np.zeros(len(self.migration.values))
        for i in range(len(self.migration.values)):
            true_y_bin_vals[i] = simpson_integration(self.migration.ybins[i],self.migration.ybins[i+1],100,dist)
        reco_y_bin_vals = true_y_bin_vals * self.efficiency.values
        return reco_y_bin_vals

    def NDim(self): return 2

def RunMinimizer(data,migration,efficiency):
    minimizer = ROOT.Math.Factory.CreateMinimizer("Minuit2")
    minimizer.SetTolerance(0.0001)
    fitter = Fitter()
    fitter.Populate(data,migration,efficiency)
    Nstart = 1000
    bstart = 0.01
    startingStep = 0.01
    minimizer.SetFunction(fitter)
    minimizer.SetVariable(0,"N",Nstart,startingStep)
    minimizer.SetVariable(1,"b",bstart,startingStep)
    minimizer.SetVariableLowerLimit(0,0.0)
    minimizer.SetVariableLowerLimit(1,0.0)
    minimizer.Minimize()
    print(minimizer.Correlation(0,1))
    result = [minimizer.X()[0],minimizer.X()[1]]
    error = [minimizer.Errors()[0],minimizer.Errors()[1]]
    return result,error

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

'''print(t_efficiency.values,t_efficiency.error,t_efficiency.bins)
print(IUE.values,IUE.error,IUE.bins)
print(t_migration.values,t_migration.error,t_migration.xbins,t_migration.ybins)'''
print(IUE.univvalues["Flux"][56])
result,error = RunMinimizer(IUE,t_migration,t_efficiency)
print(result,error)
print(x_bin_maker(result[0],result[1],t_migration,t_efficiency))
print(IUE.values)


