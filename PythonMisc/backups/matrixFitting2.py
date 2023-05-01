# Alex Ball
# Matrix Unfolding with Exponential Fit for Diffractive Pion Production

import ROOT
import numpy as np
import math
import PlotUtils
from config.AnalysisConfig import AnalysisConfig

class Variable1D:
    def __init__(self,name,filePath,skipFirst = False,isEfficiency = False):
        self.name = name
        self.values,self.error,self.bins = self.GetCVValsErrorandBins(name,filePath,skipFirst)
        self.isEfficiency = isEfficiency

    def GetCVValsErrorandBins(self,name,filePath,skipFirst):
        dataFile = ROOT.TFile.Open(filePath,"READ")
        dataHist = eval("dataFile."+str(name)+".Clone()")
        dataArray = np.zeros(dataHist.GetNbinsX()-skipFirst)
        errorArray = np.zeros(dataHist.GetNbinsX()-skipFirst)
        dataBins = np.zeros(len(dataArray)+1)
        for i in range(len(dataArray)):
            dataArray[i] = dataHist.GetBinContent(i+1+skipFirst)
            errorArray[i] = dataHist.GetBinError(i+1+skipFirst)
            dataBins[i] = dataHist.GetBinLowEdge(i+1+skipFirst)
        dataBins[i+1] = dataHist.GetBinLowEdge(i+2+skipFirst) # Gets upper edge of last bin
        return dataArray,errorArray,dataBins

class MigrationMatrix:
    def __init__(self,name,filePath,skipFirstX = False):
        self.name = name
        self.values,self.error,self.xbins,self.ybins = self.GetCVValsErrorandBins(name,filePath,skipFirstX)

    def GetCVValsErrorandBins(self,name,filePath,skipFirstX):
        dataFile = ROOT.TFile.Open(filePath,"READ")
        dataHist = eval("dataFile."+str(name)+".Clone()")
        dataArray = np.zeros((dataHist.GetNbinsY(),dataHist.GetNbinsX()-skipFirstX))
        errorArray = np.zeros((dataHist.GetNbinsY(),dataHist.GetNbinsX()-skipFirstX))
        for i in range(len(dataArray)):
            for j in range(len(dataArray[0])):
                dataArray[i,j] = dataHist.GetBinContent(j+1+skipFirstX,i+1)
                errorArray[i,j] = dataHist.GetBinError(j+1+skipFirstX,i+1)
            dataArray[i] /= np.sum(dataArray[i])
            errorArray[i] /= np.sum(dataArray[i])
        ybins = np.zeros(len(dataArray)+1)
        xbins = np.zeros(len(dataArray[0])+1)
        yaxis = dataHist.GetYaxis()
        xaxis = dataHist.GetXaxis()
        for i in range(len(ybins)):
            ybins[i] = yaxis.GetBinLowEdge(i+1)
        for j in range(len(xbins)):
            xbins[j] = xaxis.GetBinLowEdge(j+1+skipFirstX)
        return dataArray,errorArray,xbins,ybins
       
def GettData(playlist):
    migrationPath = "/minerva/data/users/ajball/nu_e/kin_dist_mc"+str(playlist)+"_collab1_fspline.root"
    IUEPath = "/minerva/data/users/ajball/nu_e/kin_dist_data"+str(playlist)+"_subtracted_collab1_fspline.root"
    efficiencyPath = "/minerva/app/users/ajball/ncdiff/CC-NuE-XSec/efficiency/"+str(playlist)+"_calculatedEfficiency.root"
    migrationFile = ROOT.TFile.Open(migrationPath,"READ")
    IUEFile = ROOT.TFile.Open(IUEPath,"READ")
    efficiencyFile = ROOT.TFile.Open(efficiencyPath,"READ")
    migrationHist = migrationFile.inline_vs_t.Clone()
    IUEHist = IUEFile.InlineUpstream_Energy.Clone()
    t_efficiencyHist = efficiencyFile.t_efficiency.Clone()
    migrationMatrix = np.zeros((migrationHist.GetNbinsY(),migrationHist.GetNbinsX()-1))
    IUE_data = np.zeros(IUEHist.GetNbinsX()-1)
    IUE_bins = np.zeros(len(IUE_data)+1)
    t_efficiency = np.zeros(t_efficiencyHist.GetNbinsX())
    t_bins = np.zeros(len(t_efficiency)+1)
    for i in range(len(migrationMatrix)):
        t_efficiency[i] = t_efficiencyHist.GetBinContent(i+1)
        for j in range(len(migrationMatrix[0])):
            if i == 0:
                IUE_data[j] = IUEHist.GetBinContent(j+2)
            migrationMatrix[i,j] = migrationHist.GetBinContent(j+2,i+1)
        migrationMatrix[i] = migrationMatrix[i]*(1/np.sum(migrationMatrix[i]))
    for i in range(len(t_bins)):
        t_bins[i] = t_efficiencyHist.GetBinLowEdge(i+1)
    for j in range(len(IUE_bins)):
        IUE_bins[j] = IUEHist.GetBinLowEdge(j+2)
    return migrationMatrix,IUE_data,IUE_bins,t_efficiency,t_bins
playlist = "Alex_tcut_nx"
blip = GettData(playlist)[0]
efficiencyPath = "/minerva/app/users/ajball/ncdiff/CC-NuE-XSec/efficiency/"+str(playlist)+"_calculatedEfficiency.root"
subtractedDataPath = "/minerva/data/users/ajball/nu_e/kin_dist_data"+str(playlist)+"_subtracted_collab1_fspline.root"
migrationPath = "/minerva/data/users/ajball/nu_e/kin_dist_mc"+str(playlist)+"_collab1_fspline.root"

t_efficiency = Variable1D("t_efficiency",efficiencyPath,isEfficiency = True)
IUE = Variable1D("InlineUpstream_Energy",subtractedDataPath,skipFirst = True)
t_migration = MigrationMatrix("inline_vs_t",migrationPath,skipFirstX = True)

pion_efficiency = Variable1D("pion_efficiency",efficiencyPath,isEfficiency = True)
Eel = Variable1D("Eel",subtractedDataPath)
pion_migration = MigrationMatrix("ee_vs_pionE",migrationPath)

print(t_efficiency.values,t_efficiency.error,t_efficiency.bins)
print(IUE.values,IUE.error,IUE.bins)
print(t_migration.values,t_migration.error,t_migration.xbins,t_migration.ybins)
print(blip)
