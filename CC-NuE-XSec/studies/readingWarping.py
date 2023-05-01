import ROOT,PlotUtils
import heapq

file = "transWrapO_LowQ2Pi0Eavail_Lepton_Pt.root"
i_uni = 2
iter_num = 1

def translate_bin(i,nx):
    y = i/nx
    x = i%nx
    return (x+1,y+1)


f = ROOT.TFile.Open(file)
chiMap =f.Get("Chi2Maps").Get("Chi2Map_Stat_{}_Iter_{}".format(i_uni,iter_num))
chi2s = []
chi2 =0
Cov = f.Get("Cov_Matrix").Get("CovMatrix_Stat_{}_Iter_{}".format(i_uni,iter_num))
print((Cov.GetNcols()))
data_truth = f.Get("Input_Hists").Get("h_data_truth")
unfolded = f.Get("Unfolded_Data").Get("Stat_{}_Iter_{}".format(i_uni,iter_num))

for i in range(data_truth.GetSize()):
    Cov[i][i]+=data_truth.GetBinError(i)**2




c =ROOT.TCanvas("c1","c1",1920,1080)
chiMap.Draw("COLZ")
c.Print("chi2.png")

Cov.Draw("COLZ")
c.Print("Cov.png")


for i in range(chiMap.GetNcols()):
    for j in range(chiMap.GetNcols()):
        chi2+=chiMap[i][j]
        heapq.heappush(chi2s,(-abs(chiMap[i][j]),i,j))


print((data_truth,unfolded))
NX = data_truth.GetNbinsX()

print (chi2)
for i in range(chiMap.GetNcols()):
    chi2_neg,i,j = heapq.heappop(chi2s)
    xi,yi = translate_bin(i,NX)
    xj,yj = translate_bin(j,NX)
    print((i,j,chiMap[i][j],Cov[yi*(NX+2)+xi][yj*(NX+2)+xj]))
    if i==j:
        print((unfolded.GetBinContent(xi,yi),data_truth.GetBinContent(xi,yi)))
