# Alex Ball
# Matrix Fitting Procedure

import numpy as np
import math
import ROOT
import PlotUtils

# Binning and Matrix Algebra and Whatnot ------------------------------------------------------------------                        
def simpson_integration(x0,xf,N,integrand): # Defines a generalized Simpson's rule algorithm that takes an integrand function
    xvals = np.linspace(x0,xf,N+1,endpoint=True) # Finds evenly spaced x values
    yvals = []
    for xval in xvals: # Creates a list of y values for each x value
        yvals.append(integrand(xval))
    h = (xf-x0)/N # Finds the distance between each x value
    integral = h/3*(yvals[0]+yvals[-1]+4*sum(yvals[1:N:2])+2*sum(yvals[2:N:2])) # Finds the integral with Simpson's rule
    return integral

def y_bin_maker(r_vec,y_bins): # r_vec is [N,b]
    n = len(y_bins)-1 # Size of the 2D Histogram Grid
    dist = lambda x: r_vec[0]*np.exp(-x/r_vec[1])/r_vec[1] # The t distribution function
    y_bin_vals = np.zeros(n)
    for i in range(n):
        y_bin_vals[i] = simpson_integration(y_bins[i],y_bins[i+1],100,dist)
    return y_bin_vals

def x_bin_maker(r_vec,y_bins,matrix): # r_vec is [N,b]
    y_bin_vals = y_bin_maker(r_vec,y_bins)
    x_bin_vals = np.dot(matrix,y_bin_vals)
    return x_bin_vals

def recover_y(x_bin_vals,matrix):
    invmatrix = np.linalg.inv(matrix)
    return np.dot(invmatrix,x_bin_vals)

def variance(r,data):
    N,b = r
    actual,y_bins,matrix = data
    fit = x_bin_maker(np.array([N,b]),y_bins,matrix)
    total_variance = 0
    for i in range(len(actual)):
        total_variance += (actual[i]-fit[i])**2/fit[i]
    return total_variance

def sens_test(true_x_dist,y_bins,matrix_elements):
    # Fitting ------------------------------------------------------------------
    N_list = []
    b_list = []
    experiments = 100
    for j in range(experiments):
        x_dist = np.zeros(len(true_x_dist))
        for i in range(len(true_x_dist)): # Adds random variation using the poisson distribution
            x_dist[i] = np.random.poisson(true_x_dist[i])
        r = np.array([1000,0.1]) # Starting guess
        results = minimize(variance,args=[x_dist,y_bins,matrix_elements],x0=r,bounds=((0,10000),(0,5)),method="Nelder-Mead")
        r = results.x

        print("Iteration",j,"The fit values were N =",r[0],"b =",r[1])
        N_list.append(r[0])
        b_list.append(r[1])

    N_avg = np.average(N_list)
    N_uncert = np.std(N_list)
    b_avg = np.average(b_list)
    b_uncert = np.std(b_list)
    rho = np.sum(np.multiply(np.array(N_list)-N_avg,np.array(b_list)-b_avg))/(N_uncert*b_uncert)/experiments
    print("N =",pres_format(N_avg,N_uncert),"b =",pres_format(b_avg,b_uncert))
    return N_avg,N_uncert,b_avg,b_uncert,rho

def dN_dt(N,b,t):
    return N/b*np.exp(-t/b)

def fit_uncert(N,b,dN,db,rho,ti,tf): # Uncertainty in dN_dt
    dF_db = N*(ti*np.exp(-ti/b)-tf*np.exp(-tf/b))/(b**2*(tf-ti))
    dF_dN = (np.exp(-ti/b)-np.exp(-tf/b))/(tf-ti)
    return np.sqrt((dF_dN*dN)**2+(dF_db*db)**2+2*rho*(dF_dN*dN)*(dF_db*db))

def pres_format(value,uncertainty): # Only works for floats
    place = int(math.floor(math.log10(uncertainty))) - 1
    final_value = round(value,-place)
    final_uncertainty = round(uncertainty,-place)
    if final_value == int(final_value) and final_uncertainty == int(final_uncertainty): 
        final_value = int(final_value)
        final_uncertainty = int(final_uncertainty)
    return str(final_value) +" ± "+ str(final_uncertainty)


# Extracting the data from the text files into rebinned histograms -------------------------------------------
matrix_elements,bin_coords,bin_error = readTXT("Histogram Data/inline_vs_t.txt",2)
t_bins = np.array([0,0.04,0.07,0.1,0.15,0.5])
U_bins = np.array([10,15,25,35,50,150])
new_matrix_elements,new_bin_coords = rebin2D(matrix_elements,bin_coords,U_bins,t_bins)
migration = np.zeros((len(U_bins)-1,len(U_bins)-1))
for i in range(len(new_matrix_elements)):
    migration[i] = new_matrix_elements[i]*1/np.sum(new_matrix_elements[i])
print(migration)
del new_matrix_elements; del matrix_elements

x_bin_size = []
x_vals_fit = []
for i in range(len(t_bins)-1):
    x_vals_fit.append(np.average([t_bins[i+1],t_bins[i]]))
    x_bin_size.append(x_vals_fit[i]-t_bins[i])
print("X binsize",x_bin_size)
U_vector,U_locs,U_error = readTXT("Histogram Data/InlineUpstream_Energy.txt",1)
nonnormalized_U_dist,U_locs = rebin1D(U_vector,U_locs,U_bins)
true_U_dist = nonnormalized_U_dist*1000/np.sum(nonnormalized_U_dist)
print(true_U_dist)

# Matrix Psedoexperiment stuff ---------------------------------------------------
t_bin_list = []
t_uncert_list = []
t_vectors = []
for i in range(100000):
    if i % 1000 == 0 : print("Matrix Pseudoexperiment #",str(i))
    U_dist = np.zeros(len(true_U_dist))
    for j in range(len(true_U_dist)): # Adds random variation using the poisson distribution
        U_dist[j] = np.random.poisson(true_U_dist[j])
    t_vectors.append(recover_y(U_dist,migration))
t_vectors = np.array(t_vectors)
avg_t_bins, uncert_t_bins = np.zeros(len(true_U_dist)), np.zeros(len(true_U_dist))
recovered_t_bins = []
for i in range(len(t_vectors[0])):
    avg_t_bins[i] = np.average(t_vectors[:,i])
    uncert_t_bins[i] = np.std(t_vectors[:,i])
for i in range(len(avg_t_bins)): # Converts to dN/d|t|
    avg_t_bins[i] /= (t_bins[i+1]-t_bins[i])
    uncert_t_bins[i] /= (t_bins[i+1]-t_bins[i])
    recovered_t_bins.append(pres_format(avg_t_bins[i],uncert_t_bins[i]))
print(recovered_t_bins)

# Fitting to Randomly Varied Truth bins ---------------------------------------
fit_data = sens_test(true_U_dist,t_bins,migration)
bin_vals_fit = y_bin_maker(np.array([fit_data[0],fit_data[2]]),t_bins)
y_vals_fit = np.zeros(len(bin_vals_fit))
y_err_fit = np.zeros(len(bin_vals_fit))
t_overall_fit = []
for j in range(len(bin_vals_fit)):
    y_vals_fit[j] = bin_vals_fit[j]/(t_bins[j+1]-t_bins[j])
    y_err_fit[j] = fit_uncert(fit_data[0],fit_data[2],fit_data[1],fit_data[3],fit_data[4],t_bins[j],t_bins[j+1])
    t_overall_fit.append(pres_format(y_vals_fit[j],y_err_fit[j]))

