"""
Module/Package Description

Author: Stuart J Goldie, Nora Campoli
Date: 25-02-2026
Version: 2.0
License: MIT
"""

__author__ = "Stuart J Goldie, Nora Campoli"
__version__ = "2.0"
__license__ = "MIT"

import os
import numpy as np
import matplotlib.pyplot as plt
from tkinter import Tk
from scipy import stats
from lmfit.models import PolynomialModel, LorentzianModel
from tkinter.filedialog import askdirectory
#from FunctionBank import multiple_output
import csv

#%% Functions
# Set-up functions and constants in relevant units
def fitting_output_search(output, parameter):
    """"search through parameter string produced by lmfit fitting module and return values in a dictionary"""
    string = str(output.params[parameter])
    index_initial = string.find('value=')+6
    index_mid = string.find('+/-')
    index_final = string.find('bounds')-2
    value = float(string[index_initial:index_mid-1])
    error = float(string[index_mid+4:index_final])
    return {parameter: value, parameter + '_error': error}

def peak_output(destination, source, peak_name):
    """Stores fitted outputs for a given peak prefix (e.g., 'si', 'e2g', 'a1g')."""
    params = ['amplitude', 'center', 'fwhm', 'height']
    combined = {}
    for param in params:
        combined.update(
            fitting_output_search(source, f'{peak_name}_{param}')
        )
    destination.update(combined)

def multiple_output(destination, source):
    """stores fitted outputs for each peak in turn"""
    peak_output(destination, source, 'si')
    peak_output(destination, source, 'e2g')
    peak_output(destination, source, 'a1g')

def estimate_noise_sliding_poly(row_data, window_size=21, poly_order=2):
    """Estimate noise by fitting a degree-2 polynomial to a sliding window and
    computing the mean absolute difference between the data and the fit."""
    n = len(row_data)
    half_window = window_size // 2
    indices = np.arange(n)
    window_errors = []
    for i in range(n):
        window_start = max(0, i - half_window)
        window_end = min(n, i + half_window + 1)
        win_idx = indices[window_start:window_end]
        win_data = row_data[window_start:window_end]
        coeffs = np.polyfit(win_idx, win_data, poly_order)
        fitted = np.polyval(coeffs, win_idx)
        window_errors.append(np.mean(np.abs(win_data - fitted)))
    return np.mean(window_errors)
   

h = 4.1356692*10**-15 #eV s
c = 3.0*10**8 #m/s
lam = 532 #nm

material = "MoS2"

peak_configs = {
        "MoS2": {"A1g_upper": 415, "A1g_lower": 395, "E2g_lower": 375},
        "WS2": {"A1g_upper": 422, "A1g_lower": 400, "E2g_lower": 375},
    }

#%% Open and parse data
# Open and parse files

Tk().withdraw()
directory = askdirectory()
filenames = [f for f in os.listdir(directory) if f.endswith('.txt')]
map_files = [f for f in filenames if '_map' in f]
if len(map_files) > 1:
    print("Warning: Multiple '_map' files found. Using the largest.")
    map_file = max(map_files, key=lambda f: os.path.getsize(os.path.join(directory, f)))
elif len(map_files) == 1:
    map_file = map_files[0]
else:
    map_file = None

#Parsing photoluminescence data sets
PL_files = [f for f in filenames if '_PL' in f]
PL_wavenumbers = []
PL_counts = []
for filename in PL_files:
    current_PL_data = np.genfromtxt(os.path.join(directory, filename), dtype=None, delimiter='\t')
    PL_wavenumbers.append(current_PL_data[:,0])
    PL_counts.append(current_PL_data[:,1])

if all(np.array_equal(arr, PL_wavenumbers[0]) for arr in PL_wavenumbers):
    #PL_data = PL_wavenumbers[0]
    PL_data = np.vstack((PL_wavenumbers[0],PL_counts))
else:
    print("Error: PL wavenumber range is not identical across all files.")

#Parsing map data (contained function stored in FunctionBank for possible later deployment)
if map_file != None:
    datalist = np.genfromtxt(os.path.join(directory, map_file), dtype=None, delimiter='\t')
    
    point_index = []
    for r in range(1,len(datalist)-1):
        if float(datalist[r+1][0]) - float(datalist[r][0]) != 0:
            point_index.append(r)
    point_index.append(len(datalist))
    map_data = np.zeros([len(point_index)+1,int(point_index[2]-point_index[1])],dtype=float)
    for r in range(1,int(point_index[2]-point_index[1])+1):
        map_data[0,r-1]=float(datalist[r-1][2])
    
    row = 1
    column = 0
    for r in range(0,len(datalist)):
        map_data[row,column] = float(datalist[r][3])
        if r in point_index:
            row += 1
            column = 0
        else:
            column +=1

else:
    print("No map file found. This error handling not included yet")
    #if other files are present, combine them into a psudo map and analyse accordingly
    #if no suitable files are present, report the false folder and re-open
    #this quality of life change for later

#%% Map Fitting
# Fit and analyse statistics of Raman peak region
raman_modes_map_data = map_data[:, map_data[0, :] < 900]
bkgd_map_data = raman_modes_map_data[
    :,
    (raman_modes_map_data[0, :] < peak_configs[material]["E2g_lower"])
    | (raman_modes_map_data[0, :] > peak_configs[material]["A1g_upper"])
]
bkgd_map_data = bkgd_map_data[
    :,
    (bkgd_map_data[0, :] < 505)
    | (bkgd_map_data[0, :] > 535)
]
si_map_data = raman_modes_map_data[
    :,
    (raman_modes_map_data[0, :] > 500)
    | (raman_modes_map_data[0, :] < 540)
]

map_x = raman_modes_map_data[0,:]
results = []
print_flag = 'none'
null_fit_list = []
null_fit_count = 0

error_array = np.array([
    estimate_noise_sliding_poly(row_data)
    for row_data in bkgd_map_data[1:]
])

#Specify peak models needed during fitting
si_peak = LorentzianModel(nan_policy='propagate',prefix='si_')
e2g_peak = LorentzianModel(nan_policy='propagate',prefix='e2g_')
a1g_peak = LorentzianModel(nan_policy='propagate',prefix='a1g_')
poly2 = PolynomialModel(1,nan_policy='propagate', prefix='Background_')

model_p = poly2
model_si = poly2 + si_peak
model_tmd = poly2 + si_peak + e2g_peak + a1g_peak

a1g_up = peak_configs[material]["A1g_upper"]
a1g_up_index = np.argmin(np.abs(map_x - a1g_up))
a1g_low = peak_configs[material]["A1g_lower"]
a1g_low_index = np.argmin(np.abs(map_x - a1g_low))
e2g_low = peak_configs[material]["E2g_lower"]
e2g_low_index = np.argmin(np.abs(map_x - e2g_low))

for spectrum_index, spectral_data in enumerate(raman_modes_map_data[1:,:]):
    pars_poly2 = poly2.guess(spectral_data, x=map_x)
    out_poly = model_p.fit(spectral_data, pars_poly2, x=map_x, weights=(1/error_array[spectrum_index]))

    pars_si = poly2.guess(bkgd_map_data[spectrum_index,:], x=bkgd_map_data[0,:])
    pars_si.update(si_peak.guess(si_map_data[spectrum_index,:], x=si_map_data[0,:]))
    pars_si['si_center'].set(min=515, max=525)
    pars_si['si_amplitude'].set(min=1.e-09)
    pars_si['si_sigma'].set(max=10)
    out_si = model_si.fit(spectral_data, pars_si, x=map_x, weights=(1/error_array[spectrum_index]))

    if ((out_si.redchi)/(out_poly.redchi)) < 0.975:
        pars_tmd = poly2.guess(bkgd_map_data[spectrum_index,:], x=bkgd_map_data[0,:])
        pars_tmd['Background_c0'].set(np.mean(bkgd_map_data[spectrum_index,:]), min=0)
        pars_tmd.update(si_peak.guess(si_map_data[spectrum_index,:], x=si_map_data[0,:]))
        pars_tmd['si_center'].set(min=505, max=535)
        pars_tmd['si_amplitude'].set(min=1.e-09)
        pars_tmd['si_sigma'].set(2, max=10)
        pars_tmd.update(e2g_peak.guess(spectral_data[a1g_low_index:e2g_low_index], x=map_x[a1g_low_index:e2g_low_index]))
        pars_tmd['e2g_center'].set(min=e2g_low, max=a1g_low)
        pars_tmd['e2g_amplitude'].set(min=1.e-09)
        pars_tmd['e2g_sigma'].set(max=10)
        pars_tmd.update(a1g_peak.guess(spectral_data[a1g_up_index:a1g_low_index], x=map_x[a1g_up_index:a1g_low_index]))
        pars_tmd['a1g_center'].set(min=a1g_low, max=a1g_up)
        pars_tmd['a1g_amplitude'].set(min=1.e-09)
        pars_tmd['a1g_sigma'].set(max=10)
        out_tmd = model_tmd.fit(spectral_data, pars_tmd, x=map_x, weights=(1/error_array[spectrum_index]))
        
        if ((out_tmd.redchi)/(out_si.redchi)) < 0.975:
            full_result = {}
            if 'report' in print_flag:
                print(f'Spectrum: {spectrum_index}')
                print(out_tmd.fit_report(show_correl=False))
            
            if 'plot' in print_flag:
                fig, ax = plt.subplots(figsize=(10, 6))
                ax.plot(map_x, spectral_data, 'o', label='Data', markersize=4)
                ax.plot(map_x, out_tmd.best_fit, '-', label='Fit', linewidth=2)
                ax.set_xlabel('Wavenumber (cm⁻¹)')
                ax.set_ylabel('Intensity (counts)')
                ax.set_title(f'Spectrum {spectrum_index + 1}')
                ax.legend()
                ax.grid(True, alpha=0.3)
                plt.tight_layout()
                plt.show()
            
            try:
                multiple_output(full_result, out_tmd)
                full_result.update({'Red Chi2':(out_tmd.redchi)})
                full_result.update({'a1g-e2g_spacing': (full_result['a1g_center'] - full_result['e2g_center'])})
                full_result.update({'a1g-e2g_spacing_error': (np.sqrt(full_result['a1g_center_error']**2 + full_result['e2g_center_error']**2))})
                results.append(full_result)
            
            except ValueError:
                null_fit_count += 1
                null_fit_list.append(spectrum_index + 1)
        else:
            null_fit_count += 1
            null_fit_list.append(spectrum_index + 1)

headings = [
    "si_amplitude",
    "si_amplitude_error",
    "si_center",
    "si_center_error",
    "si_fwhm",
    "si_fwhm_error",
    "si_height",
    "si_height_error",
    "e2g_amplitude",
    "e2g_amplitude_error",
    "e2g_center",
    "e2g_center_error",
    "e2g_fwhm",
    "e2g_fwhm_error",
    "e2g_height",
    "e2g_height_error",
    "a1g_amplitude",
    "a1g_amplitude_error",
    "a1g_center",
    "a1g_center_error",
    "a1g_fwhm",
    "a1g_fwhm_error",
    "a1g_height",
    "a1g_height_error",
    "a1g-e2g_spacing",
    "a1g-e2g_spacing_error",
    "Red Chi2"
]

result_array = np.array([[row[h] for h in headings] for row in results])
#for development, do not keep saving repeat data files
#output_filename = os.path.join(directory, 'map_fitting.csv')
#np.savetxt(output_filename, result_array, delimiter=',', header=','.join(headings), comments='')

'''
with open(output_filename, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(headings)
    writer.writerows(result_array)'''

print('Spectra processed: ' + str(spectrum_index + 1))

if null_fit_count == 0:
    print('Percentage bad fit: 0%')
else:
    print('Bad fit: ' + str(null_fit_count))
    print(null_fit_list)

print('_______________________________________________')

#Shapiro-Wilk normality test. 
#If p-value < 0.05, there is sufficient evidence to reject the null 
# hypothesis that the data comes from a normal distribution.
shapiro_test = stats.shapiro(result_array[:,24])
if shapiro_test.pvalue < 0.05:
    print('**Statistical tests suggest this data set is not normally distributed '
    'and caution should be employed if using simple summary statistics.**')
    print('____')
print(f'Mean Si peak position = {np.mean(result_array[:,2]):.2f} cm⁻¹')
print(f'Range Si peak position = {np.max(result_array[:,2]) - np.min(result_array[:,2]):.2f} cm⁻¹')
print(f'Mean TMD peak separation = {np.mean(result_array[:,24]):.2f} cm⁻¹')
print(f'SD TMD peak separation = {np.std(result_array[:,24]):.2f} cm⁻¹')

#%% Photoluminescence
# Converting and analysing PL spectra    