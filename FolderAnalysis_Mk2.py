"""
Module/Package Description

Author: Stuart J Goldie, Nora Campoli
Date: 25-02-2026
Version: 2.1
License: MIT
Comments: V2.1 assumes MoS2 is the target material. 
Other materials with different peak parameters will 
likely experience instability with fitting and analysis. 
"""

__author__ = "Stuart J Goldie, Nora Campoli"
__version__ = "2.1"
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

def multiple_output(destination, source, silicon=True):
    """stores fitted outputs for each peak in turn"""
    if silicon == True:
        peak_output(destination, source, 'si')
    peak_output(destination, source, 'e2g')
    peak_output(destination, source, 'a1g')

def multiple_output_PL(destination, source):
    """stores fitted outputs for each PL peak in turn"""
    peak_output(destination, source, 'xb')
    peak_output(destination, source, 'xa')
    xtneg_result = {}
    peak_output(xtneg_result, source, 'xtneg')
    for key in list(xtneg_result.keys()):
        destination[key.replace('xtneg', 'xt-', 1)] = xtneg_result[key]

def hex_to_RGB(hex_str):
    """ #FFFFFF -> [255,255,255]"""
    #Pass 16 to the integer function for change of base
    return [int(hex_str[i:i+2], 16) for i in range(1,6,2)]
 
def get_color_gradient(c1, c2, n):
    """
    Given two hex colors, returns a color gradient
    with n colors.
    """
    if n > 1:
        c1_rgb = np.array(hex_to_RGB(c1))/255
        c2_rgb = np.array(hex_to_RGB(c2))/255
        mix_pcts = [x/(n-1) for x in range(n)]
        rgb_colors = [((1-mix)*c1_rgb + (mix*c2_rgb)) for mix in mix_pcts]
        return ["#" + "".join([format(int(round(val*255)), "02x") for val in item]) for item in rgb_colors]
    
    elif n==1:
        return['b']

def find_raman_peaks(data, start_index, material="MoS2"):
    """
    Find Raman peaks for specified material.
    
    Parameters:
    - data: numpy array with wavenumber (col 0) and intensity (col 1)
    - start_index: offset index for peak positions
    - material: material type (default "MoS2")
    
    Returns:
    - [A1g_peak_index, E2g_peak_index]
    """
    
    peak_configs = {
        "MoS2": {"A1g_upper": 415, "A1g_lower": 395, "E2g_lower": 375},
        "WS2": {"A1g_upper": 422, "A1g_lower": 400, "E2g_lower": 375},
    }
    
    config = peak_configs.get(material, peak_configs["MoS2"])
    
    A1g_index_1 = np.where(data[:,0] < config["A1g_upper"])[0][0]
    A1g_index_2 = A1g_index_1 + np.where(data[A1g_index_1:,0] < config["A1g_lower"])[0][0]
    A1g_peak_index = A1g_index_1 + np.argmax(data[A1g_index_1:A1g_index_2, 1])

    E2g_index_1 = A1g_index_2
    E2g_index_2 = E2g_index_1 + np.where(data[E2g_index_1:,0] < config["E2g_lower"])[0][0]
    E2g_peak_index = E2g_index_1 + np.argmax(data[E2g_index_1:E2g_index_2, 1])

    return [start_index + A1g_peak_index, start_index + E2g_peak_index]

import math

def format_with_error(value, error):
    """
    Format a value and its uncertainty.

    Parameters
    ----------
    value : float
        Central (measured or fitted) value.
    error : float
        Associated uncertainty (absolute, not percent).

    Returns
    -------
    str
        Formatted string of the form "value ± uncertainty", with matched decimal precision.
    """
    if error == 0:
        return f"{value} ± 0"

    abs_err = abs(error)

    exponent = math.floor(math.log10(abs_err))
    first_digit = int(abs_err / 10**exponent)
    sig_figs = 2 if first_digit == 1 else 1

    rounded_error = round(abs_err, -(exponent - (sig_figs - 1)))
    decimal_places = max(0, -(exponent - (sig_figs - 1)))
    rounded_value = round(value, decimal_places)

    return f"{rounded_value:.{decimal_places}f} ± {rounded_error:.{decimal_places}f}"


h = 4.1356692*10**-15 #eV s
c = 3.0*10**8 #m/s
lam = 532 #nm

material = "MoS2"

peak_configs = {
        "MoS2": {"A1g_upper": 415, "A1g_lower": 395, "E2g_lower": 375},
        "WS2": {"A1g_upper": 422, "A1g_lower": 400, "E2g_lower": 375},
    }

#%% Map Analysis
# Analysing Raman Map Data
Tk().withdraw()
directory = askdirectory()
print(f'Opened folder: {directory}')

filenames = [f for f in os.listdir(directory) if f.endswith('.txt')]
map_files = [f for f in filenames if '_map' in f.lower()]
if len(map_files) > 1:
    print("Warning: Multiple '_map' files found. Using the largest.")
    map_file = max(map_files, key=lambda f: os.path.getsize(os.path.join(directory, f)))
elif len(map_files) == 1:
    map_file = map_files[0]

try:
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

    # Fit and analyse statistics of Raman peak region
    raman_modes_map_data = map_data[:, (map_data[0, :] > 250) & (map_data[0, :] < 600)]
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
    #map_flag = input('Additional Raman data outputs? [plot, report]')
    map_flag = 'none'
    null_fit_list = []
    null_fit_count = 0
    plot_xlims = (360,430)
    
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
    #Complete fitting sequence. By validating an improvement in fit with each model, low quality data is excluded
    for spectrum_index, spectral_data in enumerate(raman_modes_map_data[1:,:]):
        pars_poly2 = poly2.guess(spectral_data, x=map_x)
        out_poly = model_p.fit(spectral_data, pars_poly2, x=map_x)
        pars_si = poly2.guess(bkgd_map_data[spectrum_index,:], x=bkgd_map_data[0,:])
        pars_si.update(si_peak.guess(si_map_data[spectrum_index,:], x=si_map_data[0,:]))
        pars_si['si_center'].set(min=515, max=525)
        pars_si['si_amplitude'].set(min=1.e-09)
        pars_si['si_sigma'].set(max=10)
        out_si = model_si.fit(spectral_data, pars_si, x=map_x)
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
            out_tmd = model_tmd.fit(spectral_data, pars_tmd, x=map_x)
            
            if ((out_tmd.redchi)/(out_si.redchi)) < 0.975:
                full_result = {}
                if 'report' in map_flag:
                    print(f'Spectrum: {spectrum_index}')
                    print(out_tmd.fit_report(show_correl=False))
                
                if 'plot' in map_flag:
                    fig, ax = plt.subplots(figsize=(10, 6))
                    ax.plot(map_x, spectral_data, 'o', label='Data', markersize=4)
                    ax.plot(map_x, out_tmd.best_fit, '-', label='Fit', linewidth=2)
                    ax.set_xlabel('Wavenumber / cm⁻¹')
                    ax.set_ylabel('Intensity / counts')
                    ax.set_title(f'Spectrum {spectrum_index + 1}')
                    ax.legend()
                    ax.grid(True, alpha=0.3)
                    plt.tight_layout()
                    plt.show()
                
                try:
                    multiple_output(full_result, out_tmd)
                    full_result.update({'a1g-e2g_spacing': (full_result['a1g_center'] - full_result['e2g_center'])})
                    full_result.update({'a1g-e2g_spacing_error': (np.sqrt(full_result['a1g_center_error']**2 + full_result['e2g_center_error']**2))})
                    full_result.update({'spectrum_index': spectrum_index})
                    results.append(full_result)
                
                except ValueError:
                    null_fit_count += 1
                    null_fit_list.append(spectrum_index + 1)
            else:
                null_fit_count += 1
                null_fit_list.append(spectrum_index + 1)
        else:
            pars_tmd = poly2.guess(bkgd_map_data[spectrum_index,:], x=bkgd_map_data[0,:])
            pars_tmd['Background_c0'].set(np.mean(bkgd_map_data[spectrum_index,:]), min=0)
            pars_tmd.update(e2g_peak.guess(spectral_data[a1g_low_index:e2g_low_index], x=map_x[a1g_low_index:e2g_low_index]))
            pars_tmd['e2g_center'].set(min=e2g_low, max=a1g_low)
            pars_tmd['e2g_amplitude'].set(min=1.e-09)
            pars_tmd['e2g_sigma'].set(max=10)
            pars_tmd.update(a1g_peak.guess(spectral_data[a1g_up_index:a1g_low_index], x=map_x[a1g_up_index:a1g_low_index]))
            pars_tmd['a1g_center'].set(min=a1g_low, max=a1g_up)
            pars_tmd['a1g_amplitude'].set(min=1.e-09)
            pars_tmd['a1g_sigma'].set(max=10)
            out_tmd = model_tmd.fit(spectral_data, pars_tmd, x=map_x)

            if ((out_tmd.redchi)/(out_poly.redchi)) < 0.975:
                full_result = {}
                if 'report' in map_flag:
                    print(f'Spectrum: {spectrum_index}')
                    print(out_tmd.fit_report(show_correl=False))
                
                if 'plot' in map_flag:
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
                    multiple_output(full_result, out_tmd, silicon=False)
                    full_result.update({'a1g-e2g_spacing': (full_result['a1g_center'] - full_result['e2g_center'])})
                    full_result.update({'a1g-e2g_spacing_error': (np.sqrt(full_result['a1g_center_error']**2 + full_result['e2g_center_error']**2))})
                    full_result.update({'spectrum_index': spectrum_index})
                    full_result.update({"si_amplitude": 0.0, "si_amplitude_error": 0.0, "si_center": 0.0, "si_center_error": 0.0, 
                                        "si_fwhm": 0.0, "si_fwhm_error": 0.0, "si_height": 0.0, "si_height_error": 0.0})
                    results.append(full_result)
                
                except ValueError:
                    null_fit_count += 1
                    null_fit_list.append(spectrum_index + 1)
            else:
                null_fit_count += 1
                null_fit_list.append(spectrum_index + 1)
    
    headings = [
        "spectrum_index",
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
        "a1g-e2g_spacing_error"
    ]
    result_array = np.array([[row[h] for h in headings] for row in results])
    np.savetxt(os.path.join(directory, 'map_fitting.csv'), result_array, 
               delimiter=',', header=','.join(headings), comments='')
    
    print('Spectra processed: ' + str(spectrum_index + 1))
    if null_fit_count != 0:
        print('Bad fit: ' + str(null_fit_count))
        print(null_fit_list)
    print('_______________________________________________')
    #Shapiro-Wilk normality test. 
    #If p-value < 0.05, there is sufficient evidence to reject the null 
    # hypothesis that the data comes from a normal distribution.
    shapiro_test = stats.shapiro(result_array[:,25])
         
    silicon_position = np.mean(result_array[:,3])
    if silicon_position < 519:
        print('**Silicon peak is out of expected range**')
        print('Spectrometer may be incorrectly calibrated')
    elif silicon_position > 521:
        print('**Silicon peak is out of expected range**')
        print('Spectrometer may be incorrectly calibrated')
    
    print(f'Mean Si peak position = {silicon_position:.2f} cm⁻¹')
    print(f'Range Si peak position = {np.max(result_array[:,3]) - np.min(result_array[:,3]):.2f} cm⁻¹')
    print('\n' \
    '____TMD Stats____')
    
    if shapiro_test.pvalue < 0.05:
        print('Statistical tests suggest this data set is not normally distributed.')
        print('Caution if using mean and standard deviation.')
        print('Representative spectrum determined by modal value, not mean.')
        print(f'Shapiro-Wilk normality test result: p = {shapiro_test.pvalue:.3f}')
        print('______________________')
    
    mean_spacing = np.mean(result_array[:,25])
    histogram = np.histogram(result_array[:,25], bins = 15)
    modal_index = np.argmax(np.histogram(result_array[:,25], bins = 15)[0])
    modal_spacing = (np.histogram(result_array[:,25], bins = 15)[1][modal_index] + 
                    np.histogram(result_array[:,25], bins = 15)[1][modal_index+1])/2
    if shapiro_test.pvalue >= 0.05:
        representative_index = int(result_array[np.argmin(np.abs(result_array[:,25] - mean_spacing)),0])
    elif shapiro_test.pvalue < 0.05:
        representative_index = int(result_array[np.argmin(np.abs(result_array[:,25] - modal_spacing)),0])
    print(f'Mean TMD peak separation = {np.mean(result_array[:,25]):.2f} cm⁻¹')
    print(f'SD TMD peak separation = {np.std(result_array[:,25]):.2f} cm⁻¹')
    print(f'Modal TMD peak separation = {modal_spacing:.2f} cm⁻¹')
    print(f'Range TMD peak separation = {np.max(result_array[:,25]) - np.min(result_array[:,25]):.2f} cm⁻¹')
    print(f'Most representative spectrum = {representative_index + 1}')
    print('_______________________________________________')
    
    representative_spectrum_data = np.column_stack((raman_modes_map_data[0, :], raman_modes_map_data[representative_index + 1, :]))
    np.savetxt(os.path.join(directory, 'representative_spectrum.csv'), 
               representative_spectrum_data, delimiter=',', header='Wavenumber,Intensity', comments='')
    
    plt.figure(figsize=(12/2.54, 9/2.54))
    for i in range(1, raman_modes_map_data.shape[0]):
        y_data = raman_modes_map_data[i, :]
        x_data = raman_modes_map_data[0, :]
        mask = (x_data >= plot_xlims[0]) & (x_data <= plot_xlims[1])
        y_max = np.max(y_data[mask])
        normalised_y = y_data / y_max
        plt.plot(x_data, normalised_y, alpha=0.3, color='royalblue', linewidth=0.5)
    representative_y = raman_modes_map_data[representative_index + 1, :]
    representative_x = raman_modes_map_data[0, :]
    mask = (representative_x >= plot_xlims[0]) & (representative_x <= plot_xlims[1])
    y_max = np.max(representative_y[mask])
    normalised_representative_y = representative_y / y_max
    plt.plot(representative_x, normalised_representative_y, alpha=1.0, linewidth=1, color='blue')
    plt.xlabel('Wavenumber / $cm^{-1}$', fontname='Arial', fontsize=12)
    plt.xticks(fontname='Arial', fontsize=10)
    plt.yticks(fontname='Arial', fontsize=10)
    plt.xlim(plot_xlims)
    plt.ylabel('Normalized Intensity', fontname='Arial', fontsize=12)
    plt.ylim((0,1.05))
    plt.title('Map Spectra')
    plt.tight_layout()
    plt.savefig(os.path.join(directory, "map_spectra.png"), dpi=600)

except NameError:
    print('No map files found in selected directory. Ensure the map file is suffixed \'_map\'.')
except TypeError:
    print("There is an error in the map file format. " \
    "Likely missing data or non-numerical values.")

#%% Photoluminescence
# Converting and analysing PL spectra

#Parsing photoluminescence data sets
PL_files = [f for f in filenames if '_PL' in f or '_pl' in f]
PL_wavenumbers = []
initialised_flag = False
#pl_flag = input('Additional PL data outputs? [plot, report, summaries]')
pl_flag = 'plot'

xb_peak = LorentzianModel(nan_policy='propagate', prefix='xb_')
xa_peak = LorentzianModel(nan_policy='propagate', prefix='xa_')
xtneg_peak = LorentzianModel(nan_policy='propagate', prefix='xtneg_')
poly_PL = PolynomialModel(0, nan_policy='propagate', prefix='PLbkgd_')
model_PL = poly_PL + xb_peak + xa_peak + xtneg_peak
PL_results = []
PL_intensites = []
PL_XB_CENTER = 2.0    # eV, starting center for xb peak
PL_XA_CENTER = 1.88   # eV, starting center for xa peak
PL_XTNEG_CENTER = 1.82  # eV, starting center for xt- peak

try:
    for PL_index, filename in enumerate(PL_files):
        current_PL_data = np.genfromtxt(os.path.join(directory, filename), dtype=None, delimiter='\t')
        PL_wavenumbers.append(current_PL_data[:,0])
        PL_split_index = np.where(current_PL_data[:,0] < 1100)[0][0]
        raman_split_Lindex = PL_split_index + np.where(current_PL_data[PL_split_index:,0] < 450)[0][0]

        [A1g_peak_index, E2g_peak_index] = find_raman_peaks(current_PL_data, 0)
        normal_intensity = current_PL_data[:,1]/(current_PL_data[A1g_peak_index,1])

        E_exc = (h*c*10**9)/lam
        E_stokes = h*c*10**2*current_PL_data[:PL_split_index,0]
        PL_energy = E_exc - E_stokes
        PL_intensity = normal_intensity[:PL_split_index]

        pars_PL = poly_PL.guess(PL_intensity, x=PL_energy)
        pars_PL.update(xb_peak.guess(PL_intensity, x=PL_energy))
        pars_PL['xb_center'].set(value=PL_XB_CENTER)
        pars_PL['xb_amplitude'].set(min=1.e-09)
        pars_PL.update(xa_peak.guess(PL_intensity, x=PL_energy))
        pars_PL['xa_center'].set(value=PL_XA_CENTER)
        pars_PL['xa_amplitude'].set(min=1.e-09)
        pars_PL.update(xtneg_peak.guess(PL_intensity, x=PL_energy))
        pars_PL['xtneg_center'].set(value=PL_XTNEG_CENTER)
        pars_PL['xtneg_amplitude'].set(min=1.e-09)
        out_PL = model_PL.fit(PL_intensity, pars_PL, x=PL_energy)

        if 'report' in pl_flag:
            print(f'PL Spectrum: {filename}')
            print(out_PL.fit_report(show_correl=False))

        if 'plot' in pl_flag:
            fig, ax = plt.subplots(figsize=(12/2.54, 9/2.54))
            ax.plot(PL_energy, PL_intensity, 'x', label='Data', markersize=3, color='grey')
            ax.plot(PL_energy, out_PL.best_fit, '-', label='Fit', linewidth=1)
            ax.plot(PL_energy, out_PL.eval_components(x=PL_energy)['xb_'], '--', label='$X_B$ Peak', linewidth=1.5)
            ax.plot(PL_energy, out_PL.eval_components(x=PL_energy)['xa_'], '--', label='$X_A$ Peak', linewidth=1.5)
            ax.plot(PL_energy, out_PL.eval_components(x=PL_energy)['xtneg_'], '--', label='$X_{T^-}$ Peak', linewidth=1.5)
            ax.set_xlabel('Energy / eV', fontname='Arial', fontsize=12)
            ax.set_ylabel('Relative Intensity', fontname='Arial', fontsize=12)
            ax.set_title(f'PL Spectrum: {filename}', fontname='Arial', fontsize=12)
            ax.tick_params(labelsize=10)
            ax.legend(fontsize=10)
            plt.tight_layout()
            plt.savefig(os.path.join(directory, str(filename) + "PL fit.png"), dpi=600)

        PL_fit_result = {"sample": filename[:-4]}
        multiple_output_PL(PL_fit_result, out_PL)
        PL_fit_result.update({"Red Chi2": out_PL.redchi})
        trion_ratio = PL_fit_result['xa_amplitude'] / PL_fit_result['xt-_amplitude']
        trion_ratio_error = trion_ratio * np.sqrt((PL_fit_result['xa_amplitude_error'] / PL_fit_result['xa_amplitude'])**2 
                                                  + (PL_fit_result['xt-_amplitude_error'] / PL_fit_result['xt-_amplitude'])**2)
        PL_fit_result.update({"Xa/Xt ratio": trion_ratio, "Xa/Xt ratio_error": trion_ratio_error})
        PL_results.append(PL_fit_result)

        PL_peak_index = np.argmax(normal_intensity[:PL_split_index])
        PL_peak_intensity = round(PL_intensity[PL_peak_index], 3)
        PL_intensites.append(PL_peak_intensity)

        if initialised_flag == False:
            normal_data = {}
            PL_data = {}
            initialised_flag = True

        sample_name = filename[:-4]
        normal_data[sample_name] = {
            "wave": current_PL_data[:,0],
            "intensity": normal_intensity
        }
        PL_data[sample_name] = {
            "energy": PL_energy,
            "intensity": PL_intensity
        }
        if 'summaries' in pl_flag:
            print('\n' \
            '______PL Spectra Stats______')
            print(f'Sample Name: {sample_name}')
            formatted_ratio = format_with_error(trion_ratio, trion_ratio_error)
            print(f'Exciton/Trion Ratio = {formatted_ratio}')
            print(f'Photoluminescence Intensity = {PL_peak_intensity}')
            formatted_xa = format_with_error(PL_fit_result['xa_center'], PL_fit_result['xa_center_error'])
            print(f'A Exciton Energy = {formatted_xa} eV')
            formatted_xb = format_with_error(PL_fit_result['xb_center'], PL_fit_result['xb_center_error'])
            print(f'B Exciton Energy = {formatted_xb} eV')
            formatted_xt = format_with_error(PL_fit_result['xt-_center'], PL_fit_result['xt-_center_error'])
            print(f'Trion Energy = {formatted_xt} eV')

    plt.figure(figsize=(12/2.54, 9/2.54))
    for count, sample_name in enumerate(PL_data.keys()):
        plt.plot(
            PL_data[sample_name]["energy"],
            PL_data[sample_name]["intensity"],
            label=f'{sample_name}',
            color=get_color_gradient("#FE9800", "#000065", len(PL_files))[count],
            linestyle='solid'
        )
    plt.title("Photoluminescence Peak", fontname='Arial', fontsize=12)
    plt.xlabel("Energy / eV", fontname='Arial', fontsize=12)
    plt.ylabel("Relative Intensity", fontname='Arial', fontsize=12)
    plt.tick_params(labelsize=10)
    plt.legend(fontsize=10, frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(directory, "PL.png"), dpi=600)

    print('______________________________')
    
    #Calculate summary statistics of the PL data
    n_results = len(PL_results)
    xa_xt_ratios = [result.get("Xa/Xt ratio", 0) for result in PL_results]
    xa_centers = [result.get("xa_center", 0) for result in PL_results]
    xb_centers = [result.get("xb_center", 0) for result in PL_results]
    xt_centers = [result.get("xt-_center", 0) for result in PL_results]
    xa_amplitudes = [result.get("xa_amplitude", 0) for result in PL_results]
    xt_amplitudes = [result.get("xt-_amplitude", 0) for result in PL_results]
    total_amplitudes = [xa + xt for xa, xt in zip(xa_amplitudes, xt_amplitudes)]

    means = {
        'xa_xt': np.mean(xa_xt_ratios),
        'xa': np.mean(xa_centers),
        'xb': np.mean(xb_centers),
        'xt': np.mean(xt_centers),
        'total_amp': np.mean(total_amplitudes),
    }
    
    stds = {k: np.std(v) / np.sqrt(n_results) for k, v in 
            [('xa_xt', xa_xt_ratios), ('xa', xa_centers), ('xb', xb_centers), ('xt', xt_centers), ('total_amp', total_amplitudes)]}
    
    formatted_ratio = format_with_error(means['xa_xt'], stds['xa_xt'])
    formatted_xa = format_with_error(means['xa'], stds['xa'])
    formatted_xb = format_with_error(means['xb'], stds['xb'])
    formatted_xt = format_with_error(means['xt'], stds['xt'])
    formatted_total_amp = format_with_error(means['total_amp'], stds['total_amp'])
    
    print('___Combined Output Stats___')
    print(f'Mean Photoluminescence Intensity = {np.mean(PL_intensites):.2f}')
    print(f'Range Photoluminescence Intensity = {np.max(PL_intensites) - np.min(PL_intensites):.2f}')
    print(f'Exciton/Trion Ratio = {formatted_ratio}')
    print(f'Range Exciton/Trion Ratio = {np.max([result.get("Xa/Xt ratio", 0) for result in PL_results]) - np.min([result.get("Xa/Xt ratio", 0) for result in PL_results]):.2f}')
    print(f'A Exciton Energy = {formatted_xa}')
    print(f'Trion Energy = {formatted_xt}')
    print(f'B Exciton Energy = {formatted_xb}')
    print(f'Xa + Xt Amplitude = {formatted_total_amp}')

    #Export photoluminescence data
    if all(np.array_equal(arr, PL_wavenumbers[0]) for arr in PL_wavenumbers):
        with open(os.path.join(directory, "normalized_spectra.csv"), "w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            header = ["Wave"] + list(normal_data.keys())
            writer.writerow(header)
            wave = next(iter(normal_data.values()))["wave"]
            for i in range(len(wave)):
                row = [wave[i]] + [normal_data[sample]["intensity"][i] for sample in normal_data]
                writer.writerow(row)

        with open(os.path.join(directory, "PL_data.csv"), "w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            header = ["Energy"] + list(PL_data.keys())
            writer.writerow(header)
            energy = next(iter(PL_data.values()))["energy"]
            for i in range(len(energy)):
                row = [energy[i]] + [PL_data[sample]["intensity"][i] for sample in PL_data]
                writer.writerow(row)
    else:
        with open(os.path.join(directory, "normalized_spectra.csv"), "w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            for sample_name in normal_data.keys():
                writer.writerow([f"{sample_name}_Wave", f"{sample_name}_Intensity"])
                wave = normal_data[sample_name]["wave"]
                intensity = normal_data[sample_name]["intensity"]
                for i in range(len(wave)):
                    writer.writerow([wave[i], intensity[i]])

        with open(os.path.join(directory, "PL_data.csv"), "w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            for sample_name in PL_data.keys():
                writer.writerow([f"{sample_name}_Energy", f"{sample_name}_Intensity"])
                energy = PL_data[sample_name]["energy"]
                intensity = PL_data[sample_name]["intensity"]
                for i in range(len(energy)):
                    writer.writerow([energy[i], intensity[i]])

    PL_fit_headings = [
        "sample",
        "xb_amplitude", "xb_amplitude_error",
        "xb_center", "xb_center_error",
        "xb_fwhm", "xb_fwhm_error",
        "xb_height", "xb_height_error",
        "xa_amplitude", "xa_amplitude_error",
        "xa_center", "xa_center_error",
        "xa_fwhm", "xa_fwhm_error",
        "xa_height", "xa_height_error",
        "xt-_amplitude", "xt-_amplitude_error",
        "xt-_center", "xt-_center_error",
        "xt-_fwhm", "xt-_fwhm_error",
        "xt-_height", "xt-_height_error",
        "Xa/Xt ratio", "Xa/Xt ratio_error", "Red Chi2"
    ]

    with open(os.path.join(directory, "PL_fitting.csv"), "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(PL_fit_headings)
        for row in PL_results:
            writer.writerow([row.get(h, '') for h in PL_fit_headings])
except NameError:
    print('No PL files found in selected directory. Ensure all files are suffixed \'_PL\'.')