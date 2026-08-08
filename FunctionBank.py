#Nora Campoli, July 21, 2025
import numpy as np
import math

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

def RamanPeakInfo(data, start_index):
    #Keep in mind, the wavenumber goes from highest to lowest
    A1g_index_1 = np.where(data[:,0] < 415)[0][0]
    A1g_index_2 = A1g_index_1 + np.where(data[A1g_index_1:,0] < 395)[0][0]
    A1g_peak_index = A1g_index_1 + np.argmax(data[A1g_index_1:A1g_index_2, 1])

    #We already have a starting point from the first peak's bounds
    E2g_index_1 = A1g_index_2
    E2g_index_2 = E2g_index_1 + np.where(data[E2g_index_1:,0] < 375)[0][0]
    E2g_peak_index = E2g_index_1 + np.argmax(data[E2g_index_1:E2g_index_2, 1])

    #Finding distance between peaks.
    peak_distance = round(data[A1g_peak_index,0] - data[E2g_peak_index,0], 6)

    return [start_index + A1g_peak_index, start_index + E2g_peak_index, peak_distance]


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
    peak_output(destination, source, "lam")

def multiple_output_PL(destination, source, use_trion=True):
    """stores fitted outputs for each PL peak in turn"""
    peak_output(destination, source, 'xb')
    peak_output(destination, source, 'xa')
    if use_trion:
        xtneg_result = {}
        peak_output(xtneg_result, source, 'xtneg')
        for key in list(xtneg_result.keys()):
            destination[key.replace('xtneg', 'xt-', 1)] = xtneg_result[key]
    else:
        destination.update({
            'xt-_amplitude': np.nan, 'xt-_amplitude_error': np.nan,
            'xt-_center': np.nan, 'xt-_center_error': np.nan,
            'xt-_fwhm': np.nan, 'xt-_fwhm_error': np.nan,
            'xt-_height': np.nan, 'xt-_height_error': np.nan
        })

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
    if math.isnan(abs_err) or math.isnan(value):
        return "N/A"

    exponent = math.floor(math.log10(abs_err))
    first_digit = int(abs_err / 10**exponent)
    sig_figs = 2 if first_digit == 1 else 1

    rounded_error = round(abs_err, -(exponent - (sig_figs - 1)))
    decimal_places = max(0, -(exponent - (sig_figs - 1)))
    rounded_value = round(value, decimal_places)

    return f"{rounded_value:.{decimal_places}f} ± {rounded_error:.{decimal_places}f}"
