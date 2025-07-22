#Nora Campoli, July 21, 2025
import numpy as np

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

def hex_to_RGB(hex_str):
    """ #FFFFFF -> [255,255,255]"""
    #Pass 16 to the integer function for change of base
    return [int(hex_str[i:i+2], 16) for i in range(1,6,2)]
 
def get_color_gradient(c1, c2, n):
    """
    Given two hex colors, returns a color gradient
    with n colors.
    """
    assert n > 1
    c1_rgb = np.array(hex_to_RGB(c1))/255
    c2_rgb = np.array(hex_to_RGB(c2))/255
    mix_pcts = [x/(n-1) for x in range(n)]
    rgb_colors = [((1-mix)*c1_rgb + (mix*c2_rgb)) for mix in mix_pcts]
    return ["#" + "".join([format(int(round(val*255)), "02x") for val in item]) for item in rgb_colors]