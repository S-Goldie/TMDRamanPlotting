#Nora Campoli, July 21, 2025
import os
import numpy as np
import matplotlib.pyplot as plt
from FunctionBank import RamanPeakInfo, get_color_gradient

h = 4.1356692*10**-15 #eV s
c = 3.0*10**8 #m/s
lam = 532 #nm

#Change here to alter file locations
directory = r"C:\Users\sgold\Documents\Waterloo PostDoc\THA Intercalation\Raman\MoS2_1TControl"
filenames = [f for f in os.listdir(directory) if f.endswith('.txt')]

plt.figure("Raman")
for i, filename in enumerate(filenames):
    #Keep in mind, the wavenumber goes from highest to lowest. Indices might need to be reversed.
    data_path = os.path.join(directory, filename)
    data = np.genfromtxt(data_path, dtype=None, delimiter='\t')

    #Cutting off the PL array by the wavenumber
    PL_split_index = np.where(data[:,0] < 1100)[0][0]
    #Cutting off the Raman array by the wavenumber. Starting later to speed up the search
    raman_split_Lindex = PL_split_index + np.where(data[PL_split_index:,0] < 450)[0][0]
    raman_split_Uindex = raman_split_Lindex + np.where(data[raman_split_Lindex:,0] < 350)[0][0]

    #Finding Raman peaks using a function (passing in a truncated list):
    [A1g_peak_index, E2g_peak_index, peak_distance] = RamanPeakInfo(data[raman_split_Lindex:raman_split_Uindex], raman_split_Lindex)

    #Normalizing array based on A1g peak. Used elsewhere so it covers the full range
    normal_intensity = data[:,1]/(data[A1g_peak_index,1])
    #For the Raman graph:
    raman_wave = data[raman_split_Lindex:raman_split_Uindex,0]
    raman_intensity = normal_intensity[raman_split_Lindex:raman_split_Uindex]

    #Converting units & finding PL in units of eV:
    E_exc = (h*c*10**9)/lam
    E_stokes = h*c*10**2*data[:PL_split_index,0]
    #For the PL graph:
    PL_energy = E_exc - E_stokes
    PL_intensity = normal_intensity[:PL_split_index]
    #PL peak finding:
    PL_peak_index = np.argmax(normal_intensity[:PL_split_index])
    PL_peak_energy = round(PL_energy[PL_peak_index], 6)
    PL_peak_intensity = round(PL_intensity[PL_peak_index], 3)

    #Assembling the output files:
    if filename == filenames[0]:
        normal_data = np.vstack((data[:,0], normal_intensity))
        PL_data = np.vstack((PL_energy, PL_intensity))
        data_summary = np.array([filename, data[E2g_peak_index,0], data[A1g_peak_index,0], peak_distance, PL_peak_energy])
    else:
        normal_data = np.vstack((normal_data, data[:,0], normal_intensity))
        PL_data = np.vstack((PL_data, PL_energy, PL_intensity))
        data_summary_temp = np.array([filename, data[E2g_peak_index,0], data[A1g_peak_index,0], peak_distance, PL_peak_energy])
        data_summary = np.vstack((data_summary, data_summary_temp))
    
    #Plotting:
    plt.plot(raman_wave, raman_intensity, label=r'$ \alpha $ =' + f' {filename[:-4]}', color=get_color_gradient("#FE9800", "#000065", len(filenames))[i], linestyle='solid')



#Raman Peaks:
plt.title("Raman Peaks")
plt.legend()
plt.xlabel("Wave/$cm^{-1}$")
plt.ylabel("Intensity")

#PL Peaks:
plt.figure("Photoluminescence")
for count, i in enumerate(range(0,2*len(filenames),2)):
    plt.plot(PL_data[i,:], PL_data[i+1,:], label=r'$ \alpha $ =' + f' {filenames[count][:-4]}', color=get_color_gradient("#FE9800", "#000065", len(filenames))[count], linestyle='solid')
plt.title("Photoluminescence Peak")
plt.legend()
plt.xlabel("Wave/$eV$")
plt.ylabel("Intensity")
#Labels for the extrema. Could put these into a list during the loop. Maybe Later.
# plt.text(PL_peak_energy, PL_peak_intensity, f"({PL_peak_energy}eV, {PL_peak_intensity})", fontsize=8)

#Full Spectra:
plt.figure("Full Spectra")
for count, i in enumerate(range(0,2*len(filenames),2)):
    plt.plot(normal_data[i,:], normal_data[i+1,:], label=r'$ \alpha $ =' + f' {filenames[count][:-4]}', color=get_color_gradient("#FE9800", "#000065", len(filenames))[count], linestyle='solid')
plt.title("Normalized Spectra")
plt.legend()
plt.xlabel("Wave/$cm^{-1}$")
plt.ylabel("Intensity")
#Labels for the extrema. Could put these into a list during the loop. Maybe Later.
# plt.text(E2g_peak_wavenumber, normal_intensity[E2g_peak_index], str(E2g_peak_wavenumber) + "$cm^{-1}$", fontsize=8)
# plt.text(A1g_peak_wavenumber, normal_intensity[A1g_peak_index], str(A1g_peak_wavenumber) + "$cm^{-1}$", fontsize=8)
# plt.text(data[PL_peak_index,0], PL_peak_intensity, str(round(data[PL_peak_index,0], 2)) + "$cm^{-1}$", fontsize=8)

plt.show()

#Writing the outputs:
# os.chdir('../')
# try:
#     os.mkdir("Raman Spectra CSV Files")
#     print(f"Directory '{directory_name}' created successfully.")
# except FileExistsError:
#     print(f"Directory '{directory_name}' already exists.")

np.savetxt(f"{directory}{filename[:-4]}_normalized.csv", normal_data.T, delimiter=',', header='Wave,Intensity,'*len(filenames), comments='')
np.savetxt(f"{directory}{filename[:-4]}_PL.csv", PL_data.T, delimiter=',', header='Energy,Intensity,'*len(filenames), comments='')
np.savetxt(f"{directory}{filename[:-4]}_Summary.csv", data_summary, fmt='%s', delimiter=',', header='Sample Name,E2g Peak/cm-1,A1g Peak/cm-1,Raman Peak Difference/cm-1,Photoluminescence Peak/eV', comments='')
