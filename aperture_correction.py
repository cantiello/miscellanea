#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 10:31:44 2023

This Python script performs various data selection and analysis tasks on a FITS 
table and generates several plots to visualize the data.

The main interest is to analyze the curve of growth using a sextractor FITS catalog.

@author: mik
"""
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import time
from astropy.stats.funcs import mad_std



####################################
#####GENERIC SETTING FOR PLOTS
font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 26}


legendfont = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 20}


plt.rc('font', **font)


#Start the time and see how much time does it take 
start_time = time.time()
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
#SELECTION PARAMETERS
minmag=18.25 #brightest mag to be used
maxmag=23.5 #faintest magnitude to be used
binac=0.4 #magnitude binning
cslim=0.95 #class star limit for selecting STAR candidatess
nac=13 #num ber of apertures
# Load the table
table_path = '/Users/mik/Desktop/Astrowork/EUCLID/EUCLID_ERO/FORNAX/VIS/AL/resid_v1_cat.fits'
tab = Table.read(table_path)


# Select objects
selected_objects = (tab['MAG_AUTO'] >= minmag) & (tab['MAG_AUTO'] <= maxmag) & (tab['CLASS_STAR'] >= cslim) & (tab['FLAGS'] == 0) & (tab['FWHM_IMAGE'] <= 2.4) & (tab['FWHM_IMAGE'] >= 1.2)


#Preliminary plots
f = plt.figure(figsize=(22, 22))
f.subplots_adjust(wspace=0.3, left=0.15, right=0.85, hspace=0.25, bottom=0.15, top=  0.85)


### FWHM
f.add_subplot(221)
# Plot FWHM_IMAGE versus MAG_AUTO for all objects
plt.scatter(tab['MAG_APER'][:,1], tab['FWHM_IMAGE'], marker='.', s=1, label='All Objects', alpha=0.5)
# Plot FWHM_IMAGE versus MAG_AUTO for the selected objects
plt.scatter(tab['MAG_APER'][:,1][selected_objects], tab['FWHM_IMAGE'][selected_objects], marker='.', s=10,color='red', label='Selected Objects (MAG_AUTO < 23 and MAGERR_AUTO > 0)', alpha=0.5)
# Add labels and legend
plt.xlabel('$m_{8pix}~[mag]$')
plt.ylabel('FWHM [pix]')
plt.xlim(17.5,25)
plt.ylim(0,10)
fw_med=np.median(tab['FWHM_IMAGE'][selected_objects])
fw_mad=mad_std(tab['FWHM_IMAGE'][selected_objects])
fw_std=np.std(tab['FWHM_IMAGE'][selected_objects])

print('Median/STDDEV/MAD FWHM = %5.2f, %5.2f, %5.2f ' %(fw_med,fw_std,fw_mad))
#plt.legend()

### CLASS_STAR
f.add_subplot(222)
# Plot FWHM_IMAGE versus MAG_AUTO for all objects
plt.scatter(tab['MAG_APER'][:,1], tab['CLASS_STAR'], marker='.', s=1, label='All Objects', alpha=0.5)
# Plot FWHM_IMAGE versus MAG_AUTO for the selected objects
plt.scatter(tab['MAG_APER'][:,1][selected_objects], tab['CLASS_STAR'][selected_objects], marker='.', color='red', label='Selected Objects (MAG_AUTO < 23 and MAGERR_AUTO > 0)', alpha=0.5)

# Add labels and legend
plt.xlabel('$m_{8pix}~[mag]$')

plt.ylabel('CLASS_STAR')
plt.xlim(17.5,25)
plt.ylim(-.05,1.05)
#plt.legend()


# #Aperture correction
ac=[]
print('Warning the a.c. diameter need to be inserted by hand')
rac=[4,8,16,24,32,40,48,56,72,88,104,120,152]


# Initialize a list to store ac values for each MAG_AUTO value
ac_values = []


# Iterate over MAG_AUTO values between 19 and 25
for mag_auto_value in np.arange(minmag, maxmag, binac):
    selected_objects = (tab['MAG_AUTO'] >= minmag) & (tab['MAG_AUTO'] <= mag_auto_value) & (tab['CLASS_STAR'] >= cslim) & (tab['FLAGS'] == 0)  & (tab['FWHM_IMAGE'] <= fw_med+3*fw_mad) & (tab['FWHM_IMAGE'] >= fw_med-3*fw_mad)

    m8 = tab['MAG_APER'][:, 1][selected_objects]

    ac = []  # Initialize ac for the current MAG_AUTO value

    for i in range(nac):
        mag_aper_i = tab['MAG_APER'][:, i][selected_objects]
        diff = m8 - mag_aper_i
        median_diff = np.median(diff)
        ac.append(median_diff)
    print('%5.2f<= mag<= %5.2f Nsel=%5i' % (minmag, mag_auto_value, len(mag_aper_i)))

    # Append the ac values for the current MAG_AUTO value to the list
    ac_values.append(ac)


# Now you have a list of ac values for each MAG_AUTO value

f.add_subplot(223)

# Plot ac values for each MAG_AUTO value
for i, mag_auto_value in enumerate(np.arange(minmag, maxmag,binac)):
    formatted_mag_auto_value = round(mag_auto_value, 2)  # Rounds to 2 decimal places
    plt.scatter(rac, ac_values[i], label=f'm_faint={formatted_mag_auto_value}')
    plt.plot(rac, ac_values[i])

plt.xlabel('diameter (pixels)')
plt.ylabel('$m_{8pix}-m_{diam}$')
plt.legend(fontsize=12)
plt.ylim(-0.3,.3)
plt.axhline(y=.13, color='black', linestyle='--', label='Zero Line')
plt.fill_between(rac, 0.11, 0.15, color='gray', alpha=0.25, label='Shaded Area')



# Initialize a list to store ac values for each MAG_AUTO value
ac2_values = []


# Iterate over MAG_AUTO values between 19 and 25
for mag_auto_value in np.arange(minmag, maxmag, binac):
    selected_objects = (tab['MAG_AUTO'] >= minmag) & (tab['MAG_AUTO'] <= mag_auto_value) & (tab['CLASS_STAR'] >= cslim) & (tab['FLAGS'] == 0)  & (tab['FWHM_IMAGE'] <= fw_med+3*fw_mad) & (tab['FWHM_IMAGE'] >= fw_med-3*fw_mad)

    ac2 = []  # Initialize ac for the current MAG_AUTO value

    for i in range(nac):
        mag_aper_i = tab['MAG_APER'][:, i][selected_objects]
        mag_aper_i_1 = tab['MAG_APER'][:, i-1][selected_objects]
        diff2 =  mag_aper_i_1-mag_aper_i
        median_diff = np.median(diff2)
        ac2.append(median_diff)
    print('%5.2f<= mag<= %5.2f Nsel=%5i' % (minmag, mag_auto_value, len(mag_aper_i)))

    # Append the ac values for the current MAG_AUTO value to the list
    ac2_values.append(ac2)


# Now you have a list of ac values for each MAG_AUTO value

f.add_subplot(224)

# Plot ac values for each MAG_AUTO value
for i, mag_auto_value in enumerate(np.arange(minmag, maxmag,binac)):
    plt.scatter(rac, ac2_values[i], label=f'MAG_AUTO={mag_auto_value}')
    plt.plot(rac, ac2_values[i])

plt.xlabel('diameter (pixels)')
plt.ylabel('$m_{(diam-1)}-m_{diam}$')
#plt.legend(fontsize='small')
plt.ylim(-0.035,.035)
plt.axhline(y=.00, color='black', linestyle='--', label='Zero Line')

plt.savefig('pippo.jpeg', bbox_inches='tight', dpi=300)

