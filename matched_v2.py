#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 21 00:11:08 2023

This Python script is designed for cross-matching and combining astronomical 
data from FITS tables. It provides a flexible function called 'cross_match_and_save' t
hat takes two input FITS tables, cross-matches them based on specified coordinate
columns, and saves the combined results to an output FITS table. The script
also allows custom subscripts, matching radii, and optional third input tables 
for more advanced use cases.

@author: mik
"""
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table, hstack

def cross_match_and_save(table1_path, table2_path, output_path, matching_radius=0.5 * u.arcsecond, subscript1="_X1", subscript2="_X2", table3=None, match_columns1=None, match_columns2=None):
    # Read the two input FITS tables
    table1 = Table.read(table1_path)
    table2 = Table.read(table2_path)


    # Create SkyCoord objects for each table using the specified column names for matching
    coords1 = SkyCoord(table1[match_columns1[0]], table1[match_columns1[1]], unit=(u.deg, u.deg))
    coords2 = SkyCoord(table2[match_columns2[0]], table2[match_columns2[1]], unit=(u.deg, u.deg))

    # Perform the cross-matching to find the best matches
    idx1, d2d, _ = coords1.match_to_catalog_sky(coords2)
    best_matches = d2d < matching_radius

    # Create tables of best matches for each table
    matched_table1 = table1[best_matches]
    matched_table2 = table2[idx1[best_matches]]

    # Calculate mean alpha and delta for the output table3
    mean_alpha = 0.5 * (matched_table1[match_columns1[0]] + matched_table2[match_columns2[0]])
    mean_delta = 0.5 * (matched_table1[match_columns1[1]] + matched_table2[match_columns2[1]])

    # Rename columns with user-defined subscripts
    for colname in matched_table1.colnames:
        matched_table1[colname].name = colname + subscript1

    for colname in matched_table2.colnames:
        matched_table2[colname].name = colname + subscript2

    # Calculate angular separation in arcseconds
    separation = coords1[best_matches].separation(coords2[idx1[best_matches]]).to(u.arcsecond)

    # Add the angular separation column to matched_table1
    matched_table1['Separation' + subscript2] = separation
    
    # Stack the matched tables horizontally to keep all columns
    combined_table = hstack([matched_table1, matched_table2])

    # Add mean alpha and delta columns to the combined table
    combined_table['RA'] = mean_alpha
    combined_table['Dec'] = mean_delta

    
    # If a table3 is provided, append the results to it
    if table3 is not None:
        table3 = hstack([table3, combined_table])
    else:
        table3 = combined_table

    # Save the combined table to the output file
    table3.write(output_path, format='fits', overwrite=True)

    # Return the updated table3 for further use
    return table3

# Example usage with custom subscripts and match column names
table3 = cross_match_and_save('/Users/mik/Desktop/Astrowork/EUCLID/EUCLID_ERO/FORNAX/VIS/AL/resid_v1_cat.fits',
                             '/Users/mik/Desktop/Astrowork/EUCLID/EUCLID_ERO/FORNAX/NISP/y_v1_cat.fits',
                             'output_table.fits',
                             subscript1="_VIS",
                             subscript2="_Y",
                             match_columns1=['ALPHA_J2000', 'DELTA_J2000'],
                             match_columns2=['ALPHA_J2000', 'DELTA_J2000'])



# Example usage with custom subscripts and match column names
table4 = cross_match_and_save('output_table.fits',
                             '/Users/mik/Desktop/Astrowork/EUCLID/EUCLID_ERO/FORNAX/NISP/j_v1_cat.fits',
                             'output_table2.fits',
                             subscript1="",
                             subscript2="_J",
                             match_columns1=['RA', 'Dec'],
                             match_columns2=['ALPHA_J2000', 'DELTA_J2000'])



# Example usage with custom subscripts and match column names
table5 = cross_match_and_save('output_table2.fits',
                             '/Users/mik/Desktop/Astrowork/EUCLID/EUCLID_ERO/FORNAX/NISP/h_v1_cat.fits',
                             'VIS_YJH.fits',
                             subscript1="",
                             subscript2="_H",
                             match_columns1=['RA', 'Dec'],
                             match_columns2=['ALPHA_J2000', 'DELTA_J2000'])




# Example usage with custom subscripts and match column names
table5 = cross_match_and_save('VIS_YJH.fits',
                             'mastergc_FDS.fits',
                             'VIS_YJH_ref.fits',
                             subscript1="",
                             subscript2="_ref",
                             match_columns1=['RA', 'Dec'],
                             match_columns2=['RAJ2000', 'DEJ2000'])




# Example usage with custom subscripts and match column names
table6 = cross_match_and_save('/Users/mik/Desktop/Astrowork/EUCLID/EUCLID_ERO/FORNAX/VIS/AL/resid_v1_cat.fits',
                             'mastergc_FDS.fits',
                             'VIS_ref.fits',
                             subscript1="_VIS",
                             subscript2="_ref",
                             match_columns1=['ALPHA_J2000', 'DELTA_J2000'],
                             match_columns2=['RAJ2000', 'DEJ2000'])


#MATCH ME WITH TEYMOOR
# Example usage with custom subscripts and match column names
table3 = cross_match_and_save('/Users/mik/Desktop/Astrowork/EUCLID/EUCLID_ERO/FORNAX/VIS/AL/resid_v1_cat.fits',
                             '/Users/mik/Desktop/Astrowork/EUCLID/EUCLID_ERO/FORNAX/teymur_forced_merged_v2.fits',
                             'mik_ts.fits',
                             subscript1="_VIS",
                             subscript2="_TS",
                             match_columns1=['ALPHA_J2000', 'DELTA_J2000'],
                             match_columns2=['RA', 'DEC'])




#MATCH ME WITH TEYMOOR
# Example usage with custom subscripts and match column names
table3 = cross_match_and_save('/Users/mik/Desktop/Astrowork/EUCLID/EUCLID_ERO/FORNAX/VIS/AL/VIS_YJH.fits',
                             '/Users/mik/Desktop/Astrowork/EUCLID/EUCLID_ERO/FORNAX/teymur_forced_merged_v2.fits',
                             'mik_ts_nisp.fits',
                             subscript1="_MC",
                             subscript2="_TS",
                             match_columns1=['RA', 'Dec'],
                             match_columns2=['RA', 'DEC'])






#MATCH LSB and Flattened h
# Example usage with custom subscripts and match column names
table3 = cross_match_and_save('/Users/mik/Desktop/Astrowork/EUCLID/EUCLID_ERO/FORNAX//NISP/h_v1_cat.fits',
                             '/Users/mik/Desktop/Astrowork/EUCLID/EUCLID_ERO/FORNAX/NISP/h_v1_cat_lsb.fits',
                             'h_flat_lsb.fits',
                             subscript1="_F",
                             subscript2="_L",
                             match_columns1=['ALPHA_J2000', 'DELTA_J2000'],
                             match_columns2=['ALPHA_J2000', 'DELTA_J2000'])



#MATCH LSB and flattened j
# Example usage with custom subscripts and match column names
table3 = cross_match_and_save('/Users/mik/Desktop/Astrowork/EUCLID/EUCLID_ERO/FORNAX//NISP/j_v1_cat.fits',
                             '/Users/mik/Desktop/Astrowork/EUCLID/EUCLID_ERO/FORNAX/NISP/j_v1_cat_lsb.fits',
                             'j_flat_lsb.fits',
                             subscript1="_F",
                             subscript2="_L",
                             match_columns1=['ALPHA_J2000', 'DELTA_J2000'],
                             match_columns2=['ALPHA_J2000', 'DELTA_J2000'])


#MATCH lsb and flattened y
# Example usage with custom subscripts and match column names
table3 = cross_match_and_save('/Users/mik/Desktop/Astrowork/EUCLID/EUCLID_ERO/FORNAX//NISP/y_v1_cat.fits',
                             '/Users/mik/Desktop/Astrowork/EUCLID/EUCLID_ERO/FORNAX/NISP/y_v1_cat_lsb.fits',
                             'y_flat_lsb.fits',
                             subscript1="_F",
                             subscript2="_L",
                             match_columns1=['ALPHA_J2000', 'DELTA_J2000'],
                             match_columns2=['ALPHA_J2000', 'DELTA_J2000'])


# GAIA

#Select the output I need
tab_gaia=Table.read('/Users/mik/Desktop/Astrowork/EUCLID/EUCLID_ERO/FORNAX/gaia_RA54.023_DEEC-35.269_W0.5_H0.5.fits')
selected_columns = tab_gaia['source_id', 'ra', 'dec', 'parallax', 'parallax_error', 'phot_g_mean_mag', 'phot_g_mean_flux_error', 'bp_rp', 'bp_g', 'a_g_val','e_bp_min_rp_val' ]
selected_columns.write('gaia_sel.fits', format='fits', overwrite=True)


table5 = cross_match_and_save('/Users/mik/Desktop/Astrowork/EUCLID/EUCLID_ERO/FORNAX/gaia_sel.fits',
                             '/Users/mik/Desktop/Astrowork/EUCLID/EUCLID_ERO/FORNAX/VIS/AL/VIS_YJH.fits',
                             'gaia_euclid.fits',
                             subscript1="",
                             subscript2="",
                             match_columns1=['ra', 'dec'],
                             match_columns2=['RA', 'Dec'])

