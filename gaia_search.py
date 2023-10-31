#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 10:43:37 2023


This is a simple Python script designed to query 
the Gaia Data Release 2 (DR2) database for astronomical data.
It performs two types of searches: a cone search and a box search. 
The script allows users to specify central coordinates, search radius
(for the cone search), and search width and height (for the box search)
to retrieve astronomical data from the Gaia DR2 catalog.

@author: mik
"""

# Import required libraries
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
import astropy.units as u

# Define the central coordinates, width, and height
print('You need to change the input coordinates for the search')
cra = 54.023  # Center RA
cdec = -35.269  # Center Dec
srad = 0.7  # Search radius
swi = 0.5  # Box search width
shei = 0.5  # Box search height

# Central coordinates and the search radius
central_coord = SkyCoord(ra=cra, dec=cdec, unit=(u.degree, u.degree), frame='icrs')
search_radius = srad * u.deg

# Cone search query
query = f"""
SELECT *
FROM gaiadr2.gaia_source
WHERE 1=CONTAINS(POINT(ra, dec), 
                CIRCLE({central_coord.ra.deg}, {central_coord.dec.deg}, {search_radius.to(u.degree).value}))
"""

# Modify the output filename to include RA, DEC, and RAD
output_cone = f"gaia_RA{cra}_DEC{cdec}_RAD{srad}.fits"

# Launch the query job 
job = Gaia.launch_job_async(query, name="ConeSearchQuery", dump_to_file=True, output_file=output_cone)

# Wait for the job to complete
job.get_results()

# Access the query results
results = job.results

# You can now work with the query results
print(results)

print('IN CASE ONE WANTS THE BOX SEARCH, use the lines below')
'''
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
import astropy.units as u

# Define the central coordinates, width, and height
central_coord = SkyCoord(ra=cra, dec=cdec, unit=(u.degree, u.degree), frame='icrs')
width = swi * u.deg
height = shei * u.deg

# Create a box search query
query = f"""
SELECT *
FROM gaiadr2.gaia_source
WHERE (ra >= {central_coord.ra.value - width.value}) AND (ra <= {central_coord.ra.value + width.value})
  AND (dec >= {central_coord.dec.value - height.value}) AND (dec <= {central_coord.dec.value + height.value})
"""

# Modify the output filename to include cra, cdec, and srad
output_box = f"gaia_RA{cra}_DEEC{cdec}_W{swi}_H{shei}.fits"


# Launch the query job asynchronously
job = Gaia.launch_job_async(query, name="BoxSearchQuery",  dump_to_file=True, output_file=output_box)

# Wait for the job to complete
job.get_results()

# Access the query results
results = job.results

# You can now work with the query results
print(results)
'''
