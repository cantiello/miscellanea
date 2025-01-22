#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 15:28:44 2025

@author: mik

README:
This script simulates the globular cluster luminosity function (GCLF) using a Gaussian distribution. 
It includes:

1. A Gaussian distribution representing the GCLF.
2. Two scenarios for counting globular clusters (GCs):
   - Scenario 1: Count GCs brighter than a sharp magnitude cutoff (m_lim).
   - Scenario 2: Count GCs using a combination of the Gaussian distribution and a Pritchet function.

Key parameters:
- `mtom` (turnover magnitude): Mean of the Gaussian distribution.
- `sigma`: Width of the Gaussian distribution.
- `ngc`: Total number of GCs in the population.
- `m_lim`: Sharp magnitude cutoff.
- `alpha`: Parameter controlling the shape of the Pritchet function.

Outputs:
- Two plots: 
  1. Gaussian distribution with cutoff.
  2. Gaussian distribution, Pritchet function, and their product.
- Total counts for each scenario printed to the console.

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson

# Define the Gaussian function
def gaussian(x, mtom, sigma):
    """Compute the Gaussian function values."""
    prefactor = 1 / (sigma * np.sqrt(2 * np.pi))
    return prefactor * np.exp(-0.5 * ((x - mtom) / sigma) ** 2)

# Define the Pritchet function
def pritchet(alpha, m, m_lim):
    """Compute the Pritchet function values."""
    return 0.5 * (1 - alpha * (m - m_lim) / np.sqrt(1 + alpha**2 * (m - m_lim)**2))

# Parameters
mtom = 26.5       # Turnover magnitude (mean of the Gaussian)
sigma = 1.4       # Standard deviation of the Gaussian
ngc = 14000       # Total number of globular clusters
m_lim = 26        # Sharp magnitude cutoff
alpha = 1.5       # Pritchet function parameter
m = np.linspace(20, 33, 2500)  # Magnitude range for calculations

# Compute Gaussian values
y_unscaled = gaussian(m, mtom, sigma)

# Scale the Gaussian to achieve the target total population (ngc)
scaling_factor = ngc / simpson(y=y_unscaled, x=m)
y_scaled = y_unscaled * scaling_factor

# Compute total counts below the magnitude limit (m_lim) for the Gaussian
total_counts_below_mlim = int(simpson(y=y_scaled[m <= m_lim], x=m[m <= m_lim]))

# Compute Pritchet function values
f_values = pritchet(alpha, m, m_lim)

# Compute the product of the Pritchet function and the scaled Gaussian
total_product = f_values * y_scaled

# Calculate the integral under the product curve
integral_product = int(simpson(y=total_product, x=m))

# Plot settings
fig, axs = plt.subplots(1, 2, figsize=(16, 6), sharey=True)

# Plot 1: Gaussian Distribution
axs[0].plot(m, y_scaled, color='green', label='Gaussian')
axs[0].axvline(x=m_lim, color='red', linestyle='--', label='$m_{lim}$')
axs[0].axvline(x=mtom, color='blue', linestyle=':', label='$m_{TOM}$')
axs[0].fill_between(m[m <= m_lim], y_scaled[m <= m_lim], color='orange', alpha=0.3, label='Counted Region')
axs[0].text(22, 0.9 * np.max(y_scaled), f'$N_{{GC}}$ = {ngc}', color='black', fontsize=14)
axs[0].text(22, 0.8 * np.max(y_scaled), f'$m_{{TOM}}={mtom}$', color='black', fontsize=14)
axs[0].text(22, 0.7 * np.max(y_scaled), f'$\sigma_{{GCLF}}={sigma}$', color='black', fontsize=14)
axs[0].text(22, 0.6 * np.max(y_scaled), f'Tot. Counts = {total_counts_below_mlim}', fontsize=14)
axs[0].set_title('Gaussian Distribution', fontsize=14)
axs[0].set_xlabel('$I_E$ [mag]', fontsize=14)
axs[0].set_ylabel('Counts', fontsize=14)
axs[0].set_xlim(21, 32.5)
axs[0].legend()

# Plot 2: Pritchet, Gaussian, and Their Product
axs[1].plot(m, y_scaled, color='green', label='Gaussian')
axs[1].plot(m, np.max(y_scaled) * f_values, color='blue', label='Pritchet')
axs[1].plot(m, total_product, color='orange', linestyle='--', label='Product')
axs[1].axvline(x=m_lim, color='red', linestyle='--', label='$m_{lim}$')
axs[1].fill_between(m, total_product, color='orange', alpha=0.3, label='Counted Region')
axs[1].text(22, 1.2 * np.max(total_product), f'Tot. Counts = {integral_product}', fontsize=14)
axs[1].set_title('Pritchet Function, Gaussian, & Their Product', fontsize=14)
axs[1].set_xlabel('$I_E$ [mag]', fontsize=14)
axs[1].set_xlim(21, 32.5)
axs[1].legend()

# Save and show the figure
plt.tight_layout()
plt.savefig('gaussian_pritchet.jpg', dpi=300)
plt.show()

# Print the results
print(f"Total counts below m_lim={m_lim} for Gaussian: {total_counts_below_mlim}")
print(f"Integral under the product curve: {integral_product}")
