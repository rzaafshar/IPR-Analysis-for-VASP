# IPR-Analysis-for-VASP
Automated calculation of spin-resolved normalized inverse participation ratios (IPR) from VASP PARCHG files for quantitative analysis of electronic-state localization in materials and defects.
# Batch IPR Analysis for VASP

## Overview

`batch_ipr.py` is a Python script developed to automate the calculation of the **normalized inverse participation ratio (IPR)** from VASP partial charge density (`PARCHG`) files.

The IPR provides a quantitative measure of the spatial localization of electronic states, making it useful for investigating charge localization in pristine and defect-containing materials.

The script supports spin-polarized calculations and automatically processes multiple PARCHG files across different materials and defect configurations.

## Features

- Automatically searches for `PARCHG` and `PARCHG.*` files in directories named `parch`.
- Calculates normalized real-space IPR values from three-dimensional charge density grids.
- Supports spin-resolved IPR calculations using the total charge density and magnetization density.
- Extracts selected band information from PARCHG filenames and `IBAND` settings in INCAR.
- Identifies the spin channel exhibiting stronger localization.
- Calculates the ratio between spin-up and spin-down IPR values.
- Exports the results to a CSV file for subsequent analysis.
- Records processing errors without interrupting the entire batch calculation.

## Methodology

The normalized real-space inverse participation ratio is defined as:

\[
\mathrm{IPR}^{*}=\frac{N_{\mathrm{grid}}\sum_i\rho_i^2}{(\sum_i\rho_i)^2}
\]

where \(N_{\mathrm{grid}}\) is the total number of grid points and \(\rho_i\) represents the partial charge density at grid point \(i\).

For a uniformly distributed density, the normalized IPR equals 1, whereas larger values indicate stronger spatial localization.

For spin-polarized calculations, the spin-up and spin-down densities are reconstructed from the total charge density and magnetization density.

## Requirements

- Python 3
- NumPy
- pymatgen

Install the required packages using:

`pip install numpy pymatgen`

## Usage

Run the script from the command line:

`python batch_ipr.py /path/to/calculations`

To specify an output filename:

`python batch_ipr.py /path/to/calculations -o results.csv`

If no directory is provided, the script searches the current working directory.

By default, the calculated results are saved as `ipr_results.csv`.

## Output

The generated CSV file includes:

- Material and defect names
- PARCHG filename and selected band information
- Charge density grid dimensions
- Total normalized IPR
- Spin-up and spin-down IPR values
- Maximum IPR and corresponding spin channel
- Spin-dependent IPR ratio
- Integrated spin-density sums
- Calculation status

## Applications

This script is particularly useful for high-throughput analysis of electronic-state localization in defect-containing semiconductors, enabling systematic comparisons of spatial charge localization across multiple materials, defect configurations, and spin channels.

**Note:** The calculated IPR characterizes spatial localization of the supplied partial charge density. It does not independently establish polaron formation, which requires additional physical and energetic evidence.
