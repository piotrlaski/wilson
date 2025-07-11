# PhotoWilson

PhotoWilson is a tool for calculating photoinduced temperature changes in crystals after photoexcitation based on laser ON/laser OFF intensity ratios of indexed X-ray reflections. 
It is the public version of the project "TRL_Photowilson", which is referenced in my thesis. 


## Getting Started

Clone the repository into your local space (or download all .py files into an empty directory). When you are ready, run `wilson.py` to calculate the temperature changes and generate analysis plots.

## Required Input Files

You will need two input files:

1. **A .cif file** containing crystallographic information including:
   - Unit cell parameters (a, b, c, α, β, γ)
   - Temperature of the measurement
   - Anisotropic atomic displacement parameters (Uani values)

2. **A .hkl file** containing reflection intensity ratios with columns for:
   - Miller indices (h, k, l)
   - Intensity ratios (laser ON/laser OFF) for each triplet of indices

## Usage

The main calculation is performed by calling the `calc_kb_del_T()` function. Simply modify the file paths in `wilson.py` to point to your input files:

```python
cif_file_path = r'path_to_your_file.cif'
ratios_file_path = r'path_to_your_file.hkl'

kb, del_T = calc_kb_del_T(cif_file_path, ratios_file_path, plot=True)
```

## Output

The tool calculates and returns:

1. **Kb value** - The temperature scaling factor related to the structural changes
2. **ΔT (del_T)** - The photoinduced temperature change in Kelvin

When `plot=True` is specified, the tool generates a PhotoWilson scatter plot showing the relationship between (sin θ/λ)² and ln(intensity ratios), saved as a .png file in the same directory as your input ratios file.

## Methodology

The calculation uses the Wilson plot method to analyze the temperature-dependent changes in X-ray reflection intensities. The tool:

- Extracts crystallographic parameters from the .cif file
- Calculates d-spacings for each reflection using Miller indices
- Performs regression fitting of ln(intensity ratios) vs. (sin θ/λ)² data
- Determines temperature changes based on changes in atomic displacement parameters

An in-depth explanation of the procedure can be found in the following paper:
Schmokel, M. S., Kaminski, R., Benedict, J. B. & Coppens, P. (2010). Acta Cryst. A66, 632-636.
https://doi.org/10.1107/S0108767310029429

## File Structure

- `wilson.py` - Main execution file
- `utils.py` - Contains all calculation functions including CIF parsing, d-spacing calculations, and temperature analysis

For any additional help or bug reports, please contact pa.laski@uw.edu.pl
