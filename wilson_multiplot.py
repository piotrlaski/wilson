import math
import os
import numpy as np
import matplotlib.pyplot as plt
import re
from utils import *

cif_file_path = r'C:\Users\piotr\Documents\VS_Code\wilson\example_files\example_multiwilson\snp_220K.cif'
ratios_folder_path = r'C:\Users\piotr\Documents\VS_Code\wilson\example_files\example_multiwilson\ratios'  ## this has to contain only ratio files, but can be split under multiple subfolders (hkls)
output_csv_file = r'C:\Users\piotr\Documents\VS_Code\wilson\example_files\example_multiwilson\output.csv'

## get all ratios files absolute paths (ignore possible .png from previous runs)
ratio_paths = []
for root, dirs, files in os.walk(ratios_folder_path):
    for file in files:
        abs_file = os.path.join(root,file)
        if abs_file[-4:] != '.png' and os.path.isfile(abs_file) and abs_file != output_csv_file:
            ratio_paths.append(abs_file)
        
## run wilson for all ratios in the selected folder, write the results in output_csv_file
for file in ratio_paths:
    kb, del_T = calc_kb_del_T(cif_file_path, file, plot = True)
    print(f'Kb = {kb:.2f}, dT = {del_T:.2f} K')
    with open(output_csv_file, 'a') as f:
        f.write(f'{os.path.split(file)[1]} {kb:.2f} {del_T:.2f}\n')
    print(f'PhotoWilson scatter plot created at {os.path.abspath(os.path.splitext(file)[0]) + ".png"}')