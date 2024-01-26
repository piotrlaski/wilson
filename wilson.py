import math
import os
import numpy as np
import matplotlib.pyplot as plt
import re
from utils import *

cif_file_path = r'.\example_files\[Rh(4-Br-SA)(CO2)].cif'
ratios_file_path = r'.\example_files\\glued_sortaved.hkl'

kb, del_T = calc_kb_del_T(cif_file_path, ratios_file_path, plot = True)
print(f'Kb = {kb:.2f}, dT = {del_T:.2f} K')
print(f'PhotoWilson scatter plot created at {os.path.abspath(os.path.splitext(ratios_file_path)[0]) + ".png"}')