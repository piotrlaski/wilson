import math
import os
import numpy as np
import matplotlib.pyplot as plt
import re
from utils import *
import configparser


config = configparser.ConfigParser()
config.read('config.ini')

cif_file_path = config['data_files']['cif_file_path']
ratios_file_path = config['data_files']['ratios_file_path']

plot = config['plot_settings']['plot']
image_output_path = config['plot_settings']['image_output_path']
plot_X_axis_max = float(config['plot_settings']['plot_X_axis_max'])

export_csv = config['exported_csv_settings']['export_csv']
export_csv_output_dir = config['exported_csv_settings']['export_csv_output_dir']

if not os.path.isdir(image_output_path) and plot:
    os.mkdir(image_output_path)
kb, del_T = calc_kb_del_T(cif_file_path, ratios_file_path, plot, image_output_path, plot_X_axis_max, export_csv, export_csv_output_dir)
print(f'Kb = {kb:.2f}, dT = {del_T:.2f} K')
