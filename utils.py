import math
import os
import numpy as np
import matplotlib.pyplot as plt
import re

def d_calc_helpers(angles: list, lengths: list) -> list:
    ## unwrap cell params
    alpha, beta, gamma = [math.radians(i) for i in angles]
    a, b, c = lengths
    ## calculate helper expressions
    S11 = (b**2) * (c**2) * (math.sin(alpha)**2)
    S22 = (a**2) * (c**2) * (math.sin(beta)**2)
    S33 = (a**2) * (b**2) * (math.sin(gamma)**2)
    S12 = a * b * (c**2) * (math.cos(alpha) * math.cos(beta) - math.cos(gamma))
    S23 = (a**2) * b * c * (math.cos(beta) * math.cos(gamma) - math.cos(alpha))
    S13 = a * (b**2) * c * (math.cos(gamma) * math.cos(alpha) - math.cos(beta))
    V = a * b * c * math.sqrt(1 - math.cos(alpha)**2 - math.cos(beta)**2 - math.cos(gamma)**2 + (2 * math.cos(alpha) * math.cos(beta) * math.cos(gamma)))
    return [S11, S22, S33, S12, S23, S13, V]
    
def d_calc(helpers: list, millers: list) -> float:
    ## unwrap millers and helpers
    h, k, l = millers
    S11, S22, S33, S12, S23, S13, V = helpers
    ## calculate interplanar distance
    inv_d_sqr = (1 / V**2 ) * ((S11 * h**2) + (S22 * k**2) + (S33 * l**2) + (2 * S12 * h * k) + (2 * S23 * k * l) + (2 * S13 * h * l))
    d = 1 / math.sqrt(inv_d_sqr)
    return d

def import_ratios(file: os.PathLike) -> np.array:
    R = np.genfromtxt(file)
    return R

def hkl_ratios_to_sinthetasq_logr(angles: list, lengths: list, file: os.PathLike) -> np.array:
    ## prepare helpers, import data, prepare output array
    helpers = d_calc_helpers(angles, lengths)
    in_array = import_ratios(file)
    out_array = np.empty([0,2])
    ## loop over all relfections, calculate (sinth/lam)^2 and ln(R) for each
    for i in in_array[:,:4]:
        millers = i[:3]
        ratio = i[3]
        log_ratio = math.log(ratio)
        d = d_calc(helpers, millers)
        sint_lam_sq = (1 / (2*d))**2
        out_array = np.append(out_array, [[sint_lam_sq, log_ratio]], axis = 0)
    return out_array

def linfit(angles: list, lengths: list, ratio_file: os.PathLike, plot: bool=False, plot_output_dir: os.PathLike=None, plot_X_axis_max = 0.5, export_csv_graph=False, export_csv_output_dir: os.PathLike=None) -> tuple:
    ratio_file_name = os.path.splitext(os.path.split(ratio_file)[1])[0]
    if plot_output_dir is None and plot is True:
        print (f'!Plot_output_path not specified, plot created in ratio file location!')
        plot_output_dir = os.path.split(ratio_file)[0]
    if export_csv_output_dir is None and export_csv_graph is True:
        print (f'!export_csv_graph not specified, plot exported data created in ratio file location!')
        export_csv_output_dir = os.path.split(ratio_file)[0]

    ## do a linear regression of ln(R) = a * (sinth/lam)^2 + b
    out_array = hkl_ratios_to_sinthetasq_logr(angles, lengths, ratio_file)
    out_array_trans = np.transpose(out_array)
    x, y = out_array_trans
    a, b = np.polyfit(x, y, deg = 1)
    ## export a csv file of the plot
    if export_csv_graph:
        csv_name = os.path.join(export_csv_output_dir, ratio_file_name) + '.csv'
        np.savetxt(csv_name, out_array, delimiter=' ')
    ## scatter plot with the fitted line
    if plot:
        plt.scatter(x,y, s=1)
        plt.axline(xy1=(0, b), slope=a, color='r', label=f'$y = {a:.2f}x {b:+.2f}$')
        plt.legend(loc='upper right')
        plt.xlim(0, plot_X_axis_max)
        plt.ylim(-1, 1)
        plot_name = os.path.join(plot_output_dir, ratio_file_name) + '.png'
        plt.savefig(plot_name)
        plt.close()
    return (a,b)

def read_cif(file: os.PathLike) -> tuple[list, list, float]:
    with open(file) as f:
        cif_txt = f.read()

    ## regex extraction of all Uani vals (U_ani has to have uncertainty!!!)
    uani_vals_rx = re.findall(r'\S*\(\d\) Uani', cif_txt)
    uani_vals = [float(i.split(r'(')[0]) for i in uani_vals_rx]
    uani_avg = sum(uani_vals) / len(uani_vals)
    
    ## twofold regex extraction of unit cell params
    a = re.findall(r'cell\_length\_a.*', cif_txt)
    a = re.findall(r'[ \t]+[\.0-9]*', a[0])[0].split()[0]

    b = re.findall(r'cell\_length\_b.*', cif_txt)
    b = re.findall(r'[ \t]+[\.0-9]*', b[0])[0].split()[0]

    c = re.findall(r'cell\_length\_c.*', cif_txt)
    c = re.findall(r'[ \t]+[\.0-9]*', c[0])[0].split()[0]

    alpha = re.findall(r'cell\_angle\_alpha.*', cif_txt)
    alpha = re.findall(r'[ \t]+[\.0-9]*', alpha[0])[0].split()[0]

    beta = re.findall(r'cell\_angle\_beta.*', cif_txt)
    beta = re.findall(r'[ \t]+[\.0-9]*', beta[0])[0].split()[0]

    gamma = re.findall(r'cell\_angle\_gamma.*', cif_txt)
    gamma = re.findall(r'[ \t]+[\.0-9]*', gamma[0])[0].split()[0]

    angles = [float(alpha), float(beta), float(gamma)]
    lengths = [float(a), float(b), float(c)]

    ## twofold regex extraction of temperature
    temp = re.findall(r'cell\_measurement\_temperature.*', cif_txt)
    temp = re.findall(r'[ \t]+[\.0-9]*', temp[0])[0].split()[0]
    temp = float(temp)

    return(angles, lengths, uani_avg, temp)

def kb_calc(avg_Uani: float, del_B: float) -> float:
    kb = 1 + del_B / (8 * (math.pi**2) * avg_Uani)
    return (kb)

def del_t_calc(avg_Uani: float, del_B: float, T_off: float) -> float:
    B_off = 8 * (math.pi**2) * avg_Uani
    del_T = (((B_off + del_B) / B_off) - 1) * T_off
    return (del_T)

def calc_kb_del_T(cif_file: os.PathLike, ratios_file: os.PathLike, plot: bool=False, plot_output_dir: os.PathLike=None, plot_X_axis_max = 0.5, export_csv_graph = False, export_csv_output_dir: os.PathLike=None) -> tuple[float, float]:
    angles, lengths, avg_Uani, temp = read_cif(cif_file)
    a, b = linfit(angles, lengths, ratios_file, plot, plot_output_dir, plot_X_axis_max, export_csv_graph, export_csv_output_dir)
    del_B = -a / 2
    kb = kb_calc(avg_Uani, del_B)
    del_T = del_t_calc(avg_Uani, del_B, temp)
    return (kb, del_T)
