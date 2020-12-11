import pandas as pd
import math, os, sys

"""
Usage: ./cd_spectra_mpl.py <file_or_folder>

Plots either a single CD spectrum or a folder of spectra.
"""

def parse_ascii(filename):

    start = 0
    xunits = None
    yunits = None
    y2units = None
    enzyme_conc = None
    with open(filename, 'r') as f:
        print('reading file ', filename)


        for index, line in enumerate(f):
            if line.startswith('XUNITS'):
                xunits = line.split()[1]
            elif line.startswith('YUNITS'):
                yunits = line.split()[1]
            elif line.startswith('Y2UNITS'):
                y2units = line.split()[1]
            elif line.startswith('XYDATA'):
                start = index + 1
            elif line.startswith('enzyme') or line.startswith('ENZYME'):
                enzyme_conc = line.split()[1]
        col_list = []
        for col in [xunits, yunits, y2units]:
            if col:
                col_list.append(col)

    data = pd.read_csv(filename,names=col_list,sep='\t',skiprows=start)

    if enzyme_conc:
        print('Normalizing to molar elipticity for ', str(filename))
        #data[yunits] = 100 * (data[yunits]/float(1000)) / ((float(enzyme_conc) *
            #float(10**-6)) * (2) )
        # coef = 0.001 / 1000 * 1000 / 10 # Coefficient that convert mDeg*L*/mol/cm to 10^3*Deg*cm^2/dmol
        path_length = 0.2 # cm
        num_aa = len_dict[name_dict[os.path.basename(filename)]]
        # data['Molar Elipticity'] = coef * data[yunits] / (float(enzyme_conc) * 10**-6 ) / float(0.2)
        data['Molar Elipticity'] = data[yunits] / (float(enzyme_conc) *
                10**-6 * num_aa * path_length * 10 * 1000)
    else:
        data['Molar Elipticity'] = data[yunits]

    return pd.melt(data,id_vars=[yunits,y2units,'Molar Elipticity'])

def collect_spectra(folder):
    filepaths = []
    if os.path.isdir(folder):
        for file in os.listdir(folder):
            if file.split('.')[-1] == 'txt':
                filepaths.append(os.path.join(folder,file))
    elif os.path.isfile(folder):
        filepaths.append(folder)

    data = pd.DataFrame()
    labels = []
    for f in filepaths:
        if f.endswith('.txt'):
            df = parse_ascii(f)
            df['filename'] = f
            labels.append(f.split('/')[-1])
            data = pd.concat([data,df])

    return data, labels


def theta(T, Tm, dH, R):
    # Assume molecularity of 1 for now
    R = .001987203611
    x = (dH / R) ((1 / T) - (1 / Tm))
    
    psi = 1 / (1 + math.exp(x))

    """
    For molecularity of 2, the equation would be
    1 - (e**x)/4) (sqrt(1 + 8 e**-x) - 1)
    """

    return psi

from uuid import uuid4
from matplotlib import pyplot as plt
import scipy.optimize as opt

name_dict = {
        '2018-12-05_wt_spectrum_corrected.txt': 'Wild-Type',
        '2018-12-06_B3_spectrum_corrected.txt': 'V2D9r',
        '2019-03-27_e38d_spectrum_25c_corrected': 'V2D9r E38D',
        '2019-04-25_e38a_spectrum_25c_corrected': 'V2D9r E38A',
        '2019-08-17_lima_e38d_corrected': 'V1D8r E38D',
        '2018-12-05_wt_spectrum_corrected.txt': 'Wild-Type KSI',
        'lima_25c_corrected.txt': 'V1D8r',
        'lima_e38d_25c_corrected.txt': 'V1D8r E38D',
        '2018-12-06_e38d_spectrum_corrected.txt': 'V2D9r E38D',
        '2019-04-25_e38a_spectrum_25c_corrected.txt': 'V2D9r E38A',
        '2020-08-12_lima_25c_2.txt': 'V1D8r',
        'lima_melt_222_corrected.txt': 'V1D8r',
        'lima_e38d_melt_222_corrected.txt': 'V1D8r E38D',
        'lima_melt_corrected.txt': 'V1D8r',
        '2018-12-06_B3_melt_222_corrected.txt': 'V2D9r',
        '2018-12-06_e38d_melt_222_corrected.txt': 'V2D9r E38D',
        '2019-03-27_e38d_melt_222_corrected.txt': 'V2D9r E38D',
        '2019-03-21_e38a_melt_corrected.txt': 'V2D9r E38A',
        '2019-04-25_e38a_tempscan_corrected.txt': 'V2D9r E38A',
        '2020-11-19_lima_melt_222_corrected.txt': 'V1D8r',
        '2020-12-04-b3_scan_corrected.txt': 'V2D9r',
        '2020-11-19_lima_melt_222_corrected.txt': 'V1D8r',
        '2020-12-04_b3_melt_222_corrected.txt': 'V2D9r (new)',
        '2020-11-19_lima_scan_corrected.txt': 'V1D8r',
        }

color_dict = {
        'Wild-Type KSI': 'green',
        'V2D9r': 'darkorange',
        'V2D9r (new)': 'black',
        'V1D8r': 'blue',
        'V2D9r E38D': 'sandybrown',
        'V2D9r E38A': 'peachpuff',
        'V1D8r E38D': 'skyblue',
        }

len_dict = {
        'Wild-Type KSI': 127,
        'V2D9r': 127,
        'V2D9r (new)': 127,
        'V1D8r': 126,
        'V2D9r E38D': 127,
        'V2D9r E38A': 127,
        'V1D8r E38D': 126,
        }

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']


def update_spectra_graph(data):
    df = data[data['variable']=='NANOMETERS']

    traces = []
    i = 0
    for name, group in df.groupby(['filename']):
        #print(group['Molar Elipticity'])
        print(name)
        points = plt.plot(
                group['value'],
                group['Molar Elipticity'],
                label=name_dict[name.split('/')[-1]],
                # color=color_dict[name_dict[name.split('/')[-1]]],
                color='black'
        )
        # plt.legend()
        traces.append(points)
        i += 1
    return traces

def update_melt_graph(data):
    df = data[data['variable']=='Temperature']

    traces =[]
    for name, group in df.groupby(['filename']):
        points = go.Scatter(
                x = group['value'],
                y = group['Molar Elipticity'],
                mode = 'markers',
                name=name.split('/')[-1],
                color=color_dict[name_dict[name.split('/')[-1]]],
        )
        #optimizedParameters, pcov = opt.curve_fit(theta,
                #group['variable'], group['Molar Elipticity'])
        traces.append(points)
    return traces


if __name__ == '__main__':
    data, labels = collect_spectra(sys.argv[1])
    update_spectra_graph(data)
    plt.xlabel('Wavelength ($nm$)')
    plt.ylabel('Mean residue ellipticity ($10^3$ $deg$ $cm^2$ $dmol^{-1}$)')
    plt.show()
