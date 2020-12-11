import pandas as pd
import math, os, sys

"""
Usage: ./cd_spectra.py <folder>

Plots a folder of CD spectra and melting curves and sets up a local
server to display them in an interactive browser window.
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
        coef = 0.001 / 1000 * 1000 / 10 # Coefficient that convert mDeg*L*/mol/cm to 10^3*Deg*cm^2/dmol
        data['Molar Elipticity'] = coef * data[yunits] / (float(enzyme_conc) * 10**-6 ) / float(2)
    else:
        data['Molar Elipticity'] = data[yunits]

    return pd.melt(data,id_vars=[yunits,y2units,'Molar Elipticity'])

def collect_spectra(folder):
    filepaths = []
    for file in os.listdir(folder):
        if file.split('.')[-1] == 'txt':
            filepaths.append(os.path.join(folder,file))

    data = pd.DataFrame()
    for f in filepaths:
        if f.endswith('.txt'):
            df = parse_ascii(f)
            df['filename'] = f
            data = pd.concat([data,df])

    return data


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

import dash, dash_table
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
from itertools import cycle
from flask_caching import Cache
from uuid import uuid4
import scipy.optimize as opt


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

CACHE_CONFIG= {
        'CACHE_TYPE': 'simple',
        }
cache= Cache()
cache.init_app(app.server, config=CACHE_CONFIG)

data = collect_spectra(sys.argv[1])

def update_spectra_graph(data):
    df = data[data['variable']=='NANOMETERS']

    traces = []
    for name, group in df.groupby(['filename']):
        points = go.Scatter(
                x = group['value'],
                y = group['Molar Elipticity'],
                mode='lines',
                name=name.split('/')[-1]
        )

        traces.append(points)
    return traces

def update_melt_graph(data):
    df = data[data['variable']=='Temperature']

    traces =[]
    for name, group in df.groupby(['filename']):
        points = go.Scatter(
                x = group['value'],
                y = group['Molar Elipticity'],
                mode = 'markers',
                name=name.split('/')[-1]
        )
        #optimizedParameters, pcov = opt.curve_fit(theta,
                #group['variable'], group['Molar Elipticity'])
        traces.append(points)
    return traces

def serve_layout():
    session_id = str(uuid4())

    return html.Div([
        html.Div([

            dcc.Graph(id='spectra',
                figure={
                    'data':update_spectra_graph(data),
                    'layout':go.Layout(
                        xaxis={
                            'title':'Wavelength (nm)'
                        },
                        yaxis={'title':'Molar elipticity'},
                        margin={'l':40,'b':40,'t':100,'r':10},
                        hovermode='closest',
                        title='CD Spectra',
                    )
                }
                ),

            dcc.Graph(id='melt',
                figure={
                    'data':update_melt_graph(data),
                    'layout':go.Layout(
                        xaxis={'title':'Temperature (C)'},
                        yaxis={'title':'Molar elipticity'},
                        margin={'l':40,'b':40,'t':100,'r':10},
                        hovermode='closest',
                        title='Melting curves',
                    )
                }
                ),
        ])
    ])

app.layout = serve_layout

if __name__ == '__main__':
    app.run_server(debug=True)
