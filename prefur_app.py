import dash
from dash import html, dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output
from prefur import kinetics
from prefur import thermo
import numpy as np
import warnings

# create app object
warnings.filterwarnings('ignore')
external_stylesheets = [dbc.themes.BOOTSTRAP]
mathjax = 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/MathJax.js?config=TeX-MML-AM_CHTML'
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.scripts.append_script({ 'external_url' : mathjax })
server = app.server

# set app layout
app.layout = html.Div([
    html.H1('PREFUR Webserver'),
    dcc.Markdown('''
        PREFUR (PREdiction of Folding and Unfolding Rates) is an 
        open-source tool for finding folding and unfolding rates
        given a protein's structural class and its length in residues.

        The three structural classes are:
        * a - alpha class meaning that majority of the protein structure is alpha helix.
        * b - beta class meaning that majority of the potein structure is beta sheets.
        * ab - mixed class meaning there is a mix of both topologies in the protein.

        Sequence length will be determined by the amino acid sequence given to the server.
    '''),
    html.Br(),
    html.Label('Structure Class'),
    html.Br(),
    dcc.Dropdown(
        id='struct-class',
        options=[{'label': type, 'value': type} for type in ['a', 'b', 'ab']],
    ),
    html.Br(),
    html.Label('Sequence'),
    html.Br(),
    dcc.Input(id='protein-seq', type='text', min=1),
    html.Br(),
    html.Div(id='kf_ku'),
    html.Br(),
    dcc.Markdown('''
        ## References
        Please cite the original PREFUR paper:

        * De Sancho, D. and Muñoz, V. Integrated prediction of protein 
        folding and unfolding rates from only size and structural class, 
        Phys. Chem. Chem. Phys. 13, 17030-17043 (2011).
          [DOI](http://dx.doi.org/10.1039/C1CP20402E)
    '''),
    html.Footer(
        'Algorithm developed by David De Sancho and Victor Muñoz.\nApp developed by Mohammad Abdulqader'
    ),
], style={'padding': '10px'})

# generate call back app


@app.callback(
    Output('kf_ku', 'children'),
    Input('protein-seq', 'value'),
    Input('struct-class', 'value'))
def prefur_predict(seq, struct, temp=298.):
    """
    Function for estimating folding and unfolding rates from protein size
    and structural type
    """
    #### set up fucntions ####
    DH = {'a': [2.15, 4.82], 'b': [1.31, 5.3], 'ab': [1.5, 5.21]}
    seq_alpha = 'ACDEFGHIKLMNPQRSTVWY'

    ### Handle invalid input. ###
    if not seq:
        return 'Please input an amino acid sequence.'
    if not struct:
        return 'Please select a structure class.'
    if any(a not in seq_alpha for a in seq.upper()):
        return 'Seqeunce contains invalid amino acids.'

    ### Once inputs valid calculate ###
    nres = len(seq)
    if nres == 0:
        return 'Please input an amino acid sequence.'

    FES = thermo.FES(nres)
    FES.gen_enthalpy_global(DHloc=DH[struct][0], DHnonloc=DH[struct][1])
    FES.gen_free(temp=temp)
    bf, bu = thermo.barrier(FES.DG)
    kf = kinetics.rates(barrier=bf, nres=nres, temp=temp)
    ku = kinetics.rates(barrier=bu, nres=nres, temp=temp)
    return f'''
        Folding Rate: ${np.round(kf, 4)}$ $s^{-1}$
        
        Unfolding Rate: ${np.round(ku, 4)}$ $s^{-1}$

        '''


if __name__ == '__main__':
    app.run_server(debug=False)
