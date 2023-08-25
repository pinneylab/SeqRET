import urllib.request as urlreq
import dash_bootstrap_components as dbc
from dash import Dash, dcc, html, Input, Output, State, ALL, callback, callback_context
import dash_bio as dashbio
from dash_bio.utils import protein_reader

from layouts import submission_box, sidebar, content
from filters import filters_to_apply
#from callbacks import run_filters, update_highlighting_and_suggestions, handle_button_clicks
from callbacks import run_filters, update_highlighting_and_suggestions, handle_submit_button, handle_suggestion_buttons

from Bio import SeqIO

import json


def create_sequence_viewer_app(starting_sequence, filter_list):
    """
    Creates a dash app that displays a number of stacked sequence viewers, given the number of filters.
    Returns the app object and a list of IDs for the sequence viewers.
    """

    #app = Dash(__name__)
    app = Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

    app.layout = html.Div(
            [dcc.Location(id="url"), sidebar, content, submission_box,
            # store the working sequence
            dcc.Store(id='sequence'),
            # store the previous sequence 
            dcc.Store(id='previous_sequence'),
            #dcc.Store(id='freq_list_per_filter'),
            #dcc.Store(id='suggestion_list_per_filter'),
            dcc.Store(id='annotations_per_filter'),
            dcc.Store(id='clicked_nucleotide'),
            dcc.Store(id='clicked_filter'),
            dcc.Store(id='secondary_structure')
            ]
        )
    
    return app
    

def runner(fasta):

    # Read in the fasta file
    with open('AM231703.1.fasta', 'r') as f:
        fasta_str = f.read()

    seq = protein_reader.read_fasta(datapath_or_datastring=fasta_str, is_datafile=False)[0]['sequence']

    # detect open reading frames to preserve:
    # assume ORF 0
    #ORF = 0

    app = create_sequence_viewer_app(seq, filters_to_apply)

    app.run(debug=True, port=8050)

runner('AM231703.1.fasta')