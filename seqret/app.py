import dash_bootstrap_components as dbc
from dash import Dash, dcc, html

from .layouts import submission_box, sidebar, content
from .filters import get_filters
from .callbacks import run_filters, update_highlighting_and_suggestions, handle_submit_button, handle_suggestion_buttons


def create_sequence_viewer_app(filter_list):
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
    

def start_app():
    ### Run the app ###
    app = create_sequence_viewer_app(get_filters())
    from waitress import serve
    print('Server started.')
    print('Please navigate to http://127.0.0.1:8050 in your web browser.')
    serve(app.server, host="0.0.0.0", port=8050)

if __name__ == "__main__":
    start_app()