import dash_bootstrap_components as dbc
from dash import Dash, dcc, html
import os

from .layouts import submission_box, sidebar, content
from .filters import get_filters
from .callbacks import register_callbacks

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

def create_sequence_viewer_app():
    """
    Creates a Dash app that displays sequence viewers for each filter.
    """
    app = Dash('SeqRET', title='SeqRET', external_stylesheets=[dbc.themes.BOOTSTRAP], assets_folder=os.path.join(MODULE_DIR, 'assets'))

    app.layout = html.Div(
        [dcc.Location(id="url"), sidebar, content, submission_box,
         dcc.Store(id='sequence'),
         dcc.Store(id='previous_sequence'),
         dcc.Store(id='annotations_per_filter'),
         dcc.Store(id='clicked_nucleotide'),
         dcc.Store(id='clicked_filter'),
         dcc.Store(id='secondary_structure')]
    )
    return app

def start_app(host='127.0.0.1', port=8050):
    app = create_sequence_viewer_app()
    filters_config = get_filters()
    register_callbacks(app, filters_config)
    print(f'Server started at http://{host}:{port}')
    app.run_server(host=host, port=port, debug=True)

if __name__ == "__main__":
    start_app()