# callbacks.py #####################
import json
from dash.dependencies import Input, Output, State, ALL
from dash.exceptions import PreventUpdate
from dash import callback_context, html

from .filters import get_filters
from .optimizer import optimize_sequence

def register_callbacks(app, filters_config):
    num_filters = len(filters_config)

    ## If we have secondary structure, show it!
    @app.callback(
        Output('my-default-forna', 'sequences'),
        Input('secondary_structure', 'data'),
        State('sequence', 'data'),
    )
    def show_selected_sequences(secondary_structure, current_sequence):
        if secondary_structure is None:
            raise PreventUpdate
        return [{
            'sequence': current_sequence,
            'structure': secondary_structure
            }]

    # We include the secondary structure in the annotations_per_filter
    # So we don't need a separate callback for updating it

    @app.callback(
        Output("annotations_per_filter", "data"),
        [Input("sequence", "data")],
        prevent_initial_call=True  # We don't want this to run on load
    )
    def run_filters(seq):
        # Initialize filters
        filters_to_apply = []
        for filter_conf in filters_config:
            filter_class = filter_conf['class']
            params = filter_conf['params']
            title = filter_conf.get('title', '')
            filter_instance = filter_class(seq, title=title, **params)
            filters_to_apply.append(filter_instance)

        # Get annotations per filter
        annotations_per_filter = []
        secondary_structure = None  # To store secondary structure if any filter provides it
        for filter in filters_to_apply:
            filter.process()
            annotations = filter.get_annotations()
            annotations_per_filter.append(annotations)
            # Check if filter provides secondary structure
            if hasattr(filter, 'get_secondary_structure'):
                secondary_structure = filter.get_secondary_structure()
        return annotations_per_filter

    # New sequence -> update submission box
    @app.callback(
        Output("submission-box", "value"),
        [Input("sequence", "data")],
        prevent_initial_call=True  # We don't want this to run on load
    )
    def update_submission_box(seq):
        return seq

    # Update highlighting and suggestions
    @app.callback(
        [[Output('default-sequence-viewer-{}'.format(i), 'coverage') for i in range(num_filters)]+
         [Output('default-sequence-viewer-{}'.format(i), 'sequence') for i in range(num_filters)]+
         [Output('sidebar-content', 'children')]],
        Output('clicked_nucleotide', 'data'),
        Output('clicked_filter', 'data'),
        [[Input('default-sequence-viewer-{}'.format(i), 'mouseSelection') for i in range(num_filters)],
        Input('annotations_per_filter', 'data'),
        State('sequence', 'data'),
        State('clicked_nucleotide', 'data'),
        State('clicked_filter', 'data')],
        [Input(f'toggle-switch-{i}', 'value') for i in range(num_filters)],
    )
    def update_highlighting_and_suggestions(mouseSelections, annotations_per_filter, current_sequence, prev_nucleotide, prev_filter, *toggle_states):
        
        if current_sequence is None:
            raise PreventUpdate

        # Determine which nucleotide and filter were clicked
        ctx = callback_context
        chosen_nucleotide = None
        chosen_filter = None
        if ctx.triggered:
            triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
            if triggered_id.startswith('default-sequence-viewer'):
                i = int(triggered_id.split('-')[-1])
                mouseSelection = mouseSelections[i]
                if mouseSelection is not None:
                    chosen_nucleotide = mouseSelection['start'] - 1
                    chosen_filter = i
            else:
                # Something else was clicked
                chosen_nucleotide = prev_nucleotide
                chosen_filter = prev_filter

        # Generate coverage for each filter
        coverages = []
        for i in range(num_filters):
            annotations = annotations_per_filter[i]
            filter_conf = filters_config[i]
            filter_class = filter_conf['class']
            params = filter_conf['params']
            filter_instance = filter_class('', **params)
            coverage = sequence_coverage_from_annotations(annotations, filter_instance, chosen_nucleotide)
            coverages.append(coverage)

        # Generate sidebar suggestions
        if chosen_filter is not None:
            annotations = annotations_per_filter[chosen_filter]
            filter_conf = filters_config[chosen_filter]
            filter_class = filter_conf['class']
            params = filter_conf['params']
            filter_instance = filter_class('', **params)
            sidebar_children = sidebar_children_from_annotations(annotations, filter_instance, chosen_nucleotide)
        else:
            sidebar_children = [html.P("Suggestions will be shown here.", className="lead")]

        # Update sequences based on toggle states
        seq_list = []
        filters_enabled = [bool(state) for state in toggle_states]
        for state in filters_enabled:
            if state:
                seq_list.append(current_sequence)
            else:
                seq_list.append(' ')  # Empty string hides the sequence viewer

        return coverages + seq_list + [sidebar_children], chosen_nucleotide, chosen_filter

    def sequence_coverage_from_annotations(annotations, filter_instance, chosen_nucleotide):
        coverage = []
        for annotation in annotations:
            # Check if the chosen nucleotide is in this annotation
            underscore = False
            if chosen_nucleotide is not None:
                if annotation['start'] <= chosen_nucleotide < annotation['end']:
                    underscore = True
            coverage.append({
                'start': annotation['start'],
                'end': annotation['end'],
                'bgcolor': filter_instance.score_to_color(annotation['score']),
                'underscore': underscore
            })
        return coverage

    def sidebar_children_from_annotations(annotations, filter_instance, chosen_nucleotide):
        if chosen_nucleotide is None or filter_instance is None:
            return [html.P("Suggestions will be shown here.", className="lead")]
        else:
            # Find the annotation for the chosen nucleotide
            annotation = None
            for a in annotations:
                if a['start'] <= chosen_nucleotide < a['end']:
                    annotation = a
                    break
            if annotation is None:
                raise PreventUpdate

            suggestions = annotation['suggestions']

            # Create buttons for suggestions
            buttons = []
            for i, suggestion in enumerate(suggestions):
                suggested_string, score = suggestion
                color = filter_instance.score_to_color(score)
                buttons.append(html.Button(
                    suggested_string,
                    id={'type': 'suggestion-button', 'index': i},
                    className='btn btn-primary',
                    style={'background-color': color, 'color': 'black'}
                ))
            return buttons

    # Submit button -> sequence
    @app.callback(
        Output("sequence", "data", allow_duplicate=True),
        [Input("submit-button", 'n_clicks')],
        [State("submission-box", "value")],
        prevent_initial_call=True
    )
    def handle_submit_button(submit_button_nclicks, submitted_sequence):
        ctx = callback_context
        if not ctx.triggered:
            raise PreventUpdate
        triggered_id = ctx.triggered[0]['prop_id']
        if triggered_id == 'submit-button.n_clicks':
            # Format the sequence
            submitted_sequence = ''.join(submitted_sequence.upper().split())
            return submitted_sequence
        else:
            raise PreventUpdate

    # Suggestion buttons -> sequence
    @app.callback(
        Output("sequence", "data", allow_duplicate=True),
        [Input({'type': 'suggestion-button', 'index': ALL}, 'n_clicks')],
        [State({'type': 'suggestion-button', 'index': ALL}, 'id')],
        [State('sequence', 'data')],
        [State('annotations_per_filter', 'data')],
        [State('clicked_nucleotide', 'data')],
        [State('clicked_filter', 'data')],
        prevent_initial_call=True
    )
    def handle_suggestion_buttons(n_clicks_list, id_list, current_sequence, annotations_per_filter, chosen_nucleotide, chosen_filter):
        ctx = callback_context

        if not ctx.triggered:
            raise PreventUpdate

        # Get the index of the clicked button
        triggered = ctx.triggered[0]['prop_id'].split('.')[0]
        button_id = json.loads(triggered)['index']

        if chosen_filter is None:
            raise ValueError('No filter chosen!?')

        if n_clicks_list[button_id]:
            # Find the chosen annotation
            chosen_annotation = None
            for a in annotations_per_filter[chosen_filter]:
                if a['start'] <= chosen_nucleotide < a['end']:
                    chosen_annotation = a
                    break

            if chosen_annotation is None:
                raise PreventUpdate

            start_index = chosen_annotation['start']
            end_index = chosen_annotation['end']
            suggestion, _ = chosen_annotation['suggestions'][button_id]

            new_sequence = current_sequence[:start_index] + suggestion + current_sequence[end_index:]

            return new_sequence

        raise PreventUpdate

    # Run Filters button -> optimize sequence
    @app.callback(
        Output('sequence', 'data'),
        Input('run-filters-button', 'n_clicks'),
        State('sequence', 'data'),
        *[State(f'toggle-switch-{i}', 'value') for i in range(num_filters)]
    )
    def run_filters_button(n_clicks, current_sequence, *toggle_states):
        if n_clicks is None:
            raise PreventUpdate

        filters_enabled = [bool(state) for state in toggle_states]
        optimized_sequence = optimize_sequence(
            current_sequence,
            filters_config=filters_config,
            filters_enabled=filters_enabled
        )
        return optimized_sequence