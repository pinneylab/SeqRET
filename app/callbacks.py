import json
from dash.dependencies import Input, Output, State, ALL
from dash.exceptions import PreventUpdate
from filters import filters_to_apply
from dash import callback, callback_context, html
from dash.exceptions import PreventUpdate

## If we have secondary structure, show it!
@callback(
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

#if we have a filter that does secondary structure, update ours:
@callback(
    Output('secondary_structure', 'data'),
    Input('annotations_per_filter', 'data')
)
def update_secondary_structure(annotations_per_filter):
    for filter in filters_to_apply:
        if filter.get_title() == 'RNA Secondary Structure':
            return filter.get_secondary_structure()
    return None

@callback(
    Output("annotations_per_filter", "data"),
    [Input("sequence", "data")],
    prevent_initial_call=True #we don't want this to run on load, because we haven't input a sequence yet.
)
def run_filters(seq):
    #filter_list is static and passed in to wrapper function
    
    #get output of each filter:
    annotations_per_filter = []
    for filter in filters_to_apply:
        filter.update_sequence(seq, force=True)
        filter.process()
        annotations = filter.get_annotations()
        annotations_per_filter.append(annotations)
    return annotations_per_filter

@callback(
    [[Output('default-sequence-viewer-{}'.format(i), 'coverage') for i in range(len(filters_to_apply))]+
    [Output('default-sequence-viewer-{}'.format(i), 'sequence') for i in range(len(filters_to_apply))]+
    [Output('sidebar-content', 'children')]],
    Output('clicked_nucleotide', 'data'),
    Output('clicked_filter', 'data'),
    [[Input('default-sequence-viewer-{}'.format(i), 'mouseSelection') for i in range(len(filters_to_apply))],
    Input('annotations_per_filter', 'data'),
    State('sequence', 'data'),
    State('clicked_nucleotide', 'data'),
    State('clicked_filter', 'data')
    ]
)
def update_highlighting_and_suggestions(mouseSelections, annotations_per_filter, current_sequence, prev_nucleotide, prev_filter):
    
    #if current sequence hasn't been assigned, don't update. Otherwise our sequences break.
    if current_sequence is None:
        raise PreventUpdate
    
    #get the nucleotide + filter that were just clicked on 
    ctx = callback_context
    chosen_nucleotide = None
    chosen_filter = None
    if ctx.triggered:
        triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
        #check the id starts with 'default-sequence-viewer'
        if triggered_id.startswith('default-sequence-viewer'):
            #get the index of the filter that was clicked on
            i = int(triggered_id.split('-')[-1])
            #get the mouse selection
            mouseSelection = mouseSelections[i]
            if mouseSelection is not None:
                chosen_nucleotide = mouseSelection['start']-1
                chosen_filter = i
        else:
            #something else was clicked on
            chosen_nucleotide = prev_nucleotide
            chosen_filter = prev_filter
    
    #for each filter, get the coverage to display in sequence viewer
    coverages = []
    for filter in filters_to_apply:
        #get the annotation
        annotations = filter.get_annotations()
        #get the coverage
        coverages.append(sequence_coverage_from_annotations(annotations, filter, chosen_nucleotide))
    
    annotations = None
    if chosen_filter is not None:
        annotations = annotations_per_filter[chosen_filter]
    sidebar_children = sidebar_children_from_annotations(annotations, filter, chosen_nucleotide)

    #update sequence as well while we're doing this
    sequences = [current_sequence]*len(filters_to_apply)

    return coverages + sequences + sidebar_children, chosen_nucleotide, chosen_filter

def sequence_coverage_from_annotations(annotations, filter, chosen_nucleotide):
    '''
    Input:
        annotations: a list of dicts. Each dict is of the form:
            { 'start': start_index,
              'end': end_index,
              'score': score,
              'suggestions': [ (suggestion_1, score_1), (suggestion_2, score_2), ... ]}
    Output:
        coverage: a list of dicts. Each dict is of the form:
            { 'start': start_index,
              'end': end_index,
              'bgcolor': color,
              'underscore': True/False }
    '''
    coverage = []
    for annotation in annotations:
        #is chosen_nucleotide in this annotation?
        underscore = False
        if chosen_nucleotide is not None:
            if annotation['start'] <= chosen_nucleotide < annotation['end']:
                #yes, so underscore it
                underscore = True
        coverage.append({
            'start': annotation['start'],
            'end': annotation['end'],
            'bgcolor': filter.score_to_color(annotation['score']),
            'underscore': underscore
        })
    return coverage

def sidebar_children_from_annotations(annotations, filter, chosen_nucleotide):
    '''
    Input:
        annotations: a list of dicts. Each dict is of the form:
            { 'start': start_index,
              'end': end_index,
              'score': score,
              'suggestions': [ (suggestion_1, score_1), (suggestion_2, score_2), ... ]}
        filter: the filter object
        chosen_filter: the index of the chosen filter
    Output:
        sidebar_children: a list of html elements to be displayed in the sidebar.
    '''
    if chosen_nucleotide is None:
        return [html.P(
            "Suggestions will be shown here.", className="lead"
        )]
    else:
        #get the annotation for our chosen nucleotide
        annotation = None
        for a in annotations:
            if a['start'] <= chosen_nucleotide < a['end']:
                annotation = a
                break
        if annotation is None:
            raise PreventUpdate
        #get the suggestions
        suggestions = annotation['suggestions']

        #make buttons
        buttons = []
        for i, suggestion in enumerate(suggestions):
            suggested_string, score = suggestion
            color = filter.score_to_color(score)
            buttons.append(html.Button(suggested_string, id={'type': 'suggestion-button', 'index': i}, className='btn btn-primary', style={'background-color': color, 'color': 'black'}))

        return [buttons]

#CGCTGGCGGCGTGCTAATACATGCAGGAGGTGCNGGATCGATGGGAGCTTGCTCCCTGCAGATCAGCGGCGGACGGGTGAGTAACACGTGGGTAACCTGCCTGTAAGACTGGGATAACTCCGGGAAACCGGGGCTAATACCGGATAATTTAGTTCCTCGCATGAGGAACTGTTGAAAGGTGGCTTC

#Submit button -> sequence
@callback(
    Output("sequence", "data", allow_duplicate=True),
    [Input("submit-button", 'n_clicks')],
    [State("submission-box", "value")],
    prevent_initial_call=True
)
def handle_submit_button(submit_button_nclicks, submitted_sequence):
    ctx = callback_context
    triggered_id = ctx.triggered[0]['prop_id']
    if triggered_id == 'submit-button.n_clicks':
        return submitted_sequence
    else:
        raise PreventUpdate

#suggestion buttons   -> sequence
# + annotations
# + current selection
@callback(
    Output("sequence", "data", allow_duplicate=True),
    Output("submission-box", "value"),
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

    button_id = json.loads(ctx.triggered[0]['prop_id'].split('.')[0])['index']

    if chosen_filter is None:
        raise ValueError('No filter chosen!?')

    if n_clicks_list[button_id] is not None:
        # Logic for chosen suggestion
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

        return new_sequence, new_sequence

    raise PreventUpdate


# @callback(
#     [Output("sequence", "data")],
#     [Input({'type': 'suggestion-button', 'index': ALL}, 'n_clicks')],
#     [Input("submit-button", 'n_clicks')],
#     [State("submission-box", "value")],
#     [State({'type': 'suggestion-button', 'index': ALL}, 'id')],
#     [State('sequence', 'data')],
#     [State('annotations_per_filter', 'data')],
#     [[State('default-sequence-viewer-{}'.format(i), 'mouseSelection') for i in range(len(filters_to_apply))]]
# )
# def handle_button_clicks(n_clicks_list, submit_button_nclicks, submitted_sequence, id_list,  current_sequence, annotations_per_filter, mouseSelections):

#     ctx = callback_context
#     triggered_id = ctx.triggered_id
#     if triggered_id == 'submit-button':
#         return [submitted_sequence]

#     #otherwise, it's the suggestion buttons
    
#     if not ctx.triggered:
#         button_id = 'No clicks yet'
#         raise PreventUpdate
#     else:
#         button_id = json.loads(ctx.triggered[0]['prop_id'].split('.')[0])['index']

#     #get index
#     chosen_nucleotide = None
#     chosen_filter = None
#     for i, mouseSelection in enumerate(mouseSelections):
#         if mouseSelection is not None:
#             chosen_nucleotide = mouseSelection['start']-1 #not 0 indexed lol
#             chosen_filter = i
#             break
#     if chosen_filter is None:
#         raise ValueError('No filter chosen!?')

#     if n_clicks_list[button_id] is not None:
#         #from the button ID, we know which suggestion was chosen.
#         #we get the start and end sequence index from that suggestion.
#         #finally, we insert the suggestion into the sequence.

#         #get the index of annotation for our chosen nucleotide:
#         chosen_annotation = None
#         for a in annotations_per_filter[chosen_filter]:
#             if a['start'] <= chosen_nucleotide < a['end']:
#                 chosen_annotation = a
#                 break
#         if chosen_annotation is None:
#             raise PreventUpdate

#         #chosen_annotation = annotations_per_filter[chosen_filter][chosen_nucleotide]
#         start_index = chosen_annotation['start']
#         end_index = chosen_annotation['end']
#         suggestion, _ = chosen_annotation['suggestions'][button_id]

#         new_sequence = current_sequence[:start_index] + suggestion + current_sequence[end_index:]

#         return [new_sequence]
#     #else it hasn't been pushed yet
    
#     raise PreventUpdate