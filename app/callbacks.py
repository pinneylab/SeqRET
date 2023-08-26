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

# #if our filter toggles are clicked, hide the corresponding sequence viewer
# @callback(
#     [Output(f'default-sequence-viewer-{i}', 'sequence') for i in range(len(filters_to_apply))],
#     State('sequence', 'data'),
#     [Input(f'toggle-switch-{i}', 'value') for i in range(len(filters_to_apply))],
#     allow_duplicate=True
# )
# def toggle_sequence_viewers(sequence, *toggle_states):
#     # Print the state of each toggle button
#     seq_list = []
#     for i, state in enumerate(toggle_states):
#         #if state is 'on':
#         if not state:
#             #hide the sequence viewer
#             seq_list.append(None)
#         else:
#             #show the sequence viewer
#             seq_list.append(sequence)
#     return seq_list


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
    State('clicked_filter', 'data')],
    [Input(f'toggle-switch-{i}', 'value') for i in range(len(filters_to_apply))],
)
def update_highlighting_and_suggestions(mouseSelections, annotations_per_filter, current_sequence, prev_nucleotide, prev_filter, *toggle_states):
    
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
    print('test1')
    #sidebar_children = None
    if chosen_filter is not None:
        annotations = annotations_per_filter[chosen_filter]
        sidebar_children = sidebar_children_from_annotations(annotations, filters_to_apply[chosen_filter], chosen_nucleotide)
    else:
        sidebar_children = sidebar_children_from_annotations(annotations, None, chosen_nucleotide)

    print('test2')
    #update sequence as well while we're doing this

    #sequences = [current_sequence]*len(filters_to_apply)
    ### TODO: make filters invisible when not toggled!
    seq_list = []
    print('toggle_states', toggle_states)
    for i, state in enumerate(toggle_states):
        seq_list.append([current_sequence])
    print(seq_list)

    return coverages + seq_list + sidebar_children, chosen_nucleotide, chosen_filter

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
    print('filter', filter)
    print('chosen_nucleotide', chosen_nucleotide)
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
        #make sequence uppercase:
        submitted_sequence = submitted_sequence.upper()
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

### When the user clicks the "Run Filters" button, run all selected filters jointly.
@callback(
    Output('sequence', 'data'),
    Input('run-filters-button', 'n_clicks'),
    State('annotations_per_filter', 'data'),
    State('sequence', 'data'),
    *[State(f'toggle-switch-{i}', 'value') for i in range(len(filters_to_apply))]
)
def run_filters(n_clicks, annotations_per_filter, current_sequence, *toggle_states):
    if n_clicks is None:
        raise PreventUpdate

    filters_to_run = []
    # Print the state of each toggle button
    for i, state in enumerate(toggle_states):
        #if state is 'on':
        if state: #an unpressed toggle button is just an empty list
            filters_to_run.append(i)
    
    #get the annotations from just the filters we care about:
    annotations_for_chosen_filters = []
    for i in filters_to_run:
        annotations_for_chosen_filters.append(annotations_per_filter[i])
    
    #for each codon in the sequence, get the list of suggestions for each chosen filter.
    #for each suggested codon, we'll add up the scores for each filter. If a suggested codon is not in the list of suggestions for a filter, we'll add a score of 0 for that one.
    #then we'll choose the codon with the highest score.

    #get the sequence
    sequence = current_sequence

    #get the list of codons
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]

    #this reorganizes our data. Now we will have a list of length len(sequence)/3,
    #where each element is a list of length len(filters_to_run), where each element is a list of suggestions for that codon.
    #the suggestions are of the form (suggestion, score)
    #we'll also add up the scores for the CURRENT codons while we're here.
    suggestions_for_each_codon = []
    current_scores_for_each_codon = []
    for i in range(0, len(sequence), 3):

        start_index = i
        end_index = i+3
        suggestions_for_codon = []
        current_codon_score = 0
        for annotations in annotations_for_chosen_filters:
            for annotation in annotations:
                if annotation['start'] == start_index and annotation['end'] == end_index:
                    suggestions_for_codon.append(annotation['suggestions'])
                    current_codon_score += annotation['score']
                    break #assuming there's only one annotation per codon
        suggestions_for_each_codon.append(suggestions_for_codon)
        current_scores_for_each_codon.append(current_codon_score)

    #now for each codon, we'll make a dictionary that contains the suggested alternate codons and their scores, by adding up the scores for each filter.
    #the keys will be the suggested codons, and the values will be the scores.
    #we'll then choose the codon with the highest score.
    best_codons = []
    for i, suggestions_for_current_codon in enumerate(suggestions_for_each_codon):
        #make a dictionary of the form {suggested_codon: score}
        codon_scores = {}
        for suggestions in suggestions_for_current_codon:
            for suggestion, score in suggestions:
                if suggestion not in codon_scores:
                    codon_scores[suggestion] = score
                else:
                    codon_scores[suggestion] += score
        #now get the codon with the highest score
        #if the dictionary is empty, just keep the current codon
        if len(codon_scores) == 0:
            best_codons.append(codons[i])
            continue
        best_codon = max(codon_scores, key=codon_scores.get)
        #if it's the same score or lower than our current codon, don't change it
        if codon_scores[best_codon] <= current_scores_for_each_codon[i]:
            best_codon = codons[i]
        best_codons.append(best_codon)

    #get the new sequence with the best codons
    new_sequence = ''
    for codon in best_codons:
        new_sequence += codon

    return new_sequence