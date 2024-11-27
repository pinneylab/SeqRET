# optimizer.py #################################
from .filters import get_filters

def optimize_sequence(sequence, filters_config=None, filters_enabled=None):
    """
    Optimizes the given sequence using the provided filters.

    Parameters:
    - sequence: the nucleotide sequence to optimize
    - filters_config: list of filter configurations. Each configuration is a dict with 'class', 'params', and 'title'.
    - filters_enabled: list of booleans indicating whether each filter is enabled

    Returns:
    - optimized_sequence: the optimized nucleotide sequence
    """
    if filters_config is None:
        filters_config = get_filters()

    if filters_enabled is None:
        filters_enabled = [True] * len(filters_config)

    # Initialize filters
    filters_to_apply = []
    for filter_conf in filters_config:
        filter_class = filter_conf['class']
        params = filter_conf['params']
        title = filter_conf.get('title', '')
        filter_instance = filter_class(sequence, title=title, **params)
        filters_to_apply.append(filter_instance)

    # Get annotations per filter
    annotations_per_filter = []
    for filter in filters_to_apply:
        #filter.process() # We already processed the filters in the constructor
        annotations = filter.get_annotations()
        annotations_per_filter.append(annotations)

    # Collect annotations from enabled filters
    enabled_annotations = []
    for enabled, annotations in zip(filters_enabled, annotations_per_filter):
        if enabled:
            enabled_annotations.append(annotations)

    # Optimize the sequence
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    suggestions_for_each_codon = []
    current_scores_for_each_codon = []

    for i in range(0, len(sequence), 3):
        start_index = i
        end_index = i + 3
        suggestions_for_codon = []
        current_codon_score = 0
        for annotations in enabled_annotations:
            for annotation in annotations:
                if annotation['start'] == start_index and annotation['end'] == end_index:
                    suggestions_for_codon.append(annotation['suggestions'])
                    current_codon_score += annotation['score']
                    break
        suggestions_for_each_codon.append(suggestions_for_codon)
        current_scores_for_each_codon.append(current_codon_score)

    best_codons = []
    for i, suggestions_for_current_codon in enumerate(suggestions_for_each_codon):
        codon_scores = {}
        for suggestions in suggestions_for_current_codon:
            for suggestion, score in suggestions:
                codon_scores[suggestion] = codon_scores.get(suggestion, 0) + score
        if not codon_scores:
            best_codons.append(codons[i])
            continue
        best_codon = max(codon_scores, key=codon_scores.get)
        if codon_scores[best_codon] <= current_scores_for_each_codon[i]:
            best_codon = codons[i]
        best_codons.append(best_codon)

    optimized_sequence = ''.join(best_codons)
    return optimized_sequence
