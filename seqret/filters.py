from Bio.Data import CodonTable
import time
import re

### Helper dictionaries ###
codon_to_AA = CodonTable.unambiguous_dna_by_name['Standard'].forward_table
# Invert the forward_table to get a mapping from amino acids to codons
AA_to_codons = {}
for codon, aa in codon_to_AA.items():
    if aa not in AA_to_codons:
        AA_to_codons[aa] = []
    AA_to_codons[aa].append(codon)

### Filters ###
# Users can define a new filter by extending the SequenceFilter class and implementing the process() method.
# The score_to_color() method can also be overridden to change the color scheme for the filter.

class SequenceFilter():
    def __init__(self, sequence, title):
        self.sequence = sequence
        self.title = title
        self.AA_sequence = self.generate_AA_sequence(self.sequence) #persistent
        self.annotations = None

        self.process()

    def process(self):
        '''
        For a stored sequence, this updates the scores and suggestions based on the implemented filter.
        Input: 
            None, but uses self.sequence
            self.sequence: a string representing the nucleotide sequence
        Output:
            None, but updates self.annotations
            self.annotations:
                This breaks the sequence into regions, and gives the scores and suggestions for each region.
                A list of dicts. Each dict is of the form:
                    { 'start': start_index,
                      'end': end_index,
                      'score': score,
                      'suggestions': [ (suggestion_1, score_1), (suggestion_2, score_2), ... ]}
        '''
        raise NotImplementedError
    
    def score_to_color(self, score):
        '''
        Given a score, returns a color. This may be re-implemented for a new filter, if it's useful.
        Input:
            score: a float between 0 and 1
        Output:
            color: a string representing a color to be used by Dash
        '''
        if score > 0.75:
            return 'green'
        elif score > 0.5:
            return 'yellow'
        elif score > 0.25:
            return 'orange'
        else:
            return 'red'

    def get_sequence(self):
        return self.sequence

    def get_AA_sequence(self):
        return self.AA_sequence

    def get_title(self):
        return self.title
    
    def update_sequence(self, new_sequence, force=False):
        if force or self.sequence in ('', None):
            s = self.generate_AA_sequence(new_sequence)
            self.AA_sequence = s
        else:
            #confirm that the new sequence has same AA as old sequence:
            new_AA_sequence = self.generate_AA_sequence(new_sequence)
            if new_AA_sequence != self.AA_sequence:
                raise ValueError('new sequence must have same AA sequence as old sequence.\nold AA sequence: {}\nnew AA sequence: {}\n To override, set force=True'.format(self.AA_sequence, new_AA_sequence))
        
        #update the sequence
        self.sequence = new_sequence
    
    def generate_AA_sequence(self, s):
        #given sequence s,
        #returns the AA sequence
        #as a string
        AA_sequence = ''
        for i in range(0, len(s), 3):
            if s[i:i+3] not in codon_to_AA.keys():
                AA_sequence += 'X'
            else:
                AA_sequence += codon_to_AA[s[i:i+3]]
        return AA_sequence
    
    def get_annotations(self):
        '''
        Returns the annotations for the current sequence.
        '''
        return self.annotations
    
    def apply_suggestion(self, start, end, suggestion, force=False):
        '''
        Given a suggestion, applies it to the sequence and returns the new sequence.
        Input:
            start: (int) start index of the region to change
            end:   (int) end index of the region to change
            suggestion: (string) suggestion
            force: (bool) if True, will apply the suggestion even if it changes the AA sequence
        Output:
            None, but updates self.sequence
            self.sequence: (string) the new sequence
        '''
        #check that the suggestion is the right length:
        if len(suggestion) != end-start:
            raise ValueError('suggestion must be the same length as the region to change.\nstart: {}\nend: {}\nsuggestion: {}'.format(start, end, suggestion))
        
        #apply the suggestion
        new_sequence = self.sequence[:start] + suggestion + self.sequence[end:]

        #update the sequence
        self.update_sequence(new_sequence, force=force)

#make a codon frequency filter
class FrequencyFilter(SequenceFilter):
    def __init__(self, sequence, frequency_dict, title='E. Coli Frequency Filter'):
        '''
        Inputs:
            sequence: a string representing the nucleotide sequence
            frequency_dict: a dictionary of dictionaries, where the first key is the amino acid, and the second key is the codon. The value is the frequency of that codon for that amino acid.
            title: a string representing the title of the filter
        '''
        self.freq_dict = frequency_dict
        super().__init__(sequence, title)
    
    def process(self):
        '''
        For a stored sequence, this updates the scores and suggestions based on the implemented filter.
        Input: 
            None, but uses self.sequence
            self.sequence: a string representing the nucleotide sequence
        Output:
            None, but updates self.annotations
            self.annotations:
                This breaks the sequence into regions, and gives the scores and suggestions for each region.
                A list of dicts. Each dict is of the form:
                    { 'start': start_index,
                      'end': end_index,
                      'score': score,
                      'suggestions': [ (suggestion_1, score_1), (suggestion_2, score_2), ... ]}
        '''
        if self.sequence in ('', None):
            self.annotations = []
            return
        
        #split string into codons:
        codon_list = [self.sequence[i:i+3] for i in range(0, len(self.sequence), 3)]

        #get the frequency for each amino acid in the sequence
        frequency_list = []
        for i in range(len(self.AA_sequence)):
            AA = self.AA_sequence[i]
            if AA == 'X': #unknown AA, from unknown codon
                frequency_list += [0]
            else:
                frequency_list += [self.freq_dict[AA][codon_list[i]]]

        #for each codon, suggest the most frequent codon for that AA
        suggestion_list = []
        for AA in self.AA_sequence:
            if AA == 'X':
                suggestion_list.append([])
            else:
                #add tuple of (codon, frequency) for current AA:
                suggestion_list.append([(codon, self.freq_dict[AA][codon]) for codon in self.freq_dict[AA].keys()])

        #make the annotations
        self.annotations = []
        for i in range(len(codon_list)):
            self.annotations.append({
                'start': i*3,
                'end': i*3+3,
                'score': frequency_list[i],
                'suggestions': suggestion_list[i]
            })

    def score_to_color(self, score):
        '''
        Given a score, returns a color.
        Input:
            score: a float between 0 and 1
        Output:
            color: a string representing a color to be used by Dash
        '''
        if score > 0.75:
            return 'green'
        elif score > 0.5:
            return 'yellow'
        elif score > 0.25:
            return 'orange'
        else:
            return 'red'

class BannedCodonFilter(SequenceFilter):
    def __init__(self, sequence, banned_codons, title='Banned Codons'):
        self.banned_codons = banned_codons
        super().__init__(sequence, title)
    
    def process(self):
        '''
        For a stored sequence, this updates the scores and suggestions based on the implemented filter.
        Input: 
            None, but uses self.sequence
            self.sequence: a string representing the nucleotide sequence
        Output:
            None, but updates self.annotations
            self.annotations:
                This breaks the sequence into regions, and gives the scores and suggestions for each region.
                A list of dicts. Each dict is of the form:
                    { 'start': start_index,
                      'end': end_index,
                      'score': score,
                      'suggestions': [ (suggestion_1, score_1), (suggestion_2, score_2), ... ]}
        '''
        if self.sequence in ('', None):
            self.annotations = []
            return
        
        #split string into codons:
        codon_list = [self.sequence[i:i+3] for i in range(0, len(self.sequence), 3)]

        #for each codon, is it in the set of banned codons?
        banned_list = []
        suggestion_list = []
        for i, codon in enumerate(codon_list):
            if codon in self.banned_codons:
                banned_list.append(True)
            else:
                banned_list.append(False)
            #get alternate codons, and make sure they aren't banned
            AA = self.AA_sequence[i]
            if AA == 'X': #unknown AA, from unknown codon
                alternate_codons = []
            else:
                alternate_codons = [(codon, 1) for codon in AA_to_codons[AA] if codon not in self.banned_codons]
            suggestion_list.append(alternate_codons)

        #make the annotations
        self.annotations = []
        for i in range(len(codon_list)):
            score = 0 if banned_list[i] else 1
            self.annotations.append({
                'start': i*3,
                'end': i*3+3,
                'score': score,
                'suggestions': suggestion_list[i]
            })

    def score_to_color(self, score):
        '''
        Given a score, returns a color.
        Input:
            score: a float between 0 and 1
        Output:
            color: a string representing a color to be used by Dash
        '''
        if score == 1:
            return 'green'
        else:
            return 'red'

class BannedSequencesFilter(SequenceFilter):
    def __init__(self, sequence, banned_sequences, title='Banned Sequence'):
        self.banned_sequences = banned_sequences
        super().__init__(sequence, title)
        
    
    def process(self):
        '''
        For a stored sequence, this updates the scores and suggestions based on the implemented filter.
        Input: 
            None, but uses self.sequence
            self.sequence: a string representing the nucleotide sequence
        Output:
            None, but updates self.annotations
            self.annotations:
                This breaks the sequence into regions, and gives the scores and suggestions for each region.
                A list of dicts. Each dict is of the form:
                    { 'start': start_index,
                      'end': end_index,
                      'score': score,
                      'suggestions': [ (suggestion_1, score_1), (suggestion_2, score_2), ... ]}
        '''
        if self.sequence in ('', None):
            self.annotations = []
            return
        
        #scan over full sequence for banned subsequences. Make note of nucleotide indices that are part of a banned subsequence.
        crap_indices = self.get_bad_indices(self.sequence)
        nucleotide_scores = [0 if crap_indices[i] else 1 for i in range(len(crap_indices))]

        #split into codons:
        codon_list = [self.sequence[i:i+3] for i in range(0, len(self.sequence), 3)]
        codon_scores = [0 if any(crap_indices[i:i+3]) else 1 for i in range(0, len(self.sequence), 3)]

        #for each bad codon, suggest a new codon that is not part of a banned subsequence
        codon_suggestions = [] #this will be a list of lists, where each sub-list contains the suggestions of new codons for a given codon
        for i in range(len(codon_list)):
            #create sublist for each codon
            current_codon_suggestions = []
            #check all alternate codons. If they are not part of a banned subsequence, add them to the list of suggestions
            #if codon_scores[i] == 0:
            AA = self.AA_sequence[i]
            if AA == 'X': #unknown AA, from unknown codon
                alternate_codons = []
            else:
                alternate_codons = [(codon, 0) for codon in AA_to_codons[AA]]
            #will the alternate codons be part of a banned subsequence?
            #for each alternate codon:
            for alt_codon in alternate_codons:
                #make temp sequence with new codon inserted:
                temp_sequence = self.sequence[:i*3] + alt_codon[0] + self.sequence[i*3+3:]
                #check if the new sequence has a banned subsequence:
                crap_indices = self.get_bad_indices(temp_sequence)[i*3:i*3+3]
                #if it does, remove the alternate codon from the list of suggestions
                if not any(crap_indices):
                    current_codon_suggestions.append((alt_codon[0], 1))
            codon_suggestions.append(current_codon_suggestions)

        #for each nucleotide that is part of a banned subsequence, give its containing nucleotide a score of 0. Suggest alternate nucleotides.
        #otherwise, give it a score of 1.
        #make the annotations
        self.annotations = []
        for i in range(len(codon_list)):
            self.annotations.append({
                'start': i*3,
                'end': i*3+3,
                'score': codon_scores[i],
                'suggestions': codon_suggestions[i]
            })
        
    def get_bad_indices(self, seq):
        # crap_indices = [False]*len(seq)
        # for banned_sequence in self.banned_sequences:
        #     for i in range(len(seq)-len(banned_sequence)):
        #         if seq[i:i+len(banned_sequence)] == banned_sequence:
        #             crap_indices[i:i+len(banned_sequence)] = [True]*len(banned_sequence)
        # return crap_indices
        pattern = '|'.join(map(re.escape, self.banned_sequences))
        matches = re.finditer(pattern, seq)
        crap_indices = [False] * len(seq)
        for match in matches:
            start, end = match.span()
            for i in range(start, end):
                crap_indices[i] = True
        return crap_indices

    def score_to_color(self, score):
        '''
        Given a score, returns a color.
        Input:
            score: a float between 0 and 1
        Output:
            color: a string representing a color to be used by Dash
        '''
        if score == 1:
            return 'green'
        else:
            return 'red'

# Define the configurations for the filters
def get_filters():
    banned_codons = ['CTG', 'GTG', 'TTG', 'GAG', 'GGA']
    banned_sequences = ['AGGAG', 'GAGGT', 'AAGGA', 'TAAGG', 'GGAGG', 'GGGGG',
                        'TAAGGA', 'AAGGAG', 'AGGAGG', 'GGAGGT', 'GGGGGG',
                        'TAAGGAG', 'AAGGAGG', 'AGGAGGT', 'GGGGGGG',
                        'TAAGGAGG', 'AAGGAGGT', 'GGGGGGGG']
    E_coli_AA_freq_dict = {
        'F': {'TTT': 0.58, 'TTC': 0.42}, 'L': {'TTA': 0.14, 'TTG': 0.13, 'CTT': 0.12,
        'CTC': 0.1, 'CTA': 0.04, 'CTG': 0.47}, 'Y': {'TAT': 0.59, 'TAC': 0.41},
        '*': {'TAA': 0.61, 'TAG': 0.09, 'TGA': 0.3}, 'H': {'CAT': 0.57, 'CAC': 0.43},
        'Q': {'CAA': 0.34, 'CAG': 0.66}, 'I': {'ATT': 0.49, 'ATC': 0.39, 'ATA': 0.11},
        'M': {'ATG': 1.0}, 'N': {'AAT': 0.49, 'AAC': 0.51}, 'K': {'AAA': 0.74, 'AAG': 0.26},
        'V': {'GTT': 0.28, 'GTC': 0.2, 'GTA': 0.17, 'GTG': 0.35}, 'D': {'GAT': 0.63, 'GAC': 0.37},
        'E': {'GAA': 0.68, 'GAG': 0.32}, 'S': {'TCT': 0.17, 'TCC': 0.15, 'TCA': 0.14,
        'TCG': 0.14, 'AGT': 0.16, 'AGC': 0.25}, 'C': {'TGT': 0.46, 'TGC': 0.54},
        'W': {'TGG': 1.0}, 'P': {'CCT': 0.18, 'CCC': 0.13, 'CCA': 0.2, 'CCG': 0.49},
        'R': {'CGT': 0.36, 'CGC': 0.36, 'CGA': 0.07, 'CGG': 0.11, 'AGA': 0.07, 'AGG': 0.04},
        'T': {'ACT': 0.19, 'ACC': 0.4, 'ACA': 0.17, 'ACG': 0.25},
        'A': {'GCT': 0.18, 'GCC': 0.26, 'GCA': 0.23, 'GCG': 0.33},
        'G': {'GGT': 0.35, 'GGC': 0.37, 'GGA': 0.13, 'GGG': 0.15}
    }

    return [
        {'class': BannedSequencesFilter, 'params': {'banned_sequences': banned_sequences}, 'title': 'Banned Sequence'},
        {'class': BannedCodonFilter, 'params': {'banned_codons': banned_codons}, 'title': 'Banned Codons'},
        {'class': FrequencyFilter, 'params': {'frequency_dict': E_coli_AA_freq_dict}, 'title': 'E. Coli Frequency Filter'}
    ]

### Define Filters to Use ###
# Eventually, this should contain all the filters available to the app, and at runtime the user can check boxes interactively to enable/disable filters.
# For now, we specify them here. The "Secondary Structure" filter is disabled because it is slow for sequences longer than ~300 nucleotides. 

# banned_codons = ['CTG', 'GTG', 'TTG', 'GAG', 'GGA']
# banned_sequences = ['AGGAG', 
#                     'GAGGT', 
#                     'AAGGA', 
#                     'TAAGG', 
#                     'GGAGG', 
#                     'GGGGG', 
#                     'TAAGGA', 
#                     'AAGGAG', 
#                     'AGGAGG', 
#                     'GGAGGT', 
#                     'GGGGGG', 
#                     'TAAGGAG', 
#                     'AAGGAGG', 
#                     'AGGAGGT', 
#                     'GGGGGGG', 
#                     'TAAGGAGG', 
#                     'AAGGAGGT', 
#                     'GGGGGGGG']
# E_coli_AA_freq_dict = {'F': {'TTT': 0.58, 'TTC': 0.42},                            'L': {'TTA': 0.14, 'TTG': 0.13, 'CTT': 0.12, 'CTC': 0.1, 'CTA': 0.04, 'CTG': 0.47}, 
#                 'Y': {'TAT': 0.59, 'TAC': 0.41},                            '*': {'TAA': 0.61, 'TAG': 0.09, 'TGA': 0.3}, 
#                 'H': {'CAT': 0.57, 'CAC': 0.43},                            'Q': {'CAA': 0.34, 'CAG': 0.66}, 
#                 'I': {'ATT': 0.49, 'ATC': 0.39, 'ATA': 0.11},               'M': {'ATG': 1.0}, 
#                 'N': {'AAT': 0.49, 'AAC': 0.51},                            'K': {'AAA': 0.74, 'AAG': 0.26}, 
#                 'V': {'GTT': 0.28, 'GTC': 0.2, 'GTA': 0.17, 'GTG': 0.35},   'D': {'GAT': 0.63, 'GAC': 0.37}, 
#                 'E': {'GAA': 0.68, 'GAG': 0.32},                            'S': {'TCT': 0.17, 'TCC': 0.15, 'TCA': 0.14, 'TCG': 0.14, 'AGT': 0.16, 'AGC': 0.25}, 
#                 'C': {'TGT': 0.46, 'TGC': 0.54},                            'W': {'TGG': 1.0}, 
#                 'P': {'CCT': 0.18, 'CCC': 0.13, 'CCA': 0.2, 'CCG': 0.49},   'R': {'CGT': 0.36, 'CGC': 0.36, 'CGA': 0.07, 'CGG': 0.11, 'AGA': 0.07, 'AGG': 0.04}, 
#                 'T': {'ACT': 0.19, 'ACC': 0.4, 'ACA': 0.17, 'ACG': 0.25},   'A': {'GCT': 0.18, 'GCC': 0.26, 'GCA': 0.23, 'GCG': 0.33}, 
#                 'G': {'GGT': 0.35, 'GGC': 0.37, 'GGA': 0.13, 'GGG': 0.15}
#                 }

# filters_to_apply = [BannedSequencesFilter('', banned_sequences), BannedCodonFilter('', banned_codons), FrequencyFilter('', E_coli_AA_freq_dict)]
