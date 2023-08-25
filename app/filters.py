from Bio.Data import CodonTable
import RNA

codon_to_AA = CodonTable.unambiguous_dna_by_name['Standard'].forward_table
AA_freq_dict = {'F': {'TTT': 0.58, 'TTC': 0.42}, 'L': {'TTA': 0.14, 'TTG': 0.13, 'CTT': 0.12, 'CTC': 0.1, 'CTA': 0.04, 'CTG': 0.47}, 'Y': {'TAT': 0.59, 'TAC': 0.41}, '*': {'TAA': 0.61, 'TAG': 0.09, 'TGA': 0.3}, 'H': {'CAT': 0.57, 'CAC': 0.43}, 'Q': {'CAA': 0.34, 'CAG': 0.66}, 'I': {'ATT': 0.49, 'ATC': 0.39, 'ATA': 0.11}, 'M': {'ATG': 1.0}, 'N': {'AAT': 0.49, 'AAC': 0.51}, 'K': {'AAA': 0.74, 'AAG': 0.26}, 'V': {'GTT': 0.28, 'GTC': 0.2, 'GTA': 0.17, 'GTG': 0.35}, 'D': {'GAT': 0.63, 'GAC': 0.37}, 'E': {'GAA': 0.68, 'GAG': 0.32}, 'S': {'TCT': 0.17, 'TCC': 0.15, 'TCA': 0.14, 'TCG': 0.14, 'AGT': 0.16, 'AGC': 0.25}, 'C': {'TGT': 0.46, 'TGC': 0.54}, 'W': {'TGG': 1.0}, 'P': {'CCT': 0.18, 'CCC': 0.13, 'CCA': 0.2, 'CCG': 0.49}, 'R': {'CGT': 0.36, 'CGC': 0.36, 'CGA': 0.07, 'CGG': 0.11, 'AGA': 0.07, 'AGG': 0.04}, 'T': {'ACT': 0.19, 'ACC': 0.4, 'ACA': 0.17, 'ACG': 0.25}, 'A': {'GCT': 0.18, 'GCC': 0.26, 'GCA': 0.23, 'GCG': 0.33}, 'G': {'GGT': 0.35, 'GGC': 0.37, 'GGA': 0.13, 'GGG': 0.15}}

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
        Given a score, returns a color.
        Input:
            score: a float between 0 and 1
        Output:
            color: a string representing a color to be used by Dash
        '''
        raise NotImplementedError

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
    def __init__(self, sequence, title='Frequency Filter'):
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
                frequency_list += [AA_freq_dict[AA][codon_list[i]]]

        #for each codon, suggest the most frequent codon for that AA
        suggestion_list = []
        for AA in self.AA_sequence:
            if AA == 'X':
                suggestion_list.append([])
            else:
                #add tuple of (codon, frequency) for current AA:
                suggestion_list.append([(codon, AA_freq_dict[AA][codon]) for codon in AA_freq_dict[AA].keys()])

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
    def __init__(self, sequence, title='Banned Codons'):
        super().__init__(sequence, title)
        self.banned_codons = ['GTG', 'TTG', 'GAG', 'GGA']
    
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
                alternate_codons = [(codon, 1) for codon in AA_freq_dict[AA].keys() if codon not in self.banned_codons]
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
    def __init__(self, sequence, title='Banned Sequence'):
        self.banned_sequences = ['AGGAG', 
                                 'GAGGT', 
                                 'AAGGA', 
                                 'TAAGG', 
                                 'GGAGG', 
                                 'GGGGG', 
                                 'TAAGGA', 
                                 'AAGGAG', 
                                 'AGGAGG', 
                                 'GGAGGT', 
                                 'GGGGGG', 
                                 'TAAGGAG', 
                                 'AAGGAGG', 
                                 'AGGAGGT', 
                                 'GGGGGGG', 
                                 'TAAGGAGG', 
                                 'AAGGAGGT', 
                                 'GGGGGGGG'] #double check later
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
                alternate_codons = [(codon, 0) for codon in AA_freq_dict[AA].keys()]
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
        #for each nucleotide that is not part of a banned subsequence, give its containing nucleotide a score of 1. Suggest no alternate nucleotides.
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
        crap_indices = [False]*len(seq)
        for banned_sequence in self.banned_sequences:
            for i in range(len(seq)-len(banned_sequence)):
                if seq[i:i+len(banned_sequence)] == banned_sequence:
                    crap_indices[i:i+len(banned_sequence)] = [True]*len(banned_sequence)
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

class SecondaryStructureFilter(SequenceFilter):
    def __init__(self, sequence, title='RNA Secondary Structure'):
        self.secondary_structure = None
        self.temperature = 37.0
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
        def paired_positions(sequence, structure):
            """
            Given a sequence and its corresponding secondary structure in dot-bracket notation,
            return a list where each position has the paired nucleotide for the corresponding nucleotide
            in the sequence based on the secondary structure. If the nucleotide is not base paired, 
            that position in the list will contain None.
            """
            stack = []
            paired_list = [None] * len(sequence)

            for i, (nt, symbol) in enumerate(zip(sequence, structure)):
                if symbol == '(':
                    stack.append(i)
                elif symbol == ')':
                    start = stack.pop()
                    paired_list[start] = nt
                    paired_list[i] = sequence[start]

            return paired_list
        
        if self.sequence in ('', None):
            self.annotations = []
            return
        
        #first, generate the secondary structure using viennaRNA
        md = RNA.md()
        # change temperature to room temperature
        md.temperature = self.temperature 
        
        #convert DNA to RNA
        rna_sequence = self.sequence.replace('T', 'U')
        # compute minimum free energy (MFE) and corresponding structure
        fc = RNA.fold_compound(rna_sequence, md)
        (secondary_sequence, mfe) = fc.mfe()

        self.secondary_structure = secondary_sequence
        #split into codons:
        codon_list = [self.sequence[i:i+3] for i in range(0, len(self.sequence), 3)]
        ss_list = [secondary_sequence[i:i+3] for i in range(0, len(secondary_sequence), 3)]

        #get paired sequence:
        paired_list = paired_positions(self.sequence, secondary_sequence)
        paired_codon_list = [paired_list[i:i+3] for i in range(0, len(paired_list), 3)]

        #score each codon. If the secondary structure is not bound ('.'), give it a score of 1. Otherwise, give it a score of 0.
        codon_scores = []
        for i in range(len(codon_list)):
            #count the number of '.'s in ss_list[i]
            score = ss_list[i].count('.')/3
            codon_scores.append(score)
        
        # #for nucleotides that are base pairing, we want to swap them with ones that we know won't bind with their partner.
        # non_binding_nucleotides = {
        #     'A': ['U', 'C'],
        #     'U': ['A', 'G'],
        #     'G': ['U', 'C'],
        #     'C': ['A', 'G']
        # }
        #if we want to change a nucleotide that already isn't pairnig, we need to ensure we won't cause it to pair.
        #this dict holds acceptable alternatives:
        acceptable_alternatives = {
            'A': ['A'],
            'U': ['U', 'C'],
            'G': ['A', 'G'],
            'C': ['C']
        }
        non_pairing_partners = {
            'A': ['A', 'C', 'G'],
            'U': ['U', 'C'],
            'G': ['G', 'A'],
            'C': ['C', 'A', 'U']
        }

        codon_suggestions = []
        for i in range(len(codon_list)):
            # #if no pairing, make no suggestions:
            # if codon_scores[i] == 1:
            #     codon_suggestions.append([])
            # #if pairing, suggest alternate codons:
            # else:
            #get alternate codons, and make sure they aren't banned
            AA = self.AA_sequence[i]
            if AA == 'X': #unknown AA, from unknown codon
                alternate_codons = []
            else:
                alternate_codons = [(codon, 0) for codon in AA_freq_dict[AA].keys()]

            #for each alternate codon, count up how many bound nucleotides it resolves!
            alternate_codon_scores = []
            for alt_codon in alternate_codons:
                current_alt_codon_score = 0
                for j in range(3):
                    #current_ss = ss_list[i][j]
                    current_nucleotide = codon_list[i][j]
                    current_paired_nucleotide = paired_codon_list[i][j]
                    alternate_nucleotide = alt_codon[0][j]

                    #for nucleotides that aren't bound, suggest a reasonable alternative.
                    if current_paired_nucleotide is None:
                        if alternate_nucleotide.replace('T', 'U') in acceptable_alternatives[current_nucleotide.replace('T', 'U')]:
                            current_alt_codon_score += 1
                    else:
                        #check if the alternate nucleotide proposed will not bind with its partner
                        if alternate_nucleotide.replace('T', 'U') in non_pairing_partners[current_paired_nucleotide.replace('T', 'U')]:
                            current_alt_codon_score += 1
                current_alt_codon_score /= 3
                alternate_codon_scores.append(current_alt_codon_score)
            
            #update scores
            for j in range(len(alternate_codons)):
                alternate_codons[j] = (alternate_codons[j][0], alternate_codon_scores[j])
            
            #remove suggestions with score lower than current codon
            alternate_codons = [alt_codon for alt_codon in alternate_codons if alt_codon[1] >= codon_scores[i]]
            
            #add the alternate codons to the list of suggestions
            codon_suggestions.append(alternate_codons)

        self.annotations = []
        for i in range(len(codon_list)):
            self.annotations.append({
                'start': i*3,
                'end': i*3+3,
                'score': codon_scores[i],
                'suggestions': codon_suggestions[i]
            })

    def score_to_color(self, score):
        '''
        Given a score, returns a color.
        Input:
            score: a float between 0 and 1
        Output:
            color: a string representing a color to be used by Dash
        '''
        if score > 0.67:
            return 'green'
        elif score > 0.34:
            return 'yellow'
        elif score > 0.01:
            return 'orange'
        else:
            return 'red'

    def get_secondary_structure(self):
        return self.secondary_structure


#filters_to_apply = [BannedSequencesFilter(''), FrequencyFilter(''), SecondaryStructureFilter('')]
filters_to_apply = [BannedSequencesFilter(''), BannedCodonFilter(''), FrequencyFilter(''), SecondaryStructureFilter('')]
#filters_to_apply = [SecondaryStructureFilter('')]
'''

No_No_Codons	
GTG	
TTG	
GAG	
GGA	

No_No_Sequences
AGGAG
GAGGT
AAGGA
TAAGG
GGAGG
GGGGG
TAAGGA
AAGGAG
AGGAGG
GGAGGT
GGGGGG
TAAGGAG
AAGGAGG
AGGAGGT
GGGGGGG
TAAGGAGG
AAGGAGGT
GGGGGGGG

CGCTGGCGGCGTGCTAATACATGCCCGTCGAGCNGGATCGATGGGAGCTTGCTCCCTGCAGATCAGCGGCGGACGGGTGAGTAACACGTGGGTAACCTGCCTGTAAGACTGGGATAACTCCGGGAAACCGGGGCTAATACCGGATAATTTAGTTCCTCGCATGAGGAACTGTTGAAAGGTGGCTTC
CGCTGGCGGCGTGCTAATACATGCAGGAGGTGCNGGATCGATGGGAGCTTGCTCCCTGCAGATCAGCGGCGGACGGGTGAGTAACACGTGGGTAACCTGCCTGTAAGACTGGGATAACTCCGGGAAACCGGGGCTAATACCGGATAATTTAGTTCCTCGCATGAGGAACTGTTGAAAGGTGGCTTC

'''