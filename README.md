# SeqRET: Sequence Refinement for Enhanced Translation
*backronyms are still cool, right? we're still doing backronyms?*

SeqRET is a python package to optimize the codons of a DNA or RNA sequence for improved translation.


## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Filters](#filters)
- [Adding a Filter](#adding-a-filter)
- [License](#license)

## Installation

First, clone the directory anywhere on your machine using the following command:

```git clone https://github.com/pinneylab/SeqRET.git```

Then, install the package to your python or conda environment using:

```pip install -e SeqRET/```


## Usage

In your conda or python environment, you can now run the app using the command `seqret`.

This will start a web app that can be accessed at `localhost:8050` in your browser.

### Inputting a sequence:
Paste a DNA sequence into the text box and click "Submit".

### Editing the sequence:
Each filter will score every codon in your sequence, and show the results in the center display.

Clicking on a codon on a certain filter, you can see a suggestion for codons that improve the score for that filter while maintaining the same amino acid.

### Autorun:
Next to each filter is a check box. To automatically optimize the sequence, click the box next to each you filter you'd like to apply.

Then, hit the "Run Selected Filters" button to apply the filters. For each codon with an imperfect score, this will find a codon that satisfies all (or most) of the selected filters.

### Exporting the sequence:
Each change you make to the sequence will be reflected in the text box at the left of the page. You may copy this sequence as output.

*N.B.: Currently there is an error that the output sequence in the text box is not updated after hitting the 'Run Selected Filters' button. This can be circumvented by clicking on any codon and changing it (to the same codon) to trigger the update.*


## Filters

Currently, there are four filters:

### Banned Sequence Filter
This filter will give a score of 0 to any codon that is part of a banned sequence. The banned sequences are listed below, and can be altered in the `seqret/filters.py` file:
    
```python
banned_sequences = ['AGGAG', 
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
                    'GGGGGGGG']
```

### Banned Codon Filter
This filter will give a score of 0 to any codon that is part of a banned codon. The banned codons are listed below, and can be altered in the `seqret/filters.py` file:

```python

banned_codons = ['GTG', 'TTG', 'GAG', 'GGA']
```

### Frequency Filter
This filter will score each codon based on its relative frequency in E. Coli. The frequency data is stored in the `seqret/filters.py` file, and can be altered there.

## Adding a Filter

Users can define a new filter in `seqret/filters.py` by extending the SequenceFilter class and implementing the process() method.

The process() method takes no input, but may use the nucleotide sequence stored in `self.sequence` to calculate a score for each codon.

It must store its output in `self.annotations`, which breaks the sequence up into a list of regions (usually codons) with scores and suggestions for replacement codons.
`self.annotations` is a list of dicts, each of the form:

```python
{ 'start': start_index,
  'end': end_index,
  'score': score,
  'suggestions': [ (suggestion_1, score_1), (suggestion_2, score_2), ... ]}
```

The score_to_color() method can also be overridden to change the color scheme for the filter.

## License

This package is licensed under the MIT License.