# DNA Motif Search
A Python-based bioinformatics tool for identifying motifs in DNA sequences using configurable 
sequence patterns, gaps, and deviation penalties.

## Features
* Reads DNA sequences from FASTA and plain-text files
* Searches DNA sequences for specified motifs
* Supports alternative nucleotides at specific positions
* Supports configurable gaps between motif elements
* Allows a maximum deviation/penalty threshold
* Reports the position, sequence, and penalty of each detected motif

## How It Works
The program takes three inputs:
1. A DNA sequence file
2. A motif/signal description file
3. A maximum allowed deviation

It then searches the DNA sequences for regions that match the specified motif and reports 
any matches within the allowed deviation.

## Usage
bash
python3 motif_search.py <sequence_file> <signal_description_file> <deviation>

Example:
bash
python3 motif_search.py sequences.fasta signal_description.txt 2


## Output
For each detected motif, the program reports:

* Start and end position
* Matching DNA sequence
* Penalty/deviation score

## Technologies
* Python 3
* Regular expressions
* FASTA sequence handling
* Bioinformatics sequence analysis
