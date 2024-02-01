"""
Genomic Sequence Comparison Code (GSCC)
Author: Daniel Eze
License: MIT

Description:
This script compares DNA sequences using pairwise alignment.

Usage:
- Provide a DNA sequence (sequence_a) for comparison.
- Prepare a list of DNA sequences (sequence_b) to compare with.

Results:
The script returns a list of best alignments with simplified scores.

Note: DNA sequences are to be compared with only DNA sequences (same goes for RNA).
This code is primarily designed for comparing suspected sequences with known reference sequences.
"""

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

def compare_dna_sequences_list(sequences_a, sequences_b):
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    
    best_alignments = []

    for sequence_a in sequences_a:
        for sequence_b in sequences_b:
            alignment = aligner.align(Seq(sequence_a), Seq(sequence_b))[0]
            normalized_score = (alignment.score / max(len(sequence_a), len(sequence_b))) * 10
            best_alignment_info = (str(alignment[0]), str(alignment[1]), round(normalized_score, 1),
                                   alignment.aligned[0], alignment.aligned[1])
            best_alignments.append(best_alignment_info)

    return best_alignments

# Provided sequences
sequence_a = ""                     # Your target sequence here
sequence_b = ["", "", ""]           # List of Known sequences here

results = compare_dna_sequences_list([sequence_a], sequence_b)

# Assuming 'results' contains the alignment results

simplified_results = [(a, b, score) for a, b, score, _, _ in results]
print(simplified_results)
