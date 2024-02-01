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
The script returns a list of best alignments with normalized scores.

Note: DNA sequences are to be compared with only DNA sequences (same goes for RNA).
This code is primarily designed for comparing suspected sequences with known reference sequences.
"""

from Bio import pairwise2
from Bio.Seq import Seq

def compare_dna_sequences_list(sequences_a, sequences_b):
    best_alignments = []

    for sequence_a in sequences_a:
        for sequence_b in sequences_b:
            alignments = pairwise2.align.globalxx(Seq(sequence_a), Seq(sequence_b), one_alignment_only=True)
            best_alignment = alignments[0] if alignments else None
            if best_alignment:
                normalized_score = (best_alignment[2] / max(len(sequence_a), len(sequence_b))) * 10  # Normalize to a scale from 0 to 10
                best_alignment = (*best_alignment[:2], normalized_score, *best_alignment[3:])
            best_alignments.append(best_alignment)

    return best_alignments

# Provided sequences
sequence_a = ""                         
sequence_b = ["", "", ""]            # Sequences are to be compared with only sequences (same goes for RNA)

results = compare_dna_sequences_list([sequence_a], sequence_b)

print(results)
