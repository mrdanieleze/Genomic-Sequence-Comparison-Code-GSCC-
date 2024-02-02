# Genomic Sequence Comparison Code (GSCC)
# Author: Daniel Eze
# License: MIT

# Description:
# This script compares genomic sequences using pairwise alignment.

# Usage:
# - Provide a target genomic sequence (target_sequence) for comparison.
# - Prepare a list of known genomic sequences (known_sequences) to compare with.

# Results:
# The script returns a list of best alignments with simplified scores.

# Note: Genomic sequences can be either DNA or RNA. This code is primarily designed for comparing suspected sequences with known reference sequences.

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

# Constants
NORMALIZATION_FACTOR = 10 
# Used to normalize the alignment score to a range of 0-10, providing a more interpretable and standardized score.

def compare_genomic_sequences_list(target_sequence, known_sequences):
    """
    Compare a target genomic sequence with a list of known genomic sequences using pairwise alignment.

    Parameters:
    - target_sequence (str): The target genomic sequence for comparison.
    - known_sequences (list): List of known genomic sequences to compare with.

    Returns:
    List of best alignments with simplified scores.
    """
    # Basic error handling
    if not target_sequence or not known_sequences:
        raise ValueError("Target sequence and known sequences must not be empty.")

    aligner = PairwiseAligner()
    aligner.mode = 'global'
    
    best_alignments = []

    for known_sequence in known_sequences:
        alignment = aligner.align(Seq(target_sequence), Seq(known_sequence))[0]
        normalized_score = (alignment.score / max(len(target_sequence), len(known_sequence))) * NORMALIZATION_FACTOR
        best_alignment_info = (str(alignment[0]), str(alignment[1]), round(normalized_score, 1),
                               alignment.aligned[0], alignment.aligned[1])
        best_alignments.append(best_alignment_info)

    return best_alignments

# Example usage:
target_sequence = "TCTTGATCAT"  # Replace with your target sequence
known_sequences = ["TTATCCACA", "TTGTCCACA", "TTATCCATA"]  # Replace with your list of known sequences

results = compare_genomic_sequences_list(target_sequence, known_sequences)

# Assuming 'results' contains the alignment results
simplified_results = [(a, b, score) for a, b, score, _, _ in results]
print(simplified_results)
