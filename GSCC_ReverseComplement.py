# Genomic Sequence Comparison Code (GSCC) - Reverse Complement
# Author: Daniel Eze
# License: MIT

# Description:
# This script compares genomic sequences using pairwise alignment, with a focus on reverse complements.

# Usage:
# - Provide a target genomic sequence (target_sequence) for comparison.
# - Prepare a list of known genomic sequences (known_sequences) to compare with.

# Results:
# The script returns both normal and reverse complement alignments, each with a list of best alignments and simplified scores.
# The positions of the normal alignments for the reverse complement alignments are the same, indicating a reverse complement relationship.

# Note: 
# Genomic sequences can be either DNA or RNA. This code is designed to compare suspected sequences with known reference sequences.
# Reverse complements of sequences are considered in the results, comparing both normal and reverse complement alignments and displaying both sets of results.

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

# Constants
NORMALIZATION_FACTOR = 10 
# Used to normalize the alignment score to a range of 0-10, providing a more interpretable and standardized score.

def compare_genomic_sequences_list(target_sequence, known_sequences):
    """
    Compare a target genomic sequence with a list of known genomic sequences using pairwise alignment, giving emphasis to reverse complements.

    Parameters:
    - target_sequence (str): The target genomic sequence for comparison.
    - known_sequences (list): List of known genomic sequences to compare with.

    Returns:
    Tuple containing:
    - List of best alignments with simplified scores for normal alignment.
    - List of best alignments with simplified scores for reverse complement alignment.
    """
    # Basic error handling
    if not target_sequence or not known_sequences:
        raise ValueError("Target sequence and known sequences must not be empty.")

    aligner = PairwiseAligner()
    aligner.mode = 'global'
    
    normal_alignments = []
    reverse_complement_alignments = []

    for known_sequence in known_sequences:
        # Compare with normal alignment
        alignment_normal = aligner.align(Seq(target_sequence), Seq(known_sequence))[0]
        normal_alignments.append((str(alignment_normal[0]), str(alignment_normal[1]), round(normalization_score(alignment_normal.score, target_sequence, known_sequence), 1),
                                  alignment_normal.aligned[0], alignment_normal.aligned[1]))

        # Compare with reverse complement alignment
        alignment_reverse_complement = aligner.align(Seq(target_sequence), Seq(known_sequence).reverse_complement())[0]
        reverse_complement_alignments.append((str(alignment_reverse_complement[0]), str(alignment_reverse_complement[1]), round(normalization_score(alignment_reverse_complement.score, target_sequence, known_sequence), 1),
                                              alignment_reverse_complement.aligned[0], alignment_reverse_complement.aligned[1]))

    return normal_alignments, reverse_complement_alignments

def normalization_score(score, sequence_a, sequence_b):
    """
    Calculate and normalize the alignment score.

    Parameters:
    - score (float): The alignment score.
    - sequence_a (str): The first sequence in the alignment.
    - sequence_b (str): The second sequence in the alignment.

    Returns:
    Normalized alignment score.
    """
    return (score / max(len(sequence_a), len(sequence_b))) * NORMALIZATION_FACTOR

# Example usage:
target_sequence = "TCTTGATCAT"  # Replace with your target sequence
known_sequences = ["TTATCCACA", "TTGTCCACA", "TTATCCATA"]  # Replace with your list of known sequences

normal_results, reverse_complement_results = compare_genomic_sequences_list(target_sequence, known_sequences)

# Assuming 'normal_results' and 'reverse_complement_results' contain the alignment results
print("Normal Alignments:")
simplified_normal_results = [(a, b, score) for a, b, score, _, _ in normal_results]
print(simplified_normal_results)

print("\nReverse Complement Alignments:")
simplified_reverse_complement_results = [(a, b, score) for a, b, score, _, _ in reverse_complement_results]
print(simplified_reverse_complement_results)