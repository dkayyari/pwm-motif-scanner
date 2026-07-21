#!/usr/bin/env python3
"""
Position Weight Matrix (PWM) Motif Scanner for argR binding sites.

Part 1: Compute F(b,j), F'(b,j), and W(b,j) from the counts matrix.
         Output -> Part1_Matrices_Output.txt

Part 2: Scan upstream sequences with PWM and report top 30 gene IDs.
         Output -> Part2_Top30_Binding_Sites.txt

Usage
-----
Run with no arguments to use the default files in ../data/:

    python pwm_scan.py

Or point it at different input files / output directory:

    python pwm_scan.py --counts my_counts.txt --sequences my_seqs.txt -o my_results/
"""

import argparse
import math
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
DEFAULT_COUNTS = os.path.join(REPO_ROOT, "data", "argR-counts-matrix.txt")
DEFAULT_SEQUENCES = os.path.join(REPO_ROOT, "data", "E_coli_K12_MG1655.400_50")
DEFAULT_OUTPUT_DIR = os.path.join(REPO_ROOT, "results")


# ============================================================================
#  HELPER FUNCTIONS
# ============================================================================

def read_counts_matrix(filepath):
    """
    Read the counts matrix from the provided text file.
    Expected format per line: base | count1 count2 ... countN
    """
    counts = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('|')
            if len(parts) != 2:
                continue
            base = parts[0].strip().lower()
            count_values = list(map(int, parts[1].split()))
            counts[base] = count_values
    return counts


def compute_frequency_matrix(counts, total_sites):
    """Compute F(b, j) = count(b, j) / total_sites."""
    freq = {}
    for base in counts:
        freq[base] = [c / total_sites for c in counts[base]]
    return freq


def compute_adjusted_frequency_matrix(counts, total_sites_adjusted):
    """Compute F'(b, j) = (count(b, j) + 1) / total_sites_adjusted."""
    freq_prime = {}
    for base in counts:
        freq_prime[base] = [(c + 1) / total_sites_adjusted for c in counts[base]]
    return freq_prime


def compute_weight_matrix(freq_prime, background=0.25):
    """Compute W(b, j) = log2[ F'(b, j) / background ]."""
    weight = {}
    for base in freq_prime:
        weight[base] = [math.log2(f / background) for f in freq_prime[base]]
    return weight


def format_matrix(matrix, label, bases=('a', 'c', 'g', 't')):
    """Format a matrix as a string for display and file output."""
    lines = []
    lines.append(f"{'='*100}")
    lines.append(f"  {label}")
    lines.append(f"{'='*100}")
    n_positions = len(matrix[bases[0]])
    header = "Base |  " + "  ".join(f"{j+1:>7}" for j in range(n_positions))
    lines.append(header)
    lines.append("-" * len(header))
    for base in bases:
        row = f"  {base.upper()}  |  " + "  ".join(f"{v:7.4f}" for v in matrix[base])
        lines.append(row)
    lines.append("")
    return "\n".join(lines)


def format_counts_matrix(counts, bases=('a', 'c', 'g', 't')):
    """Format the counts matrix as a string."""
    lines = []
    lines.append(f"{'='*100}")
    lines.append("  Counts Matrix (original, from 27 binding sites)")
    lines.append(f"{'='*100}")
    n_positions = len(counts[bases[0]])
    header = "Base |  " + "  ".join(f"{j+1:>4}" for j in range(n_positions))
    lines.append(header)
    lines.append("-" * len(header))
    for base in bases:
        row = f"  {base.upper()}  |  " + "  ".join(f"{v:4d}" for v in counts[base])
        lines.append(row)
    lines.append("")
    return "\n".join(lines)


def read_fasta_sequences(filepath):
    """
    Read upstream regulatory sequences. Supports two formats:
    1. Custom format: gene_id \\ sequence (one per line)
    2. Standard FASTA: >header lines followed by sequence lines
    """
    sequences = []
    with open(filepath, 'r') as f:
        first_line = f.readline().strip()
        f.seek(0)

        if first_line.startswith('>'):
            gene_id = None
            seq_parts = []
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if gene_id is not None:
                        sequences.append((gene_id, ''.join(seq_parts).lower()))
                    gene_id = line[1:].split()[0]
                    seq_parts = []
                else:
                    seq_parts.append(line)
            if gene_id is not None:
                sequences.append((gene_id, ''.join(seq_parts).lower()))
        else:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if '\\' in line:
                    parts = line.split('\\', 1)
                    gene_id = parts[0].strip()
                    sequence = parts[1].strip().lower()
                    sequences.append((gene_id, sequence))
    return sequences


def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    complement = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}
    return ''.join(complement.get(base, 'n') for base in reversed(sequence))


def score_sequence(subsequence, weight_matrix):
    """Score a subsequence using the PWM."""
    score = 0.0
    for j, base in enumerate(subsequence):
        if base in weight_matrix:
            score += weight_matrix[base][j]
        else:
            score += 0.0
    return score


def scan_sequence(sequence, weight_matrix, motif_length):
    """
    Slide the PWM across both forward and reverse complement strands.
    Returns (best_score, best_position, best_subsequence, strand).
    """
    best_score = float('-inf')
    best_pos = 0
    best_subseq = ""
    best_strand = "+"

    # Scan forward strand
    for i in range(len(sequence) - motif_length + 1):
        subseq = sequence[i:i + motif_length]
        sc = score_sequence(subseq, weight_matrix)
        if sc > best_score:
            best_score = sc
            best_pos = i
            best_subseq = subseq
            best_strand = "+"

    # Scan reverse complement strand
    rev_seq = reverse_complement(sequence)
    for i in range(len(rev_seq) - motif_length + 1):
        subseq = rev_seq[i:i + motif_length]
        sc = score_sequence(subseq, weight_matrix)
        if sc > best_score:
            best_score = sc
            best_pos = i
            best_subseq = subseq
            best_strand = "-"

    return best_score, best_pos, best_subseq, best_strand


# ============================================================================
#  MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Build a PWM from argR binding sites and scan upstream regions for matches."
    )
    parser.add_argument(
        "--counts", default=DEFAULT_COUNTS,
        help=f"Counts matrix file (default: {DEFAULT_COUNTS})",
    )
    parser.add_argument(
        "--sequences", default=DEFAULT_SEQUENCES,
        help=f"Upstream sequences file (default: {DEFAULT_SEQUENCES})",
    )
    parser.add_argument(
        "-o", "--output-dir", default=DEFAULT_OUTPUT_DIR,
        help=f"Directory to write output files to (default: {DEFAULT_OUTPUT_DIR})",
    )
    args = parser.parse_args()

    # ---- File paths (input) ----
    counts_file = args.counts
    sequences_file = args.sequences

    # ---- Output file paths ----
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    part1_output = os.path.join(output_dir, "Part1_Matrices_Output.txt")
    part2_output = os.path.join(output_dir, "Part2_Top30_Binding_Sites.txt")

    bases = ('a', 'c', 'g', 't')

    # ========================================================================
    #  PART 1: Compute F(b,j), F'(b,j), and W(b,j)
    # ========================================================================

    # Step 1: Read counts matrix
    counts = read_counts_matrix(counts_file)
    motif_length = len(counts['a'])
    total_sites = 27  # Given: 27 binding sites

    # Step 2: Compute frequency matrix F(b, j) = count / 27
    freq = compute_frequency_matrix(counts, total_sites)

    # Step 3: Compute adjusted frequency matrix F'(b, j) = (count + 1) / 31
    total_adjusted = total_sites + 4  # 27 + 4 = 31
    freq_prime = compute_adjusted_frequency_matrix(counts, total_adjusted)

    # Step 4: Compute weight matrix W(b, j) = log2[ F'(b,j) / 0.25 ]
    weight_matrix = compute_weight_matrix(freq_prime, background=0.25)

    # --- Build Part 1 output ---
    part1_text = []
    part1_text.append("PART 1 — PWM CONSTRUCTION")
    part1_text.append("Computation of Frequency Matrix, Adjusted Frequency Matrix, and Weight Matrix")
    part1_text.append(f"Motif length: {motif_length} bases")
    part1_text.append(f"Total binding sites: {total_sites}")
    part1_text.append(f"Total adjusted sites (with pseudocounts +1 per base): {total_adjusted}")
    part1_text.append(f"Background frequency F(b, o): 0.25 for all bases")
    part1_text.append("")

    part1_text.append(format_counts_matrix(counts))
    part1_text.append(format_matrix(freq, "Frequency Matrix F(b, j) = count(b, j) / 27"))
    part1_text.append(format_matrix(freq_prime, "Adjusted Frequency Matrix F'(b, j) = (count(b, j) + 1) / 31"))
    part1_text.append(format_matrix(weight_matrix, "Weight Matrix (PWM) W(b, j) = log2[ F'(b, j) / 0.25 ]"))

    part1_content = "\n".join(part1_text)

    # Print Part 1 to console
    print(part1_content)

    # Write Part 1 to file
    with open(part1_output, 'w') as f:
        f.write(part1_content)
    print(f"\n>>> Part 1 output saved to: {part1_output}")

    # ========================================================================
    #  PART 2: Scan sequences and report top 30 gene IDs
    # ========================================================================

    # Step 5: Read upstream regulatory sequences
    sequences = read_fasta_sequences(sequences_file)
    print(f"\nNumber of sequences loaded: {len(sequences)}")

    # Step 6: Score each sequence (both forward and reverse complement strands)
    results = []
    for gene_id, seq in sequences:
        if len(seq) < motif_length:
            continue
        best_score, best_pos, best_subseq, strand = scan_sequence(seq, weight_matrix, motif_length)
        results.append((gene_id, best_score, best_pos, best_subseq, strand))

    # Sort by score descending
    results.sort(key=lambda x: x[1], reverse=True)

    # --- Build Part 2 output ---
    part2_text = []
    part2_text.append("PART 2 — MOTIF SCANNING")
    part2_text.append("Top 30 Predicted argR Binding Sites from Upstream Regulatory Regions")
    part2_text.append(f"Total sequences scanned: {len(sequences)}")
    part2_text.append(f"Both forward (+) and reverse complement (-) strands scanned")
    part2_text.append(f"Motif length: {motif_length} bases")
    part2_text.append("")
    part2_text.append(f"{'='*90}")
    part2_text.append(f"{'Rank':<6}{'Gene ID':<15}{'Score':>10}{'Position':>10}{'Strand':>8}   {'Binding Site Sequence'}")
    part2_text.append(f"{'-'*90}")

    for rank, (gene_id, score, pos, subseq, strand) in enumerate(results[:30], 1):
        part2_text.append(f"{rank:<6}{gene_id:<15}{score:>10.4f}{pos:>10}{strand:>8}   {subseq.upper()}")

    part2_text.append(f"{'='*90}")

    part2_content = "\n".join(part2_text)

    # Print Part 2 to console
    print(f"\n{part2_content}")

    # Write Part 2 to file
    with open(part2_output, 'w') as f:
        f.write(part2_content)
    print(f"\n>>> Part 2 output saved to: {part2_output}")

    # --- Summary ---
    print(f"\n{'='*70}")
    print("  OUTPUT FILES GENERATED")
    print(f"{'='*70}")
    print(f"  Part 1: {part1_output}")
    print(f"  Part 2: {part2_output}")
    print(f"{'='*70}")
    print("  Done.")


if __name__ == "__main__":
    main()
