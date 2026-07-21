# PWM Motif Scanner — argR Binding Site Prediction

![Python](https://img.shields.io/badge/python-3.x-blue)
![License](https://img.shields.io/badge/license-MIT-green)

Builds a Position Weight Matrix (PWM) from known transcription-factor binding sites
and uses it to scan a genome's upstream regulatory regions for new candidate binding
sites — a classic motif-discovery approach in bioinformatics.

## Overview

Given 27 known binding sites for **argR** (an arginine-responsive transcriptional
regulator in *E. coli* K-12 MG1655), this project:

1. Builds a **PWM** (log-odds weight matrix) describing the base preference at each
   of the motif's 18 positions
2. Scans the upstream regulatory region (400 bp upstream to 50 bp downstream of the
   start codon) of 4,319 *E. coli* genes for the strongest match to that motif, on
   **both** the forward and reverse-complement strand
3. Ranks all genes by their best score and reports the top 30 — the genes most
   likely to be directly regulated by argR

## Approach

**Part 1 — build the PWM**

- Frequency matrix: `F(b,j) = count(b,j) / 27`
- Pseudocount-adjusted frequency matrix (avoids `log(0)`):
  `F'(b,j) = (count(b,j) + 1) / 31`
- Log-odds weight matrix against a uniform 25% background:
  `W(b,j) = log2[F'(b,j) / 0.25]`

**Part 2 — scan for matches**

- Slide the 18 bp PWM window across every position of each upstream sequence, on
  both strands (transcription factors can bind either orientation)
- Score a window as the sum of the matrix values for the bases observed
- Keep each gene's single best-scoring window (position, strand, sequence)
- Rank all genes by that score and report the top 30

## Results

The top predicted argR binding sites, ranked by PWM score:

| Rank | Gene ID | Score | Position | Strand | Sequence |
|---:|---|---:|---:|:---:|---|
| 1 | b3171 | 20.98 | 102 | + | TCACTGAATTTTTATGCA |
| 2 | 16131063 | 20.98 | 101 | − | TCACTGAATTTTTATGCA |
| 3 | 16128258 | 20.76 | 353 | + | AAAGTGAATTTTAATTCA |
| 4 | 16128828 | 20.53 | 119 | − | CAAGTGAATTTATATGCA |
| 5 | 16132076 | 19.92 | 355 | + | AAATTGAATTTTAATTCA |

The top hits share a strong `TGAAT...ATGCA`-like core, consistent with the expected
argR consensus. Several top pairs (e.g. ranks 1–2, ranks 5–6) are adjacent genes in
the genome whose upstream regions overlap — each picks up the same physical binding
site from its own strand, which is exactly what you'd expect from a shared, divergent
promoter sitting between two genes transcribed in opposite directions.

Full matrices and the complete top-30 list are in `results/`.

## Project structure

```
.
├── data/
│   ├── argR-counts-matrix.txt        # Counts matrix from 27 known binding sites
│   └── E_coli_K12_MG1655.400_50      # Upstream regions for 4,319 genes
├── src/
│   └── pwm_scan.py                   # Full pipeline: build PWM -> scan -> rank
├── results/
│   ├── Part1_Matrices_Output.txt     # Counts, frequency, and weight matrices
│   └── Part2_Top30_Binding_Sites.txt # Top 30 predicted binding sites
├── LICENSE
```

## Getting started

No external dependencies — just the Python standard library (`math`, `os`, `argparse`).

```bash
git clone https://github.com/<your-username>/pwm-motif-scanner.git
cd pwm-motif-scanner
python src/pwm_scan.py
```

With no arguments it uses the files in `data/` and writes results to `results/`.
To use different input files or output location:

```bash
python src/pwm_scan.py --counts path/to/counts.txt --sequences path/to/seqs.txt -o path/to/output_dir
```

## Method notes

- Both strands are scanned since a transcription factor can bind in either
  orientation relative to the gene it regulates.
- Unknown bases (e.g. `n`) contribute a score of 0 rather than being excluded.
- The sequence parser supports both the provided custom format
  (`gene_id \ sequence`) and standard FASTA.

## Tech stack

- Python (standard library only)

## License

This project is licensed under the [MIT License](LICENSE).
