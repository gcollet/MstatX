# MstatX

A small, easy-to-extend command-line tool that computes per-column (and
optionally global) conservation/variability statistics from a multiple
sequence alignment.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Contents

- [About](#about)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Available statistics](#available-statistics)
- [Command-line options](#command-line-options)
- [Scoring matrices](#scoring-matrices)
- [Background distributions (jensen)](#background-distributions-jensen)
- [Running the tests](#running-the-tests)
- [Roadmap](#roadmap)
- [Citing MstatX](#citing-mstatx)
- [References](#references)
- [License](#license)

## About

A multiple alignment can be built with tools like Clustal W [[Larkin,
2007]](#references) or MUSCLE [[Edgar, 2004]](#references), or by hand.
Alignment tools usually give a single, global alignment score, but you
often need something more specific: the conservation of one or a few
columns, or a score computed a different way. MstatX fills that gap: a
simple command-line tool dedicated to computing per-column statistics
from an already-built alignment.

MstatX is not an alignment tool and does not validate its input - it
assumes the multiple alignment it's given is correct.

## Installation

Requirements: a C++17 compiler (`g++`).

```sh
git clone https://github.com/gcollet/MstatX.git
cd MstatX
make
```

This builds the `mstatx` binary at the repository root.

## Quick start

```sh
./mstatx -i example/valdar.mali -s trident -o result.txt
```

Computes the `trident` statistic for each column of the alignment in
`example/valdar.mali`, one score per line in `result.txt` (line 1 = score
of column 1, etc). Add `-g` to get a single global score (the mean of
the per-column scores) instead.

The default statistic (if `-s` is omitted) is `wentropy`.

## Available statistics

| Name | What it measures | Reference |
|---|---|---|
| `gap` | Fraction of gaps in the column: `nb_gaps(x) / N` | - |
| `kabat` | Variability: `k / n1`, with `k` the number of distinct symbol types in the column and `n1` the count of the most frequent one | [Wu & Kabat, 1970](#references) |
| `wentropy` | Shannon entropy of the column, weighted by Henikoff & Henikoff sequence weights | [Valdar, 2002](#references), [Henikoff & Henikoff, 1994](#references) |
| `trident` | Combines three factors: weighted entropy, a stereochemical divergence score built from a substitution matrix, and the gap fraction | [Valdar, 2002](#references) |
| `jensen` | Jensen-Shannon divergence between the column's (weighted) amino acid composition and a background distribution | [Capra & Singh, 2007](#references) |
| `mvector` | Mean normalized substitution score of the column against every amino acid of the scoring matrix's alphabet (one value per column *and* per amino acid, not a single score) | MstatX-specific |

For a broader comparison of conservation/variability scores in general,
see [Johansson & Toh, 2010](#references).

## Command-line options

```
./mstatx -h
```

prints the full, current list. Summary:

| Flag | Description | Default |
|---|---|---|
| `-i`, `--input` | MSA input file name (required) | - |
| `-s`, `--statistic` | Statistic to compute (see table above) | `wentropy` |
| `-o`, `--output` | Output file name | `output.txt` |
| `-g`, `--global` | Output a single global score (mean of column scores) instead of one per column | off |
| `-m`, `--matrix` | Substitution matrix file (AAindex format), used by `trident` and `mvector` | `HENS920102.mat` (BLOSUM62-derived) |
| `-k`, `--background` | Background distribution for `jensen`: `uniform`, `legacy`, or a file path | `legacy` |
| `-n`, `--nb_seq` | Maximum number of sequences read from the input | 500 |
| `-t`, `--threshold` | Threshold used when printing correlations | 0.8 |
| `-a`, `--trident_a` | Factor applied to `t(x)` in `trident` | 1.0 |
| `-b`, `--trident_b` | Factor applied to `r(x)` in `trident` | 0.5 |
| `-c`, `--trident_c` | Factor applied to `g(x)` in `trident` | 3.0 |
| `-v`, `--verbose` | Verbose mode | off |
| `-h`, `--help` | Print usage and exit | - |

`-w`/`--window` also exists but currently has no effect on any
statistic - see [TODO.md](TODO.md).

## Scoring matrices

`trident` and `mvector` compare residues using a substitution matrix,
given in the [AAindex](https://www.genome.jp/aaindex/) format via `-m`.
A BLOSUM62-derived matrix (`data/aaindex/HENS920102.mat`) is bundled and
used by default. AAindex itself hasn't seen a substantive update to its
substitution-matrix section since around 2008; any file in the same
format works with `-m`, so you're not limited to what's bundled.

## Background distributions (jensen)

`jensen` compares each column's composition to a background
distribution via `-k`/`--background`:

- `uniform` - 1/20 for each of the 20 standard amino acids.
- `legacy` (default) - the background table from [Capra & Singh,
  2007](#references), used by `jensen` before `-k` existed. Kept as the
  default so existing results aren't silently changed.
- a file path - one `<letter> <frequency>` pair per line, e.g.:
  ```
  A 0.073
  C 0.025
  ...
  ```

Choice of background is not a neutral detail - see [Johansson & Toh,
2010](#references) for a discussion of how it affects conservation
scores. Named, literature-sourced presets (e.g. a BLOSUM62-derived
background) are on the [roadmap](#roadmap).

## Running the tests

```sh
make test
```

Builds and runs the unit test suite (one binary per module, under
`tests/`), covering the alignment reader, the scoring matrix and
background distribution readers, the argument parser, and all six
statistics against synthetic, hand-verifiable alignments.

## Roadmap

See [TODO.md](TODO.md) for planned additions (currently: a `sumofpairs`
statistic, and named background/matrix presets).

## Citing MstatX

If you use MstatX itself in your research, see
[CITATION.cff](CITATION.cff). If you use one of its statistics, please
cite the corresponding original method below rather than (or in
addition to) MstatX.

## References

- Larkin MA, Blackshields G, Brown NP, et al. **Clustal W and Clustal X
  version 2.0.** *Bioinformatics*. 2007;23(21):2947-2948.
- Edgar RC. **MUSCLE: multiple sequence alignment with high accuracy and
  high throughput.** *Nucleic Acids Research*. 2004;32(5):1792-1797.
- Valdar WSJ. **Scoring residue conservation.** *Proteins: Structure,
  Function, and Bioinformatics*. 2002;48(2):227-241.
- Henikoff S, Henikoff JG. **Position-based sequence weights.** *Journal
  of Molecular Biology*. 1994;243(4):574-578.
- Wu TT, Kabat EA. **An analysis of the sequences of the variable
  regions of Bence Jones proteins and myeloma light chains and their
  implications for antibody complementarity.** *The Journal of
  Experimental Medicine*. 1970;132(2):211-250.
- Capra JA, Singh M. **Predicting functionally important residues from
  sequence conservation.** *Bioinformatics*. 2007;23(15):1875-1882.
- Johansson F, Toh H. **A comparative study of conservation and
  variation scores.** *BMC Bioinformatics*. 2010;11:388.
- Kawashima S, Pokarowski P, Pokarowska M, Kolinski A, Katayama T,
  Kanehisa M. **AAindex: amino acid index database, progress report
  2008.** *Nucleic Acids Research*. 2008;36(suppl 1):D202-D205.

(Full BibTeX entries are in [`doc/biblio.bib`](doc/biblio.bib).)

## License

MIT - see [LICENSE](LICENSE). Questions or bug reports: open a GitHub
issue, or contact guillaume[AT]gcollet.fr.
