# Shepherd

[![CI](https://github.com/Nik-Tavakolian/Shepherd/actions/workflows/ci.yml/badge.svg)](https://github.com/Nik-Tavakolian/Shepherd/actions/workflows/ci.yml)
[![Paper](https://img.shields.io/badge/Bioinformatics-10.1093%2Fbioinformatics%2Fbtac395-blue)](https://doi.org/10.1093/bioinformatics/btac395)
[![License: GPL v3](https://img.shields.io/badge/license-GPLv3-blue)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.13%20%7C%203.14-blue)](pyproject.toml)

## Getting Started

Shepherd is a Python program for correcting substitution errors and single insertion and deletion errors in DNA barcode reads. These errors occur during PCR amplification and sequencing of the DNA barcodes. Shepherd is cross-platform and runs on any computer with Python 3.10 or later.

### Installation

Clone the repository and install the package, which also installs its only dependency, SciPy:

```bash
git clone https://github.com/Nik-Tavakolian/Shepherd.git
cd Shepherd
pip install .
```

This provides the `shepherd` command with two subcommands, `shepherd cluster` and `shepherd track`, described below. `python -m shepherd` works as well.

**Upgrading from version 1:** `shepherd cluster` replaces `python3 shepherd_t0.py` and `shepherd track` replaces `python3 shepherd_multi.py`. All options are unchanged. The original scripts are available under the `v1.0.0` tag.

### shepherd cluster

This command is designed to cluster the sequencing reads from a single time point to correct substitution errors and single insertion and deletion errors.

**IMPORTANT NOTE:** Shepherd will try to estimate the error rate from the input file. However, since we are estimating a small probability we need a large number of input sequences to do so accurately. If your data has fewer than 10 000 sequences or if you observe unrealistic error rate estimates we suggest providing Shepherd with an error rate estimate using the optional input parameter **-e**.

### Inputs

#### Required Inputs

These inputs must be provided to run the command.

**-l, --length:** (integer) The correct barcode length.

**-f, --input:** (.txt file) The input file with a sequence and a sequence count in each row, separated by whitespace. Currently this is the only input file format supported by Shepherd. Only sequences in the file with lengths l (correct barcode length), l + 1 (single insertion errors) and l - 1 (single deletion errors) will be processed by Shepherd.

    Example file:   testdata_t0.txt

                    TCCCTTACTAATCGAAGAAG	5
                    ATAGTATGGATCTGGACCGC	10
                    ATCCAGTGCTAGTTCAACTC	3
                    AATTTTGGAACAGGCCGTAG	200

#### Optional Inputs

These inputs are optional and we recommend using the default values determined by Shepherd.

**-e, --error-rate:** (float) An estimate of the substitution error rate of the sequencing protocol used to generate the input data. This is a floating point number, e.g. 0.01 if the estimated error rate is 1%. If not provided this parameter is automatically determined based on the input data.

**-eps, --epsilon:** (integer) The maximum Hamming distance considered for merging two sequences into the same cluster. If not provided this parameter is automatically determined based on the input data.

**-k, --kmer-length:** (integer) The substring length used to divide the sequences into partitions. If not provided this parameter is automatically determined based on the input data.

**-tau, --tau:** (integer) A distance threshold for frequency 1 sequences that determines if they should be merged with another sequence. If a frequency 1 sequence has Hamming distance less than or equal to this threshold to a candidate sequence it will be merged. If not provided this parameter is automatically determined based on the input data.

**-ft, --count-threshold:** (integer) A count threshold for defining true barcodes. Any sequence with at least this many reads is a putative barcode. If not provided this parameter is automatically determined based on the input data.

**-bft, --log-bf-threshold:** (float) The threshold for log Bayes factor. The default value is -4.

**-Nh, --n-top:** (integer) Number of high-count sequences used for estimation of the substitution error rate. The default value is 500 (see Implementation notes).

### Outputs

**_seq_clust.csv:** A .csv file where the unique sequences are in the first column and the cluster labels are in the second column.

**_pb_freq.csv:** A .csv file where the putative barcodes are in the first column and the estimated counts are in the second column.

**_params.json:** The parameters used for the run (error rate, epsilon, tau, count threshold, k-mer length, ...) in JSON format. `shepherd track` reads them from this file.

### Usage

**Command line usage example:** <code>shepherd cluster -f testdata_t0.txt -l 20 -e 0.01</code>

### shepherd track

This command is designed to use the the clustering from the first time point, i.e., the outputs of shepherd cluster, to estimate the counts of the putative barcodes at later time points, given the sequencing reads from each time point. If new barcodes that did not appear in the first time point emerge in later time points, the program is capable of identifying and tracking them. Note that shepherd cluster must be run on the first time point in the same folder before running shepherd track.

### Inputs

**-f0, --first:** (.txt file) The same input file used to run shepherd cluster containing the sequences and the sequence counts.

**-fn, --later:** (.txt files) Space separated list of .txt files containing the sequences and sequence counts for each time point. These files should have the same format as the input file to shepherd cluster (see testdata_t0.txt) and should be ordered by time point (see usage example below).\

**-o, --output:** (string) The prefix of the final output file. By default set to 'multi_freqs' which produces an output file called 'multi_freqs.csv'.

### Outputs

**multi_freqs.csv:** A .csv file where each row is a putative barcode and the columns give the estimated counts for each time point.

**seq_clust.csv:** A .csv file for each time point where each row contains a sequence and the cluster ID it was assigned.

### Usage

**Command line usage example:**\
<code>shepherd track -f0 testdata_t0.txt -fn testdata_t1.txt testdata_t2.txt</code>

## Implementation notes

The code follows the method described in the paper and its Supplementary Material, with these differences. The first four do not change the results.

- **The k-mer Index only contains putative barcodes.** The paper builds the index from all sequences and then keeps the putative barcodes among the neighbours of a sequence (Algorithm S1). Since sequences are processed in descending order of read count, `cluster_reads` instead adds each sequence to the index when it becomes a putative barcode, which gives the same neighbours directly with a much smaller index.
- **The log Bayes factor** is computed with `math.lgamma` instead of `scipy.stats.binom.logpmf`. The formula is the same; it avoids SciPy's overhead per call.
- **Ties.** A sequence at the same distance from two putative barcodes with the same read count joins the alphabetically first one. Shepherd 1.x picked one arbitrarily, so the result could change between runs.
- **Parameter selection** (Supplementary Section 3) classifies a sequence as a true barcode when ln K < 0, regardless of the `-bft` threshold used during clustering. The error rate is estimated from the 501 (`-Nh` + 1) highest-count sequences, as in Shepherd 1.x.
- **Sequences at distance 1** from their closest putative barcode are merged without the Bayesian test if they have fewer reads than the count threshold. This shortcut is not described in the paper and saves most of the tests. It can only differ from the full test when two true barcodes one substitution apart both have fewer reads than the count threshold, which is rare for random barcodes.
- **Emerging barcodes in `shepherd track`.** A cluster is not tested for emerging barcodes if one of its sequences has more reads than the putative barcode itself at that time point.

## Development

Install the package in editable mode with the development tools, and set up the pre-commit hooks:

```bash
pip install -e ".[dev]"
pre-commit install
```

Run the tests, the linter and the type checker:

```bash
pytest
ruff check .
mypy
```

The tests include golden-output tests (`tests/test_golden.py`) that check that the results on a fixed synthetic data set are unchanged. If a change of results is intended, re-record them with `python tests/golden/make_golden.py` and explain the change in the commit message.

## Citation

If you use Shepherd, please cite:

Tavakolian N, Frazão JG, Bendixsen D, Stelkens R, Li C-B. Shepherd: accurate clustering for correcting DNA barcode errors. *Bioinformatics* 38(15):3710–3716 (2022). https://doi.org/10.1093/bioinformatics/btac395

```bibtex
@article{tavakolian2022shepherd,
  title   = {Shepherd: accurate clustering for correcting {DNA} barcode errors},
  author  = {Tavakolian, Nik and Fraz{\~a}o, Jo{\~a}o Guilherme and Bendixsen, Devin and Stelkens, Rike and Li, Chun-Biu},
  journal = {Bioinformatics},
  volume  = {38},
  number  = {15},
  pages   = {3710--3716},
  year    = {2022},
  doi     = {10.1093/bioinformatics/btac395}
}
```

## License

Copyright (C) 2021-2026 Nik Tavakolian

Shepherd is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed in the hope that it will be useful, but without any warranty. See [LICENSE](LICENSE) for the full text.
