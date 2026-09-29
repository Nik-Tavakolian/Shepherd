# Changelog

## 2.0.0 (unreleased)

Shepherd is now an installable package with a single `shepherd` command. The
clustering method is unchanged; see "Implementation notes" in the README for
where the code differs from the paper's description.

### Fixed

- `shepherd track` (multiple time points):
  - Barcodes found after the first time point were stored in the k-mer Index
    under their nucleotides (`set(seq)` instead of `{seq}`), so later reads
    could often not find them.
  - Later time points used the index of the first time point, which only kept
    k-mer combinations shared by two or more sequences. With a sparse index
    (for example `-k 3`) reads close to an existing barcode were left
    unassigned.
  - Unassigned reads were clustered without an index of their own, so every
    unassigned read became a candidate barcode. Error sequences of a barcode
    that emerged after the first time point were reported as separate
    lineages, and the barcode was undercounted.
  - When an emerging barcode was separated, reads moved to it could be
    separated again, which counted them twice.
- Ties between equally close putative barcodes with the same read count are
  broken alphabetically. Before, the choice could change from run to run.

### Changed

- `shepherd cluster` and `shepherd track` replace `python3 shepherd_t0.py` and
  `python3 shepherd_multi.py`. All options are unchanged and each also has a
  long name, for example `-l/--length`.
- The parameters are written to `<input>_params.json` instead of the pickle
  `<input>_params`. The `<input>_index` pickle is no longer written.
- Rows of `multi_freqs.csv` are in the order in which the barcodes are first
  seen, instead of an order that varied from run to run.
- Errors in the input are reported with a message instead of a traceback.
- Faster and leaner: on 937 000 unique sequences clustering takes 11 s instead
  of 30 s, and with `-k 3` 1.6 GB of memory instead of 8.6 GB. `shepherd track`
  is about twice as fast.
- Requires Python 3.10 or later.

### Added

- A Python API: `read_counts`, `estimate_parameters`, `cluster` and `Tracker`.
- Tests, including golden-output tests on simulated data, and continuous
  integration on Linux, macOS and Windows.

## 1.0.0

The two scripts `shepherd_t0.py` and `shepherd_multi.py`. The code used for the
paper is tagged `paper-2022`.
