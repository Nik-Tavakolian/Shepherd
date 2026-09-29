"""Shared helpers for the Shepherd regression tests.

The tests run shepherd_t0.py and shepherd_multi.py end to end on synthetic
data generated as described in Supplementary Section 4.A: random barcodes
with exponentially distributed read counts and a constant per nucleotide
substitution error rate.
"""
import csv
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
BARCODE_LENGTH = 20


def simulate_counts(n_barcodes, error_rate, seed):
    """Return {sequence: read count} for random barcodes with substitution errors."""
    rng = np.random.default_rng(seed)
    barcodes = rng.integers(0, 4, size=(n_barcodes, BARCODE_LENGTH), dtype=np.uint8)
    counts = np.ceil(rng.exponential(100, n_barcodes)).astype(int)
    reads = np.repeat(barcodes, counts, axis=0)
    is_error = rng.random(reads.shape) < error_rate
    shift = rng.integers(1, 4, size=reads.shape, dtype=np.uint8)
    reads = np.where(is_error, (reads + shift) % 4, reads)
    alphabet = np.frombuffer(b'ACGT', dtype=np.uint8)
    seqs = alphabet[reads].view(f'S{BARCODE_LENGTH}').ravel()
    unique_seqs, seq_counts = np.unique(seqs, return_counts=True)
    return {seq.decode(): int(count) for seq, count in zip(unique_seqs, seq_counts)}


def hamming(seq_1, seq_2):
    return sum(n_1 != n_2 for n_1, n_2 in zip(seq_1, seq_2))


@pytest.fixture(scope='session')
def background():
    """About 19 000 unique sequences from 2 000 barcodes with a 0.5% error rate."""
    return simulate_counts(n_barcodes=2000, error_rate=0.005, seed=0)


def write_counts(path, background, extra):
    """Write background plus extra sequences in the Shepherd input format."""
    assert not set(extra) & set(background), 'test sequence collides with background'
    with open(path, 'w') as fh:
        for seq, count in {**background, **extra}.items():
            fh.write(f'{seq}\t{count}\n')


def run_t0(workdir, filename, *args):
    subprocess.run([sys.executable, str(REPO_ROOT / 'shepherd_t0.py'), '-f', filename,
                    '-l', str(BARCODE_LENGTH), *args],
                   cwd=workdir, check=True, capture_output=True)


def run_multi(workdir, f0, later_files):
    subprocess.run([sys.executable, str(REPO_ROOT / 'shepherd_multi.py'), '-f0', f0,
                    '-fn', *later_files],
                   cwd=workdir, check=True, capture_output=True)


def read_multi_freqs(workdir):
    """Return {barcode: [count at t0, count at t1, ...]} from multi_freqs.csv."""
    with open(Path(workdir) / 'multi_freqs.csv') as fh:
        rows = list(csv.reader(fh))[1:]
    return {row[0]: [int(count) for count in row[1:]] for row in rows}


def read_seq_clust(path):
    """Return {sequence: cluster label} from a *_seq_clust.csv file."""
    with open(path) as fh:
        return dict(list(csv.reader(fh))[1:])
