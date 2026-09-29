"""Create the golden input files and record the expected outputs.

The golden tests check that refactoring does not change Shepherd's output.
Only rerun this script when a change of output is intended, and explain the
change in the commit message.

    python tests/golden/make_golden.py            # re-record outputs
    python tests/golden/make_golden.py --inputs   # also re-simulate inputs
"""

import argparse
import gzip
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

GOLDEN_DIR = Path(__file__).resolve().parent
DATA_DIR = GOLDEN_DIR / 'data'
sys.path.insert(0, str(GOLDEN_DIR.parent))

from conftest import run_cluster, run_track  # noqa: E402

TIME_POINTS = ['t0.txt', 't1.txt', 't2.txt']

# (name, input file, extra command line arguments)
SINGLE_SCENARIOS = [
    ('default', 't0.txt', []),
    ('k3', 't0_k3.txt', ['-e', '0.01', '-k', '3']),
]
OUTPUT_FILES = [
    't0_seq_clust.csv',
    't0_pb_freq.csv',
    't0_k3_seq_clust.csv',
    't0_k3_pb_freq.csv',
    't1_seq_clust.csv',
    't2_seq_clust.csv',
    'multi_freqs.csv',
]

PARAMETER_FILES = ['t0_params.json', 't0_k3_params.json']


def simulate_time_series(seed=2022, n_barcodes=500, length=20, error_rate=0.005, indel_rate=0.001):
    """Simulate three time points of barcode read counts.

    Counts drift between time points and three barcodes emerge after t0: one
    far from all others and two at distance 2 from an existing barcode (one
    next to a small barcode, one next to the largest barcode). A small
    fraction of reads carry a single insertion or deletion.
    """
    rng = np.random.default_rng(seed)
    barcodes = rng.integers(0, 4, size=(n_barcodes, length), dtype=np.uint8)
    counts = np.ceil(rng.exponential(100, n_barcodes))

    near_small = barcodes[0].copy()
    near_small[:2] = (near_small[:2] + 1) % 4
    near_largest = barcodes[np.argmax(counts)].copy()
    near_largest[:2] = (near_largest[:2] + 1) % 4
    far_from_all = rng.integers(0, 4, size=length, dtype=np.uint8)
    barcodes = np.vstack([barcodes, near_small, near_largest, far_from_all])
    emerging_counts = [np.zeros(3), np.array([200, 200, 300]), np.array([400, 400, 600])]

    alphabet = np.frombuffer(b'ACGT', dtype=np.uint8)
    series = []
    for t in range(3):
        drift = counts if t == 0 else rng.poisson(counts * rng.lognormal(0, 0.3 * t, n_barcodes))
        reads = np.repeat(barcodes, np.concatenate([drift, emerging_counts[t]]).astype(int), axis=0)
        is_error = rng.random(reads.shape) < error_rate
        shift = rng.integers(1, 4, size=reads.shape, dtype=np.uint8)
        reads = np.where(is_error, (reads + shift) % 4, reads)
        seqs = [s.decode() for s in alphabet[reads].view(f'S{length}').ravel()]
        for i in np.nonzero(rng.random(len(seqs)) < indel_rate)[0]:
            pos = rng.integers(length)
            if rng.random() < 0.5:
                seqs[i] = seqs[i][:pos] + seqs[i][pos + 1 :]
            else:
                seqs[i] = seqs[i][:pos] + 'ACGT'[rng.integers(4)] + seqs[i][pos:]
        unique_seqs, seq_counts = np.unique(seqs, return_counts=True)
        series.append(dict(zip(unique_seqs.tolist(), seq_counts.tolist(), strict=True)))
    return series


def write_gzip(path, data):
    # A fixed timestamp keeps the file identical when its content does not change.
    with open(path, 'wb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', mtime=0) as fh:
        fh.write(data)


def write_inputs():
    DATA_DIR.mkdir(exist_ok=True)
    for filename, counts in zip(TIME_POINTS, simulate_time_series(), strict=True):
        text = ''.join(f'{seq}\t{count}\n' for seq, count in counts.items())
        write_gzip(DATA_DIR / (filename + '.gz'), text.encode())


def copy_inputs(workdir):
    """Decompress the golden inputs into workdir."""
    for filename in TIME_POINTS:
        with (
            gzip.open(DATA_DIR / (filename + '.gz'), 'rt') as src,
            open(Path(workdir) / filename, 'w') as dst,
        ):
            shutil.copyfileobj(src, dst)
    shutil.copy(Path(workdir) / 't0.txt', Path(workdir) / 't0_k3.txt')


def run_all_scenarios(workdir):
    copy_inputs(workdir)
    for _, filename, args in SINGLE_SCENARIOS:
        run_cluster(workdir, filename, *args)
    run_track(workdir, 't0.txt', TIME_POINTS[1:])


def write_outputs():
    with tempfile.TemporaryDirectory() as workdir:
        run_all_scenarios(workdir)
        for filename in OUTPUT_FILES:
            write_gzip(DATA_DIR / (filename + '.gz'), (Path(workdir) / filename).read_bytes())
        for filename in PARAMETER_FILES:
            shutil.copy(Path(workdir) / filename, DATA_DIR / filename)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--inputs', action='store_true', help='also re-simulate the input files')
    if parser.parse_args().inputs:
        write_inputs()
    write_outputs()
