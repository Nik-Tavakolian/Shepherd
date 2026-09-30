"""Check the expected results documented in examples/README.md."""

import shutil
from pathlib import Path

import pytest

from conftest import read_multi_freqs, run_cluster, run_track
from shepherd.io import read_barcode_counts

EXAMPLES = Path(__file__).resolve().parents[1] / 'examples'


@pytest.fixture
def workdir(tmp_path):
    for path in EXAMPLES.glob('*.txt'):
        shutil.copy(path, tmp_path)
    return tmp_path


@pytest.mark.parametrize(
    ('name', 'expected'),
    [
        (
            'single_basic',
            {
                'AAAAACCCCCGGGGGTTTTT': 1010,
                'CCCCCAAAAATTTTTGGGGG': 504,
                'GGGGGTTTTTAAAAACCCCC': 50,
            },
        ),
        (
            'single_close_barcodes',
            {'ACGTACGTACGTACGTACGT': 805, 'ACGTACGTACGTACGTAGGA': 15},
        ),
        (
            'single_indels',
            {'AAAAACCCCCGGGGGTTTTT': 1011, 'CCCCCAAAAATTTTTGGGGG': 200},
        ),
    ],
)
def test_single_time_point_examples(workdir, name, expected):
    run_cluster(workdir, f'{name}.txt', '-e', '0.01')
    assert read_barcode_counts(workdir / f'{name}_pb_freq.csv') == expected


def test_multiple_time_point_example(workdir):
    run_cluster(workdir, 'multi_t0.txt', '-e', '0.01')
    run_track(workdir, 'multi_t0.txt', ['multi_t1.txt', 'multi_t2.txt'])
    assert read_multi_freqs(workdir) == {
        'AAAAACCCCCGGGGGTTTTT': [1006, 905, 802],
        'CCCCCAAAAATTTTTGGGGG': [503, 603, 704],
        'GGGGGTTTTTAAAAACCCCC': [101, 0, 0],
        'AAAAACCCCCGGGGGTTTAA': [0, 300, 350],
        'TTTTTGGGGGCCCCCAAAAA': [0, 204, 405],
    }
