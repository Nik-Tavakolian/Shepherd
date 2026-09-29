"""Tests for reading and writing files."""

import csv

import pytest

from shepherd.errors import ShepherdError
from shepherd.io import output_path, read_counts, write_count_table


def test_read_counts_splits_sequences_by_length(tmp_path):
    path = tmp_path / 'reads.txt'
    path.write_text('ACGT\t5\nACG 2\n\nACGTA\t1\t\nAC\t7\nACGA\t3\n')
    reads = read_counts(path, 4)
    assert reads.barcodes == {'ACGT': 5, 'ACGA': 3}
    assert reads.deletions == {'ACG': 2}
    assert reads.insertions == {'ACGTA': 1}


def test_read_counts_reports_malformed_lines(tmp_path):
    path = tmp_path / 'reads.txt'
    path.write_text('ACGT\t5\nACGT five\n')
    with pytest.raises(ShepherdError, match='line 2'):
        read_counts(path, 4)


@pytest.mark.parametrize(
    ('input_path', 'expected'),
    [('t0.txt', 't0_pb_freq.csv'), ('data/reads.t0.txt', 'data/reads.t0_pb_freq.csv')],
)
def test_output_path(input_path, expected):
    assert str(output_path(input_path, '_pb_freq.csv')) == expected


def test_write_count_table_fills_missing_time_points_with_zero(tmp_path):
    path = tmp_path / 'multi_freqs.csv'
    write_count_table(path, [{'AAAA': 5, 'CCCC': 2}, {'CCCC': 3, 'GGGG': 7}])
    with open(path, newline='') as fh:
        assert list(csv.reader(fh)) == [
            ['barcode', 'time_point_1', 'time_point_2'],
            ['AAAA', '5', '0'],
            ['CCCC', '2', '3'],
            ['GGGG', '0', '7'],
        ]
