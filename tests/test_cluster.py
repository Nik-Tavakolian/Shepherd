"""Tests for single time point clustering (shepherd cluster)."""

import pytest

from conftest import read_seq_clust, run_cluster, write_counts

BARCODE_1 = 'AAAAAAAAAACCCCCCCCCC'
BARCODE_2 = 'AAAAAAAAAACCCCCCGGGG'  # distance 4 from BARCODE_1
BETWEEN = 'AAAAAAAAAACCCCCCGGCC'  # distance 2 from both barcodes


@pytest.mark.parametrize('first', [BARCODE_1, BARCODE_2])
def test_equally_close_barcodes_with_equal_counts_break_ties_alphabetically(
    tmp_path, background, first
):
    second = BARCODE_2 if first == BARCODE_1 else BARCODE_1
    write_counts(tmp_path / 't0.txt', background, {first: 50, second: 50, BETWEEN: 1})
    run_cluster(tmp_path, 't0.txt')

    labels = read_seq_clust(tmp_path / 't0_seq_clust.csv')
    assert labels[BETWEEN] == labels[BARCODE_1]
