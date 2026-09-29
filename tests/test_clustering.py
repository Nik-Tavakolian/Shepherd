"""Tests for clustering a single time point."""

import pytest

from conftest import read_seq_clust, run_cluster, write_counts
from shepherd.clustering import (
    Match,
    cluster,
    cluster_reads,
    find_closest_barcode,
    is_error_sequence,
)
from shepherd.io import ReadCounts
from shepherd.parameters import Parameters

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


# Unit tests on hand-made data. BARCODE_1 has 100 reads; the other sequences
# are single-substitution variants (or unrelated barcodes) with few reads.
PARAMS = Parameters(
    barcode_length=20,
    error_rate=0.005,
    max_count=100,
    epsilon=3,
    tau=2,
    count_threshold=15,
    kmer_length=4,
    n_partitions=5,
)


def test_find_closest_barcode_prefers_distance_then_count():
    counts = {'AAAA': 5, 'AATT': 9, 'CCAA': 7}
    # AAAA (distance 1) beats CCAA (distance 3) despite fewer reads.
    assert find_closest_barcode('AAAT', ['AAAA', 'CCAA'], counts, 3) == Match('AAAA', 1)
    # AAAA and AATT are both at distance 1; AATT has more reads.
    assert find_closest_barcode('AAAT', ['AAAA', 'AATT'], counts, 3) == Match('AATT', 1)
    assert find_closest_barcode('TTTT', counts, counts, 1) is None


def test_find_closest_barcode_breaks_ties_independently_of_order():
    counts = {'AACC': 9, 'CCAA': 9, 'GGGG': 20}
    for candidates in (['AACC', 'CCAA'], ['CCAA', 'AACC']):
        assert find_closest_barcode('ACAC', candidates, counts, 2) == Match('AACC', 2)


def test_is_error_sequence():
    # Arguments: read count, count of the putative barcode, distance.
    assert is_error_sequence(1, 1, 2, PARAMS)  # one read within tau
    assert is_error_sequence(10, 11, 1, PARAMS)  # distance 1, below the count threshold
    # The Bayes test decides the rest (threshold -4):
    assert is_error_sequence(2, 100, 2, PARAMS)  # ln K = 2.4
    assert not is_error_sequence(2, 100, 3, PARAMS)  # ln K = -16.8
    assert not is_error_sequence(10, 11, 2, PARAMS)  # ln K = -96.6


def test_cluster_reads_merges_error_sequences_into_their_barcode():
    error_1 = 'AAAAAAAAAACCCCCCCCCA'  # distance 1 from BARCODE_1
    error_2 = 'AAAAAAAAAACCCCCCCCAA'  # distance 2, one read (within tau)
    other = 'GGGGGGGGGGTTTTTTTTTT'
    counts = {error_1: 3, BARCODE_1: 100, other: 20, error_2: 1}

    clustering = cluster_reads(counts, PARAMS)

    assert clustering.barcode_counts == {BARCODE_1: 104, other: 20}
    assert clustering.labels == {BARCODE_1: 0, other: 1, error_1: 0, error_2: 0}


def test_correct_indels_merges_single_insertions_and_deletions():
    reads = ReadCounts(
        barcodes={BARCODE_1: 100},
        insertions={'AAAAAAAAAACCCCCCCCCCG': 2, 'TTTTTTTTTTTTTTTTTTTTT': 5},
        deletions={'AAAAAAAAACCCCCCCCCC': 3},
    )
    clustering = cluster(reads, PARAMS)

    assert clustering.barcode_counts == {BARCODE_1: 105}
    assert set(clustering.labels) == {
        BARCODE_1,
        'AAAAAAAAAACCCCCCCCCCG',
        'AAAAAAAAACCCCCCCCCC',
    }
