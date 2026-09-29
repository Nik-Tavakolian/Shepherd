"""Tests for the k-mer Index."""

import random

import pytest

from shepherd.kmer_index import KmerIndex


def mutate(seq, n_substitutions, rng):
    seq = list(seq)
    for position in rng.sample(range(len(seq)), n_substitutions):
        seq[position] = rng.choice([n for n in 'ACGT' if n != seq[position]])
    return ''.join(seq)


def hamming(seq_1, seq_2):
    return sum(n_1 != n_2 for n_1, n_2 in zip(seq_1, seq_2, strict=True))


@pytest.mark.parametrize(
    ('barcode_length', 'kmer_length', 'n_partitions', 'epsilon'),
    [(20, 4, 5, 3), (20, 3, 7, 3), (26, 4, 7, 4), (20, 5, 4, 2), (12, 2, 6, 5)],
)
def test_neighbours_contain_every_sequence_within_epsilon(
    barcode_length, kmer_length, n_partitions, epsilon
):
    # The pigeonhole guarantee of Section 2.1, checked against brute force.
    rng = random.Random(0)
    barcodes = [''.join(rng.choices('ACGT', k=barcode_length)) for _ in range(50)]
    indexed = barcodes + [mutate(b, rng.randint(1, epsilon), rng) for b in barcodes]
    index = KmerIndex(barcode_length, kmer_length, n_partitions, epsilon)
    for seq in indexed:
        index.add(seq)

    for barcode in barcodes:
        query = mutate(barcode, rng.randint(0, epsilon), rng)
        within_epsilon = {seq for seq in indexed if hamming(seq, query) <= epsilon}
        assert within_epsilon <= index.neighbours(query)


def test_neighbours_are_indexed_sequences_only():
    index = KmerIndex(20, 4, 5, 3)
    index.add('ACGTACGTACGTACGTACGT')
    assert index.neighbours('ACGTACGTACGTACGTTTTT') == {'ACGTACGTACGTACGTACGT'}
    assert index.neighbours('TTTTTTTTTTTTTTTTTTTT') == set()


def test_an_added_sequence_is_its_own_neighbour():
    index = KmerIndex(20, 4, 5, 3)
    index.add('ACGTACGTACGTACGTACGT')
    assert index.neighbours('ACGTACGTACGTACGTACGT') == {'ACGTACGTACGTACGTACGT'}


def test_epsilon_must_be_smaller_than_the_number_of_partitions():
    with pytest.raises(ValueError, match='must be larger than epsilon'):
        KmerIndex(20, 10, 2, 3)
