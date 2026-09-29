"""Helpers for working with barcode sequences and their read counts."""

from collections.abc import Iterator, Mapping

NUCLEOTIDES = 'ACGT'


def sort_by_count(counts: Mapping[str, int]) -> list[str]:
    """Sequences in descending order of read count; equal counts keep their input order."""
    return sorted(counts, key=counts.__getitem__, reverse=True)


def single_substitutions(seq: str) -> Iterator[str]:
    """All sequences that differ from seq at exactly one position."""
    for position, current in enumerate(seq):
        for nucleotide in NUCLEOTIDES:
            if nucleotide != current:
                yield seq[:position] + nucleotide + seq[position + 1 :]


def single_deletions(seq: str) -> Iterator[str]:
    """The sequences obtained by deleting one nucleotide, from the first position to the last."""
    for position in range(len(seq)):
        yield seq[:position] + seq[position + 1 :]


def single_insertions(seq: str) -> Iterator[str]:
    """The sequences obtained by inserting one nucleotide, position by position, A, C, G, T."""
    for position in range(len(seq) + 1):
        for nucleotide in NUCLEOTIDES:
            yield seq[:position] + nucleotide + seq[position:]


def truncated_hamming_distance(seq_1: str, seq_2: str, max_distance: int) -> int:
    """The Hamming distance (Eq. 1) if it is at most max_distance, otherwise len(seq_1) (Eq. S20).

    Comparison stops as soon as max_distance is exceeded, which for unrelated
    sequences happens after a few positions.
    """
    distance = 0
    for nucleotide_1, nucleotide_2 in zip(seq_1, seq_2, strict=True):
        if nucleotide_1 != nucleotide_2:
            distance += 1
            if distance > max_distance:
                return len(seq_1)
    return distance
