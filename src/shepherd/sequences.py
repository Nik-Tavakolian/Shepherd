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
