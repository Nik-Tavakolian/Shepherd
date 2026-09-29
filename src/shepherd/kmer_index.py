"""The k-mer Index (Section 2.1 of the paper)."""

from collections.abc import Iterator
from itertools import combinations

KmerCombination = tuple[tuple[int, str], ...]


class KmerIndex:
    """Find the indexed sequences that may lie within Hamming distance epsilon of a query.

    Every sequence of length ``barcode_length`` is split into non-overlapping
    k-mers of length ``kmer_length`` (the last one may be shorter), each tagged
    with its position. By the pigeonhole principle, two sequences within Hamming
    distance epsilon share at least ``n_partitions - epsilon`` of these k-mers.
    The index maps every combination of that many position-tagged k-mers to the
    indexed sequences containing it, so any indexed sequence within epsilon of a
    query shares at least one combination with the query.

    :meth:`neighbours` can also return sequences further away than epsilon (up
    to about ``kmer_length * epsilon``); callers filter them by distance.
    """

    def __init__(
        self, barcode_length: int, kmer_length: int, n_partitions: int, epsilon: int
    ) -> None:
        if n_partitions <= epsilon:
            raise ValueError(
                f'the number of partitions ({n_partitions}) must be larger than epsilon '
                f'({epsilon}); choose a smaller k-mer length'
            )
        self._kmer_length = kmer_length
        self._positions = [
            (start // kmer_length + 1, start) for start in range(0, barcode_length, kmer_length)
        ]
        self._combination_size = n_partitions - epsilon
        self._sequences_by_combination: dict[KmerCombination, set[str]] = {}

    def _combinations(self, seq: str) -> Iterator[KmerCombination]:
        k = self._kmer_length
        kmers = [(position, seq[start : start + k]) for position, start in self._positions]
        return combinations(kmers, self._combination_size)

    def add(self, seq: str) -> None:
        """Add a sequence to the index."""
        for combination in self._combinations(seq):
            sequences = self._sequences_by_combination.get(combination)
            if sequences is None:
                self._sequences_by_combination[combination] = {seq}
            else:
                sequences.add(seq)

    def neighbours(self, seq: str) -> set[str]:
        """Return the indexed sequences that share at least one k-mer combination with seq.

        This includes every indexed sequence within Hamming distance epsilon of seq.
        """
        neighbours: set[str] = set()
        for combination in self._combinations(seq):
            sequences = self._sequences_by_combination.get(combination)
            if sequences is not None:
                neighbours.update(sequences)
        return neighbours
