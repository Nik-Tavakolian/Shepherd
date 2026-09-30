"""Clustering the reads of one time point (Section 2.2)."""

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import NamedTuple

from shepherd.io import ReadCounts
from shepherd.kmer_index import KmerIndex
from shepherd.model import Parameters


def sort_by_count(counts: Mapping[str, int]) -> list[str]:
    """Sequences in descending order of read count; equal counts keep their input order."""
    return sorted(counts, key=counts.__getitem__, reverse=True)


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


@dataclass
class Clustering:
    """The putative barcodes of one time point and the cluster of every sequence.

    Every putative barcode has its own integer cluster label, and every other
    clustered sequence has the label of the barcode it was merged with.
    """

    barcode_counts: dict[str, int] = field(default_factory=dict)
    """Total read count of the cluster of each putative barcode."""
    labels: dict[str, int] = field(default_factory=dict)
    """Cluster label of each clustered sequence."""

    def add_barcode(self, barcode: str, count: int, label: int) -> None:
        self.labels[barcode] = label
        self.barcode_counts[barcode] = count

    def merge(self, seq: str, count: int, barcode: str) -> None:
        """Add a sequence and its reads to the cluster of a putative barcode."""
        self.labels[seq] = self.labels[barcode]
        self.barcode_counts[barcode] += count

    def split_off(self, seq: str, count: int, barcode: str, label: int) -> None:
        """Turn a sequence in the cluster of barcode into a putative barcode of its own."""
        self.barcode_counts[barcode] -= count
        self.add_barcode(seq, count, label)

    def move(self, seq: str, count: int, from_barcode: str, to_barcode: str) -> None:
        """Move a sequence and its reads from one cluster to another."""
        self.barcode_counts[from_barcode] -= count
        self.barcode_counts[to_barcode] += count
        self.labels[seq] = self.labels[to_barcode]


class Match(NamedTuple):
    barcode: str
    distance: int


def find_closest_barcode(
    seq: str, candidates: Iterable[str], counts: Mapping[str, int], max_distance: int
) -> Match | None:
    """The candidate closest to seq in Hamming distance, if it is within max_distance.

    Equally close candidates are ranked by their count in ``counts`` (higher
    first) and then alphabetically, so the result does not depend on the order
    of the candidates.
    """
    best: tuple[int, int, str] | None = None
    for candidate in candidates:
        distance = truncated_hamming_distance(seq, candidate, max_distance)
        if distance <= max_distance:
            rank = (distance, -counts[candidate], candidate)
            if best is None or rank < best:
                best = rank
    if best is None:
        return None
    distance, _, barcode = best
    return Match(barcode, distance)


def is_error_sequence(count: int, barcode_count: int, distance: int, params: Parameters) -> bool:
    """Whether a sequence with ``count`` reads is an error sequence of a putative barcode.

    The putative barcode has ``barcode_count`` reads and lies at ``distance``
    (at most epsilon) from the sequence. The Bayesian test of Section 2.2.2 is
    skipped where its outcome is clear (Supplementary Section 3.D).
    """
    if count == 1 and distance <= params.tau:
        return True
    if distance == 1:
        # Not described in the paper: below the count threshold, sequences at
        # distance 1 are merged without the test. See "Implementation notes" in
        # the README.
        return True
    log_k = params.error_model.log_bayes_factor(count, barcode_count, distance)
    return log_k > params.log_bf_threshold


def cluster_reads(counts: Mapping[str, int], params: Parameters) -> Clustering:
    """Cluster sequences of the barcode length around putative barcodes (Algorithm S1).

    Sequences are processed in descending order of read count. A sequence with
    fewer than ``params.count_threshold`` reads is merged with the closest
    putative barcode within epsilon if it is an error sequence of it
    (:func:`is_error_sequence`). Otherwise it becomes a putative barcode, with
    its position in the sorted order as cluster label.

    Algorithm S1 looks up the k-mer neighbourhood of each sequence in an index
    of all sequences and keeps the putative barcodes in it. Here the index only
    contains the putative barcodes found so far, which yields that set directly
    and keeps the index small: since sequences are processed in descending count
    order, every putative barcode that a sequence can be merged with has already
    been added.
    """
    index = KmerIndex.for_parameters(params)
    clustering = Clustering()
    for label, seq in enumerate(sort_by_count(counts)):
        count = counts[seq]
        if count < params.count_threshold:
            match = find_closest_barcode(seq, index.neighbours(seq), counts, params.epsilon)
            if match is not None and is_error_sequence(
                count, counts[match.barcode], match.distance, params
            ):
                clustering.merge(seq, count, match.barcode)
                continue
        clustering.add_barcode(seq, count, label)
        index.add(seq)
    return clustering


def correct_indels(reads: ReadCounts, clustering: Clustering) -> None:
    """Merge sequences with a single insertion or deletion into putative barcodes (Section 2.2.1).

    A sequence of length l + 1 (l - 1) joins the first putative barcode that is
    obtained by deleting (inserting) one nucleotide. Sequences without such a
    barcode are left unclustered.
    """
    for seq, count in reads.insertions.items():
        deletions = (seq[:i] + seq[i + 1 :] for i in range(len(seq)))
        _merge_into_first_barcode(seq, count, deletions, clustering)
    for seq, count in reads.deletions.items():
        insertions = (seq[:i] + n + seq[i:] for i in range(len(seq) + 1) for n in 'ACGT')
        _merge_into_first_barcode(seq, count, insertions, clustering)


def _merge_into_first_barcode(
    seq: str, count: int, variants: Iterable[str], clustering: Clustering
) -> None:
    barcode = next((v for v in variants if v in clustering.barcode_counts), None)
    if barcode is not None:
        clustering.merge(seq, count, barcode)


def cluster(reads: ReadCounts, params: Parameters) -> Clustering:
    """Cluster the reads of one time point, including single insertion and deletion errors."""
    clustering = cluster_reads(reads.barcodes, params)
    correct_indels(reads, clustering)
    return clustering
