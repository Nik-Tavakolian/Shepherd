"""Following the putative barcodes through later time points (Section 2.3, Suppl. Section 2)."""

from collections.abc import Mapping
from dataclasses import dataclass, field

from shepherd.clustering import (
    Clustering,
    cluster_reads,
    correct_indels,
    find_closest_barcode,
    sort_by_count,
    truncated_hamming_distance,
)
from shepherd.io import ReadCounts
from shepherd.kmer_index import KmerIndex
from shepherd.model import Parameters


@dataclass
class _TimePoint:
    """Working state while one time point is processed."""

    counts: Mapping[str, int]
    """Read counts of the sequences of the barcode length."""
    clustering: Clustering = field(default_factory=Clustering)
    members: dict[str, list[str]] = field(default_factory=dict)
    """Sequences assigned to each putative barcode in step 1, in order of decreasing count."""
    distances: dict[str, int] = field(default_factory=dict)
    """Distance of each sequence assigned in step 1 to its putative barcode."""
    unassigned: dict[str, int] = field(default_factory=dict)
    """Sequences without a putative barcode within epsilon, with their read counts."""


class Tracker:
    """Follows the putative barcodes of the first time point through later time points.

    Instead of clustering every time point anew, the reads of each time point
    are assigned to the putative barcodes of the previous one (Supplementary
    Section 2.B):

    1. Assign each sequence to the closest putative barcode within epsilon.
    2. Separate emerging barcodes from the barcodes they were assigned to,
       using the Bayesian test.
    3. Move the other sequences of such a cluster to the emerging barcode if
       they are more likely to originate from it.
    4. Assign sequences left over in step 1 to barcodes added at this time
       point.
    5. Cluster the remaining sequences. The resulting candidate barcodes are
       kept only if they are observed at the next time point.

    Each time point is added with :meth:`add_time_point`, and
    :attr:`counts_per_time_point` holds the barcode counts of all time points,
    starting with the first.
    """

    def __init__(self, first_counts: Mapping[str, int], params: Parameters) -> None:
        """``first_counts`` are the putative barcodes of the first time point and their counts."""
        self.params = params
        self.counts_per_time_point: list[dict[str, int]] = [dict(first_counts)]
        self._index = KmerIndex.for_parameters(params)
        for barcode in first_counts:
            self._index.add(barcode)
        self._candidates: dict[str, int] = {}

    def add_time_point(self, reads: ReadCounts) -> Clustering:
        """Assign the reads of the next time point to putative barcodes (steps 1 to 5)."""
        time_point = _TimePoint(reads.barcodes)
        self._assign_to_previous_barcodes(time_point)
        self._separate_emerging_barcodes(time_point)
        self._assign_to_new_barcodes(time_point)
        self._candidates = cluster_reads(time_point.unassigned, self.params).barcode_counts
        correct_indels(reads, time_point.clustering)
        self.counts_per_time_point.append(time_point.clustering.barcode_counts)
        return time_point.clustering

    def _assign_to_previous_barcodes(self, time_point: _TimePoint) -> None:
        """Step 1, and confirmation of the previous time point's candidate barcodes."""
        previous = self.counts_per_time_point[-1]
        clustering = time_point.clustering
        for label, seq in enumerate(sort_by_count(time_point.counts)):
            count = time_point.counts[seq]
            if seq in previous:
                if seq in clustering.barcode_counts:
                    # Some of its error sequences had more reads and came first.
                    clustering.barcode_counts[seq] += count
                else:
                    clustering.add_barcode(seq, count, label)
                    time_point.members[seq] = []
            elif seq in self._candidates:
                # A candidate barcode from step 5 of the previous time point is
                # observed again. It is also reported at the previous time point,
                # with the read count of its cluster there.
                previous[seq] = self._candidates[seq]
                clustering.add_barcode(seq, count, label)
                time_point.members[seq] = []
                self._index.add(seq)
            else:
                neighbours = (b for b in self._index.neighbours(seq) if b in previous)
                match = find_closest_barcode(seq, neighbours, previous, self.params.epsilon)
                if match is None:
                    time_point.unassigned[seq] = count
                    continue
                time_point.distances[seq] = match.distance
                if match.barcode in clustering.barcode_counts:
                    clustering.merge(seq, count, match.barcode)
                    time_point.members[match.barcode].append(seq)
                else:
                    # The barcode has not been seen at this time point yet: its
                    # cluster gets the label of this sequence.
                    clustering.add_barcode(match.barcode, count, label)
                    clustering.labels[seq] = label
                    time_point.members[match.barcode] = [seq]

    def _separate_emerging_barcodes(self, time_point: _TimePoint) -> None:
        """Step 2: sequences that fail the Bayesian test become putative barcodes."""
        model = self.params.error_model
        threshold = self.params.log_bf_threshold
        label = max(time_point.clustering.labels.values(), default=0)
        moved: set[str] = set()
        for barcode, members in time_point.members.items():
            if barcode not in time_point.counts:
                continue  # the barcode itself has no reads at this time point
            barcode_count = time_point.counts[barcode]
            for position, seq in enumerate(members):
                if seq in moved:
                    continue
                count = time_point.counts[seq]
                if count > barcode_count:
                    # As in Shepherd 1.x: since members are sorted by count, a
                    # cluster is not split if any sequence has more reads than
                    # the barcode itself. See "Implementation notes" in the README.
                    break
                log_k = model.log_bayes_factor(count, barcode_count, time_point.distances[seq])
                if log_k > threshold:
                    continue  # an error sequence of the barcode
                label += 1
                time_point.clustering.split_off(seq, count, barcode, label)
                self._index.add(seq)
                self._move_to_emerging_barcode(
                    seq, barcode, members[position + 1 :], moved, time_point
                )

    def _move_to_emerging_barcode(
        self,
        new_barcode: str,
        barcode: str,
        others: list[str],
        moved: set[str],
        time_point: _TimePoint,
    ) -> None:
        """Step 3: move the sequences that are more likely error sequences of the new barcode."""
        model = self.params.error_model
        counts = time_point.counts
        for other in others:
            if other in moved:
                continue
            distance = truncated_hamming_distance(new_barcode, other, self.params.epsilon)
            log_k_new = model.log_bayes_factor(counts[other], counts[new_barcode], distance)
            log_k = model.log_bayes_factor(
                counts[other], counts[barcode], time_point.distances[other]
            )
            if log_k_new > log_k:
                time_point.clustering.move(other, counts[other], barcode, new_barcode)
                moved.add(other)

    def _assign_to_new_barcodes(self, time_point: _TimePoint) -> None:
        """Step 4: assign unassigned sequences to barcodes of this time point, if within epsilon."""
        current = time_point.clustering.barcode_counts
        for seq, count in list(time_point.unassigned.items()):
            neighbours = (b for b in self._index.neighbours(seq) if b in current)
            match = find_closest_barcode(seq, neighbours, current, self.params.epsilon)
            if match is not None:
                time_point.clustering.merge(seq, count, match.barcode)
                del time_point.unassigned[seq]
