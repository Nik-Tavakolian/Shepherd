"""Shepherd's parameters and how they are chosen from the data (Supplementary Section 3)."""

import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from functools import cached_property
from pathlib import Path

from scipy.stats import binom

from shepherd.errors import ShepherdError
from shepherd.model import ErrorModel
from shepherd.sequences import single_substitutions, sort_by_count

DEFAULT_LOG_BF_THRESHOLD = -4.0
DEFAULT_N_TOP = 500
MAX_ESTIMATED_ERROR_RATE = 0.1


@dataclass(frozen=True)
class Parameters:
    """The parameters of a Shepherd run, fixed at the first time point."""

    barcode_length: int
    """Length l of the barcodes."""
    error_rate: float
    """Substitution error rate rho per nucleotide."""
    max_count: int
    """Highest read count f_max among the barcode-length sequences of the first time point."""
    epsilon: int
    """Largest Hamming distance at which a sequence can be merged with a putative barcode."""
    tau: int
    """One-read sequences within this distance of a putative barcode are merged without the test."""
    count_threshold: int
    """Sequences with at least this many reads (f_t) are putative barcodes without the test."""
    kmer_length: int
    """Length k of the k-mers in the k-mer Index."""
    n_partitions: int
    """Number of k-mers p per sequence, used for the size p - epsilon of the k-mer combinations."""
    log_bf_threshold: float = DEFAULT_LOG_BF_THRESHOLD
    """Sequences with ln K above this threshold are error sequences."""

    @cached_property
    def error_model(self) -> ErrorModel:
        return ErrorModel(self.barcode_length, self.error_rate, self.max_count)

    def save(self, path: str | Path) -> None:
        with open(path, 'w') as fh:
            json.dump(asdict(self), fh, indent=2)
            fh.write('\n')

    @classmethod
    def load(cls, path: str | Path) -> 'Parameters':
        path = Path(path)
        if not path.exists() and path.with_suffix('').exists():
            raise ShepherdError(
                f'{path} not found, but {path.with_suffix("")} from Shepherd 1.x is. '
                'Run `shepherd cluster` on the first time point again.'
            )
        with open(path) as fh:
            return cls(**json.load(fh))


def estimate_parameters(
    counts: Mapping[str, int],
    barcode_length: int,
    *,
    error_rate: float | None = None,
    epsilon: int | None = None,
    kmer_length: int | None = None,
    tau: int | None = None,
    count_threshold: int | None = None,
    log_bf_threshold: float = DEFAULT_LOG_BF_THRESHOLD,
    n_top: int = DEFAULT_N_TOP,
) -> Parameters:
    """Choose every parameter that is not given from the read counts of the first time point.

    ``counts`` maps each unique sequence of length ``barcode_length`` to its
    read count. ``n_top`` is the number of high-count sequences used to
    estimate the error rate.
    """
    if not counts:
        raise ShepherdError(f'the input contains no sequences of length {barcode_length}')

    if error_rate is None:
        error_rate = estimate_error_rate(counts, sort_by_count(counts), barcode_length, n_top)
        if error_rate == 0 or error_rate > MAX_ESTIMATED_ERROR_RATE:
            raise ShepherdError(
                f'the error rate could not be reliably estimated from the data (estimate: '
                f'{error_rate:.4g}). Please provide an error rate estimate with -e.'
            )

    model = ErrorModel(barcode_length, error_rate, max(counts.values()))
    if epsilon is None:
        epsilon = choose_epsilon(model)
    if tau is None:
        tau = choose_tau(model)
    if count_threshold is None:
        count_threshold = choose_count_threshold(model)
    if kmer_length is None:
        kmer_length, n_partitions = choose_kmer_length(barcode_length, epsilon)
    else:
        n_partitions = math.ceil(barcode_length / kmer_length)

    return Parameters(
        barcode_length=barcode_length,
        error_rate=error_rate,
        max_count=model.max_count,
        epsilon=epsilon,
        tau=tau,
        count_threshold=count_threshold,
        kmer_length=kmer_length,
        n_partitions=n_partitions,
        log_bf_threshold=log_bf_threshold,
    )


def estimate_error_rate(
    counts: Mapping[str, int], sorted_seqs: Sequence[str], barcode_length: int, n_top: int
) -> float:
    """Estimate the substitution error rate rho (Supplementary Section 3.A).

    The highest-count sequences are almost certainly true barcodes. With n0 their
    combined read count and n1 the combined read count of their single-nucleotide
    variants, rho = n1 / (n1 + l n0) (Eq. S17).
    """
    # The n_top + 1 highest-count sequences are used, as in Shepherd 1.x.
    top = sorted_seqs[: n_top + 1]
    n0 = sum(counts[seq] for seq in top)
    n1 = sum(counts.get(variant, 0) for seq in top for variant in single_substitutions(seq))
    ratio = n1 / (n0 * barcode_length)
    return ratio / (ratio + 1)


# In the three functions below a sequence counts as a true barcode when ln K < 0,
# i.e. when M2 is more likely than M1 (Supplementary Section 3.D), independently
# of the log Bayes factor threshold used during clustering.


def choose_epsilon(model: ErrorModel) -> int:
    """The largest distance at which a sequence can be an error sequence (Suppl. Section 3.B).

    For increasing distances, test the most likely read count (Eq. S18) of an
    error sequence of the highest-count barcode.
    """
    for distance in range(1, model.barcode_length):
        count = max(model.most_likely_error_count(model.max_count, distance), 1)
        if model.log_bayes_factor(count, model.max_count, distance) < 0:
            return distance - 1
    raise ShepherdError('epsilon could not be determined; please set it with -eps')


def choose_tau(model: ErrorModel) -> int:
    """The largest distance at which a one-read sequence is an error sequence of any barcode.

    The worst case is a putative barcode with a single read (Supplementary Section 3.D).
    """
    for distance in range(1, model.barcode_length):
        if model.log_bayes_factor(1, 1, distance) < 0:
            return distance - 1
    raise ShepherdError('tau could not be determined; please set it with -tau')


def choose_count_threshold(model: ErrorModel) -> int:
    """The smallest read count at which a sequence is a true barcode even at distance 1.

    Tested against the highest-count barcode, starting from the most likely
    error count (Supplementary Section 3.D).
    """
    start = model.most_likely_error_count(model.max_count, 1)
    for count in range(start, model.max_count):
        if model.log_bayes_factor(count, model.max_count, 1) < 0:
            return count
    raise ShepherdError('the count threshold could not be determined; please set it with -ft')


def choose_kmer_length(barcode_length: int, epsilon: int) -> tuple[int, int]:
    """Choose the k-mer length k and the number of partitions p (Supplementary Section 3.C).

    Returns the largest k (fewest partitions, each more than epsilon) for which a
    random barcode is more likely than not to lie either within epsilon or beyond
    the largest distance in the k-mer neighbourhood (Eq. S19).
    """
    length = barcode_length
    for n_partitions in range(epsilon + 1, length):
        # As in Shepherd 1.x, k is rounded from length / n_partitions, so the returned
        # n_partitions can differ from the actual number of k-mers, ceil(length / k).
        # For barcodes of 14 nt or more (checked up to 60 nt and epsilon <= 5) it is
        # never larger, so the index is still complete, only bigger than needed.
        k = round(length / n_partitions)
        remainder = length % k
        if remainder == 0:
            max_neighbourhood_distance = k * epsilon
        else:
            max_neighbourhood_distance = length - (k * (n_partitions - epsilon - 1) + remainder)
        p_in_between = sum(
            binom.pmf(d, length, 3 / 4) for d in range(epsilon + 1, max_neighbourhood_distance + 1)
        )
        if p_in_between < 0.5:
            return k, n_partitions
    raise ShepherdError('the k-mer length could not be determined; please set it with -k')
