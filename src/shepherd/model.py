"""The substitution error model behind the Bayesian test (Section 2.2.2, Supplementary Section 1).

Each nucleotide of a barcode is misread independently with probability rho
(the error rate) and replaced by each of the other three nucleotides with equal
probability. The test compares two models for a sequence S_c with read count
f_c that lies at Hamming distance d from a putative barcode with read count f_p:

    M1: S_c is an error sequence of the putative barcode,
    M2: S_c is a true barcode.
"""

import math
from dataclasses import dataclass
from functools import cached_property

from scipy.stats import binom


def log_binom_pmf(k: int, n: int, p: float) -> float:
    """Log of the binomial probability mass function, ln C(n, k) p^k (1 - p)^(n - k).

    The same formula as scipy.stats.binom.logpmf, without SciPy's per-call
    overhead (about 60 times faster for scalar arguments).
    """
    log_comb = math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    return log_comb + k * math.log(p) + (n - k) * math.log1p(-p)


@dataclass(frozen=True)
class ErrorModel:
    """Substitution error model for barcodes of a fixed length.

    ``max_count`` is the highest read count f_max in the data, which bounds the
    uniform prior on the read count of a true barcode (Eq. S10).
    """

    barcode_length: int
    error_rate: float
    max_count: int

    @cached_property
    def p_no_error(self) -> float:
        """Probability that a barcode is read without errors, (1 - rho)^l (Eq. S13)."""
        return float(binom.pmf(0, self.barcode_length, self.error_rate))

    @cached_property
    def neg_log_likelihood_m2(self) -> float:
        """-ln P(f_c, S_c | M2) = l ln 4 + ln f_max (Eq. S11)."""
        return self.barcode_length * math.log(4) + math.log(self.max_count)

    def conversion_probability(self, distance: int) -> float:
        """Probability that a read of a barcode is one given sequence at this distance (Eq. S3)."""
        rho = self.error_rate
        return (rho / 3) ** distance * (1 - rho) ** (self.barcode_length - distance)

    def estimated_true_count(self, parent_count: int, count: int) -> int:
        """Estimate n of the number of reads of a barcode with parent_count error-free reads.

        The maximum likelihood estimate (Eq. S5), but at least parent_count + count,
        since both sequences originate from the same barcode under M1 (Eq. S6).
        """
        return max(int(parent_count / self.p_no_error), parent_count + count)

    def most_likely_error_count(self, parent_count: int, distance: int) -> int:
        """Most likely read count of an error sequence at this distance (Eq. S18, no max)."""
        n_mle = int(parent_count / self.p_no_error)
        return int(self.conversion_probability(distance) * (n_mle + 1))

    def log_bayes_factor(self, count: int, parent_count: int, distance: int) -> float:
        """ln K = ln P(f_c, S_c | M1) - ln P(f_c, S_c | M2) (Eq. 3 and Eq. S12).

        Positive values favour M1 (error sequence), negative values M2 (true barcode).
        """
        p = self.conversion_probability(distance)
        n = self.estimated_true_count(parent_count, count)
        return log_binom_pmf(count, n, p) + math.log(p) + self.neg_log_likelihood_m2
