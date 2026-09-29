"""Tests for the Bayesian test (shepherd.model)."""

import math

import pytest
from scipy.stats import binom

from shepherd.model import ErrorModel, log_binom_pmf


@pytest.mark.parametrize('n', [1, 10, 1_000, 1_000_000])
@pytest.mark.parametrize('p', [1e-12, 1e-6, 1e-3, 0.3])
def test_log_binom_pmf_matches_scipy(n, p):
    for k in {0, 1, min(2, n), n // 2, n}:
        assert log_binom_pmf(k, n, p) == pytest.approx(binom.logpmf(k, n, p), rel=1e-9, abs=1e-9)


def test_log_binom_pmf_is_normalised():
    assert math.fsum(math.exp(log_binom_pmf(k, 30, 0.2)) for k in range(31)) == pytest.approx(1)


def test_conversion_probabilities_sum_to_one():
    # Eq. S4: summed over all 4^l sequences, of which C(l, d) 3^d lie at distance d.
    model = ErrorModel(barcode_length=20, error_rate=0.01, max_count=1000)
    total = math.fsum(math.comb(20, d) * 3**d * model.conversion_probability(d) for d in range(21))
    assert total == pytest.approx(1)


def test_p_no_error():
    model = ErrorModel(barcode_length=20, error_rate=0.01, max_count=1000)
    assert model.p_no_error == pytest.approx(0.99**20)
    assert model.conversion_probability(0) == pytest.approx(model.p_no_error)


def test_log_bayes_factor_favours_true_barcodes_further_away():
    model = ErrorModel(barcode_length=20, error_rate=0.01, max_count=1000)
    log_k = [model.log_bayes_factor(count=2, parent_count=500, distance=d) for d in range(1, 6)]
    assert log_k == sorted(log_k, reverse=True)
    assert log_k[0] > 0 > log_k[-1]


def test_log_bayes_factor_favours_true_barcodes_with_high_counts():
    # 15 vs 14 reads at distance 1 are two barcodes, not a barcode and its error.
    model = ErrorModel(barcode_length=20, error_rate=0.005, max_count=10_000)
    assert model.log_bayes_factor(count=1, parent_count=15, distance=1) > 0
    assert model.log_bayes_factor(count=14, parent_count=15, distance=1) < -4
