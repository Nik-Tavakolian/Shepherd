"""Tests for the Bayesian test (shepherd.model)."""

import math

import pytest
from scipy.stats import binom

from shepherd.model import log_binom_pmf


@pytest.mark.parametrize('n', [1, 10, 1_000, 1_000_000])
@pytest.mark.parametrize('p', [1e-12, 1e-6, 1e-3, 0.3])
def test_log_binom_pmf_matches_scipy(n, p):
    for k in {0, 1, min(2, n), n // 2, n}:
        assert log_binom_pmf(k, n, p) == pytest.approx(binom.logpmf(k, n, p), rel=1e-9, abs=1e-9)


def test_log_binom_pmf_is_normalised():
    assert math.fsum(math.exp(log_binom_pmf(k, 30, 0.2)) for k in range(31)) == pytest.approx(1)
