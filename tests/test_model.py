"""Tests for the error model and parameter selection (shepherd.model)."""

import math

import pytest
from scipy.stats import binom

from shepherd.model import (
    ErrorModel,
    Parameters,
    ShepherdError,
    choose_kmer_length,
    estimate_error_rate,
    estimate_parameters,
    log_binom_pmf,
)


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


@pytest.mark.parametrize(
    ('epsilon', 'expected'),
    [(3, (4, 5)), (4, (3, 7))],  # Supplementary Table S1: data sets A/B and C
)
def test_kmer_length_matches_the_paper(epsilon, expected):
    assert choose_kmer_length(20, epsilon) == expected


def test_error_rate_estimate_is_close_to_the_simulated_rate(background):
    estimate = estimate_error_rate(background, 20, n_top=500)
    assert estimate == pytest.approx(0.005, rel=0.05)


def test_given_values_are_used_and_the_rest_estimated(background):
    params = estimate_parameters(background, 20, error_rate=0.01, kmer_length=3, tau=1)
    assert (params.error_rate, params.kmer_length, params.n_partitions, params.tau) == (
        0.01,
        3,
        7,
        1,
    )
    assert params.epsilon == 3
    assert params.max_count == max(background.values())


def test_no_sequences_of_the_barcode_length():
    with pytest.raises(ShepherdError, match='no sequences of length 20'):
        estimate_parameters({}, 20)


def test_parameters_round_trip_through_json(tmp_path):
    params = Parameters(20, 0.005, 600, 3, 2, 17, 4, 5)
    params.save(tmp_path / 'params.json')
    assert Parameters.load(tmp_path / 'params.json') == params


def test_parameter_files_from_version_1_are_reported(tmp_path):
    (tmp_path / 't0_params').write_bytes(b'pickled list from Shepherd 1.x')
    with pytest.raises(ShepherdError, match=r'Shepherd 1\.x'):
        Parameters.load(tmp_path / 't0_params.json')
