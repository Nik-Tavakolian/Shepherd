"""Tests for parameter selection (Supplementary Section 3)."""

import pytest

from shepherd.errors import ShepherdError
from shepherd.parameters import (
    Parameters,
    choose_kmer_length,
    estimate_error_rate,
    estimate_parameters,
)
from shepherd.sequences import sort_by_count


@pytest.mark.parametrize(
    ('epsilon', 'expected'),
    [(3, (4, 5)), (4, (3, 7))],  # Supplementary Table S1: data sets A/B and C
)
def test_kmer_length_matches_the_paper(epsilon, expected):
    assert choose_kmer_length(20, epsilon) == expected


def test_error_rate_estimate_is_close_to_the_simulated_rate(background):
    estimate = estimate_error_rate(background, sort_by_count(background), 20, n_top=500)
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
