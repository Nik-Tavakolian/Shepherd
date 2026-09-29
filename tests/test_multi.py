"""Regression tests for multiple time point mode (shepherd_multi.py)."""
from conftest import read_multi_freqs, run_multi, run_t0, write_counts


def test_pipeline_conserves_reads_of_stable_barcodes(tmp_path, background):
    write_counts(tmp_path / 't0.txt', background, {})
    write_counts(tmp_path / 't1.txt', background, {})
    run_t0(tmp_path, 't0.txt')
    run_multi(tmp_path, 't0.txt', ['t1.txt'])

    freqs = read_multi_freqs(tmp_path)
    total_t0 = sum(counts[0] for counts in freqs.values())
    total_t1 = sum(counts[1] for counts in freqs.values())
    assert total_t0 == total_t1
