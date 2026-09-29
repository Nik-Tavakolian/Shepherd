"""Regression tests for multiple time point mode (shepherd_multi.py)."""
from conftest import read_multi_freqs, run_multi, run_t0, write_counts
from shepherd_multi import add_seq_to_k_mer_dict


def test_pipeline_conserves_reads_of_stable_barcodes(tmp_path, background):
    write_counts(tmp_path / 't0.txt', background, {})
    write_counts(tmp_path / 't1.txt', background, {})
    run_t0(tmp_path, 't0.txt')
    run_multi(tmp_path, 't0.txt', ['t1.txt'])

    freqs = read_multi_freqs(tmp_path)
    total_t0 = sum(counts[0] for counts in freqs.values())
    total_t1 = sum(counts[1] for counts in freqs.values())
    assert total_t0 == total_t1


def test_add_seq_to_k_mer_dict_stores_the_sequence_under_new_keys():
    seq = 'ACGTACGTACGTACGTACGT'
    k_mer_dict = add_seq_to_k_mer_dict(seq, {}, q=4, l=20, p=5, eps=3)

    assert len(k_mer_dict) == 10  # 5 choose 2 k-mer combinations
    assert all(seqs == {seq} for seqs in k_mer_dict.values())
