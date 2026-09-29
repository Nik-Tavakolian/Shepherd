"""Regression tests for multiple time point mode (shepherd_multi.py)."""
from conftest import read_multi_freqs, read_seq_clust, run_multi, run_t0, write_counts
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


def test_later_time_points_find_barcodes_with_unshared_k_mer_combinations(tmp_path, background):
    # With k = 3 (as used for the Illumina data in the paper) most k-mer
    # combinations occur in a single t0 sequence. A barcode without error
    # reads at t0 therefore shares no combination with any other t0 sequence.
    barcode = 'TGCATGCATGCATGCATGCA'
    single_error = 'TGCATGCATGCATGCATGCT'
    write_counts(tmp_path / 't0.txt', background, {barcode: 20})
    write_counts(tmp_path / 't1.txt', background, {barcode: 20, single_error: 1})
    run_t0(tmp_path, 't0.txt', '-k', '3')
    run_multi(tmp_path, 't0.txt', ['t1.txt'])

    labels = read_seq_clust(tmp_path / 't1_seq_clust.csv')
    assert labels[single_error] == labels[barcode]
    assert read_multi_freqs(tmp_path)[barcode] == [20, 21]


def single_substitutions(seq):
    return [seq[:i] + n + seq[i + 1:] for i in range(len(seq)) for n in 'ACGT' if n != seq[i]]


def test_error_reads_of_an_emerging_barcode_are_merged_into_it(tmp_path, background):
    emerging = 'ACGTTGCAACGTAGCTAGCA'
    errors = single_substitutions(emerging)[:30]
    write_counts(tmp_path / 't0.txt', background, {})
    write_counts(tmp_path / 't1.txt', background, {emerging: 300, **{seq: 2 for seq in errors}})
    write_counts(tmp_path / 't2.txt', background, {emerging: 600, **{seq: 2 for seq in errors}})
    run_t0(tmp_path, 't0.txt')
    run_multi(tmp_path, 't0.txt', ['t1.txt', 't2.txt'])

    freqs = read_multi_freqs(tmp_path)
    assert freqs[emerging] == [0, 360, 660]
    assert not set(errors) & set(freqs), 'error sequences reported as lineages'
