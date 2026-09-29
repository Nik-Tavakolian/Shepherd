"""Tests for following barcodes through later time points (shepherd track)."""

from conftest import read_multi_freqs, read_seq_clust, run_cluster, run_track, write_counts
from shepherd.io import ReadCounts
from shepherd.model import Parameters
from shepherd.tracking import Tracker


def test_pipeline_conserves_reads_of_stable_barcodes(tmp_path, background):
    write_counts(tmp_path / 't0.txt', background, {})
    write_counts(tmp_path / 't1.txt', background, {})
    run_cluster(tmp_path, 't0.txt')
    run_track(tmp_path, 't0.txt', ['t1.txt'])

    freqs = read_multi_freqs(tmp_path)
    total_t0 = sum(counts[0] for counts in freqs.values())
    total_t1 = sum(counts[1] for counts in freqs.values())
    assert total_t0 == total_t1


def test_later_time_points_find_barcodes_with_unshared_k_mer_combinations(tmp_path, background):
    # With k = 3 (as used for the Illumina data in the paper) most k-mer
    # combinations occur in a single t0 sequence. A barcode without error
    # reads at t0 therefore shares no combination with any other t0 sequence.
    barcode = 'TGCATGCATGCATGCATGCA'
    single_error = 'TGCATGCATGCATGCATGCT'
    write_counts(tmp_path / 't0.txt', background, {barcode: 20})
    write_counts(tmp_path / 't1.txt', background, {barcode: 20, single_error: 1})
    run_cluster(tmp_path, 't0.txt', '-k', '3')
    run_track(tmp_path, 't0.txt', ['t1.txt'])

    labels = read_seq_clust(tmp_path / 't1_seq_clust.csv')
    assert labels[single_error] == labels[barcode]
    assert read_multi_freqs(tmp_path)[barcode] == [20, 21]


def single_substitutions(seq):
    return [seq[:i] + n + seq[i + 1 :] for i in range(len(seq)) for n in 'ACGT' if n != seq[i]]


def test_error_reads_of_an_emerging_barcode_are_merged_into_it(tmp_path, background):
    emerging = 'ACGTTGCAACGTAGCTAGCA'
    errors = single_substitutions(emerging)[:30]
    write_counts(tmp_path / 't0.txt', background, {})
    write_counts(tmp_path / 't1.txt', background, {emerging: 300, **{seq: 2 for seq in errors}})
    write_counts(tmp_path / 't2.txt', background, {emerging: 600, **{seq: 2 for seq in errors}})
    run_cluster(tmp_path, 't0.txt')
    run_track(tmp_path, 't0.txt', ['t1.txt', 't2.txt'])

    freqs = read_multi_freqs(tmp_path)
    assert freqs[emerging] == [0, 360, 660]
    assert not set(errors) & set(freqs), 'error sequences reported as lineages'


def test_separating_an_emerging_barcode_counts_each_read_once(tmp_path, background):
    barcode = 'CCCCCCCCCCGGGGGGGGGG'
    emerging = 'CCCCCCCCCCGGGGGGGGTT'  # distance 2 from barcode
    near_emerging = 'CCCCCCCCCCGGGGGGGTTT'  # distance 3 from barcode, 1 from emerging
    write_counts(tmp_path / 't0.txt', background, {barcode: 1000})
    write_counts(tmp_path / 't1.txt', background, {barcode: 1000, emerging: 200, near_emerging: 50})
    run_cluster(tmp_path, 't0.txt')
    run_track(tmp_path, 't0.txt', ['t1.txt'])

    freqs = read_multi_freqs(tmp_path)
    assert freqs[barcode] == [1000, 1000]
    assert freqs[emerging] == [0, 250]
    assert near_emerging not in freqs


# Unit tests of the Tracker on hand-made data.
PARAMS = Parameters(
    barcode_length=20,
    error_rate=0.005,
    max_count=100,
    epsilon=3,
    tau=2,
    count_threshold=15,
    kmer_length=4,
    n_partitions=5,
)
BARCODE = 'AAAAAAAAAACCCCCCCCCC'
NEW = 'GGGGGGGGGGTTTTTTTTTT'


def test_reads_are_assigned_to_the_barcodes_of_the_previous_time_point():
    tracker = Tracker({BARCODE: 100}, PARAMS)
    error = 'AAAAAAAAAACCCCCCCCCA'
    clustering = tracker.add_time_point(ReadCounts(barcodes={BARCODE: 80, error: 3}))

    assert clustering.barcode_counts == {BARCODE: 83}
    assert clustering.labels[error] == clustering.labels[BARCODE]


def test_a_new_barcode_counts_from_the_time_point_after_it_is_first_seen():
    tracker = Tracker({BARCODE: 100}, PARAMS)
    new_error = 'GGGGGGGGGGTTTTTTTTTA'

    # t1: NEW is far from every barcode, so it and its error form a candidate.
    tracker.add_time_point(ReadCounts(barcodes={BARCODE: 90, NEW: 40, new_error: 2}))
    assert NEW not in tracker.counts_per_time_point[1]

    # t2: NEW is observed again, so it is confirmed, also at t1.
    tracker.add_time_point(ReadCounts(barcodes={BARCODE: 95, NEW: 60, new_error: 1}))
    assert [counts.get(NEW, 0) for counts in tracker.counts_per_time_point] == [0, 42, 61]
