import csv
import time

from shepherd.kmer_index import KmerIndex
from shepherd.parameters import DEFAULT_LOG_BF_THRESHOLD, DEFAULT_N_TOP, estimate_parameters
from shepherd.sequences import sort_by_count


def trunc_ham_dist(seq_1, seq_2, d, n):

    h = 0
    for n_1, n_2 in zip(seq_1, seq_2):
        if n_1 != n_2:
            h += 1
        if h > d:
            return n

    return h


def locate_mins(a):

    smallest = min(a)
    return smallest, [index for index, element in enumerate(a) if smallest == element]


def cluster_reads(seq_list, seq_to_freq_dict, params):

    l, eps, tau, f = params.barcode_length, params.epsilon, params.tau, params.count_threshold
    model = params.error_model
    # Only putative barcodes can absorb a sequence, so the k-mer Index holds just the
    # putative barcodes found so far. Sequences are processed in descending count order
    # and added to the index when they are classified as putative barcodes.
    index = KmerIndex(l, params.kmer_length, params.n_partitions, eps)
    pb_to_freq_dict = {}
    seq_to_clust_dict = {}
    for i, S_c in enumerate(seq_list):
        f_c = seq_to_freq_dict[S_c]
        if f_c < f:
            pb_neighbors = list(index.neighbours(S_c))
            if pb_neighbors:
                min_dist, indices = locate_mins(
                    [trunc_ham_dist(S_c, pb_neighbor, eps, l) for pb_neighbor in pb_neighbors]
                )
                if len(indices) == 1:
                    S_b = pb_neighbors[indices[0]]
                else:
                    # Ties: the higher count wins, then the alphabetically first sequence,
                    # so the choice does not depend on the iteration order of a set.
                    S_b = min(
                        [pb_neighbors[j] for j in indices],
                        key=lambda x: (-seq_to_freq_dict[x], x),
                    )
                if min_dist != l:
                    if (f_c == 1 and min_dist <= tau) or min_dist == 1:
                        seq_to_clust_dict[S_c] = seq_to_clust_dict[S_b]
                        pb_to_freq_dict[S_b] += f_c
                        continue

                    f_b = seq_to_freq_dict[S_b]
                    logK = model.log_bayes_factor(f_c, f_b, min_dist)
                    if logK > params.log_bf_threshold:
                        seq_to_clust_dict[S_c] = seq_to_clust_dict[S_b]
                        pb_to_freq_dict[S_b] += f_c
                        continue

        seq_to_clust_dict[S_c] = i
        pb_to_freq_dict[S_c] = f_c
        index.add(S_c)

    return seq_to_clust_dict, pb_to_freq_dict


def correct_deletions(deletions_dict, pb_to_freq_dict, seq_to_clust_dict, l):
    for seq, freq in deletions_dict.items():
        for i in range(l):
            for n in ['A', 'C', 'G', 'T']:
                corrected_seq = seq[:i] + n + seq[i:]
                if corrected_seq in pb_to_freq_dict:
                    pb_to_freq_dict[corrected_seq] += freq
                    seq_to_clust_dict[seq] = seq_to_clust_dict[corrected_seq]
                    break
            else:
                continue
            break

    return seq_to_clust_dict, pb_to_freq_dict


def correct_insertions(insertions_dict, pb_to_freq_dict, seq_to_clust_dict, l):

    for seq, freq in insertions_dict.items():
        for i in range(l + 1):
            corrected_seq = seq[:i] + seq[i + 1 :]
            if corrected_seq in pb_to_freq_dict:
                pb_to_freq_dict[corrected_seq] += freq
                seq_to_clust_dict[seq] = seq_to_clust_dict[corrected_seq]
                break

    return seq_to_clust_dict, pb_to_freq_dict


def run(args):
    """Run ``shepherd cluster`` with arguments parsed by shepherd.cli."""

    filename = args.f
    file_prefix = filename[:-4]
    l = args.l

    deletions_dict = {}
    insertions_dict = {}
    seq_freq_dict = {}
    with open(filename) as a_file:
        for line in a_file:
            seq, freq = line.split()
            seq_len = len(seq)
            if seq_len == l:
                seq_freq_dict[seq] = int(freq)
            if seq_len == l - 1:
                deletions_dict[seq] = int(freq)
            if seq_len == l + 1:
                insertions_dict[seq] = int(freq)

    start_tot = time.time()

    seq_list = sort_by_count(seq_freq_dict)
    params = estimate_parameters(
        seq_freq_dict,
        l,
        error_rate=args.e,
        epsilon=args.eps,
        kmer_length=args.k,
        tau=args.tau,
        count_threshold=args.ft,
        log_bf_threshold=DEFAULT_LOG_BF_THRESHOLD if args.bft is None else args.bft,
        n_top=DEFAULT_N_TOP if args.Nh is None else args.Nh,
    )

    print('Shepherd Single Parameters:')
    print('\t')
    print('Sequence length: ' + str(params.barcode_length))
    print('Substitution Error Rate Estimate: ' + str(params.error_rate))
    print('epsilon: ' + str(params.epsilon))
    print('tau: ' + str(params.tau))
    print('f: ' + str(params.count_threshold))
    print('Bayes Factor Threshold: ' + str(params.log_bf_threshold))
    print('Substring Length: ' + str(params.kmer_length))
    print('Number of Partitions: ' + str(params.n_partitions))
    print('\t')

    start = time.time()
    seq_to_clust_dict, pb_to_freq_dict = cluster_reads(seq_list, seq_freq_dict, params)
    seq_to_clust_dict, pb_to_freq_dict = correct_insertions(
        insertions_dict, pb_to_freq_dict, seq_to_clust_dict, l
    )
    seq_to_clust_dict, pb_to_freq_dict = correct_deletions(
        deletions_dict, pb_to_freq_dict, seq_to_clust_dict, l
    )
    end = time.time()
    print('Clustering time: ' + str(end - start))

    params.save(file_prefix + '_params.json')

    with open(file_prefix + '_seq_clust.csv', 'w', newline='') as fh:
        writer = csv.writer(fh)
        i = 0
        for key, value in seq_to_clust_dict.items():
            if i == 0:
                writer.writerow(['sequence', 'cluster'])
            writer.writerow([key, value])
            i += 1

    with open(file_prefix + '_pb_freq.csv', 'w', newline='') as fh:
        writer = csv.writer(fh)
        i = 0
        for key, value in pb_to_freq_dict.items():
            if i == 0:
                writer.writerow(['barcode', 'frequency'])
            writer.writerow([key, value])
            i += 1

    end_tot = time.time()
    print('Total time: ' + str((end_tot - start_tot) / 60))
