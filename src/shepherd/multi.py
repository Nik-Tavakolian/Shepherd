import csv
import math
import pickle
import time

from scipy.stats import binom

from shepherd.kmer_index import add_seq_to_k_mer_dict, build_k_mer_dict, get_candidates
from shepherd.single import cluster_reads


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


def get_closest_pb(seq, pb_to_freq_dict_t0, k_mer_dict, q, l, p, eps):

    candidates = get_candidates(seq, k_mer_dict, q, l, p, eps)
    if candidates:
        pb_neighbors = [cand for cand in candidates if cand in pb_to_freq_dict_t0]
        if pb_neighbors:
            min_dist, indices = locate_mins(
                [trunc_ham_dist(seq, pb_neighbor, eps, l) for pb_neighbor in pb_neighbors]
            )
            if min_dist == l:
                return None
            if len(indices) == 1:
                closest_pb = pb_neighbors[indices[0]]
            else:
                # Ties: the higher count wins, then the alphabetically first sequence,
                # so the choice does not depend on the iteration order of a set.
                closest_pb = min(
                    [pb_neighbors[j] for j in indices], key=lambda x: (-pb_to_freq_dict_t0[x], x)
                )

            return closest_pb, min_dist


def get_log_K(f_n, f_c, p_no_err, d, l, total_err_rate, logdenom):

    n_hat = max(int(f_c / p_no_err), f_c + f_n)
    p_est = (total_err_rate / 3) ** d * (1 - total_err_rate) ** (l - d)

    return binom.logpmf(f_n, n_hat, p_est) + math.log(p_est) + logdenom


def classify_reads(
    seq_list_sorted, seq_freq_dict, k_mer_dict, pb_to_freq_dict_t0, unassigned_pb_freq_dict, params
):

    q, l, p, eps, p_no_err, total_err_rate, bft, logdenom, f, tau = params
    unassigned_seqs_dict = {}
    pb_to_freq_dict_t1 = {}
    seq_to_clust_dict_t1 = {}
    pb_to_seqs_dict_t1 = {}
    seq_to_dist_dict_t1 = {}
    for i, seq in enumerate(seq_list_sorted):
        freq = seq_freq_dict[seq]
        if seq in pb_to_freq_dict_t0:
            if seq in pb_to_freq_dict_t1:
                pb_to_freq_dict_t1[seq] += freq
                continue
            pb_to_freq_dict_t1[seq] = freq
            pb_to_seqs_dict_t1[seq] = []
            seq_to_clust_dict_t1[seq] = i
        elif seq in unassigned_pb_freq_dict:
            pb_to_freq_dict_t0[seq] = unassigned_pb_freq_dict[seq]
            pb_to_freq_dict_t1[seq] = freq
            pb_to_seqs_dict_t1[seq] = []
            seq_to_clust_dict_t1[seq] = i
            k_mer_dict = add_seq_to_k_mer_dict(seq, k_mer_dict, q, l, p, eps)
        else:
            out = get_closest_pb(seq, pb_to_freq_dict_t0, k_mer_dict, q, l, p, eps)
            if out:
                closest_pb, dist = out
                seq_to_dist_dict_t1[seq] = dist
                if closest_pb:
                    if closest_pb not in pb_to_freq_dict_t1:
                        pb_to_freq_dict_t1[closest_pb] = freq
                        pb_to_seqs_dict_t1[closest_pb] = [seq]
                        seq_to_clust_dict_t1[closest_pb] = i
                        seq_to_clust_dict_t1[seq] = i
                    else:
                        pb_to_freq_dict_t1[closest_pb] += freq
                        pb_to_seqs_dict_t1[closest_pb].append(seq)
                        seq_to_clust_dict_t1[seq] = seq_to_clust_dict_t1[closest_pb]
            else:
                unassigned_seqs_dict[seq] = freq

    return (
        pb_to_freq_dict_t0,
        pb_to_freq_dict_t1,
        seq_to_clust_dict_t1,
        pb_to_seqs_dict_t1,
        seq_to_dist_dict_t1,
        unassigned_seqs_dict,
        k_mer_dict,
    )


def separate_emerging(
    pb_to_freq_dict_t1,
    seq_freq_dict,
    seq_to_clust_dict_t1,
    pb_to_seqs_dict_t1,
    seq_to_dist_dict_t1,
    k_mer_dict,
    params,
):

    q, l, p, eps, p_no_err, total_err_rate, bft, logdenom, f, tau = params
    id_count = max(seq_to_clust_dict_t1.values())
    reassigned_seqs = set()
    for pb, seqs in pb_to_seqs_dict_t1.items():
        for i, seq in enumerate(seqs):
            if pb not in seq_freq_dict:
                break
            if seq in reassigned_seqs:
                continue
            f_c = seq_freq_dict[seq]
            f_b = seq_freq_dict[pb]
            if f_b < f_c:
                break
            d = seq_to_dist_dict_t1[seq]
            logK = get_log_K(f_c, f_b, p_no_err, d, l, total_err_rate, logdenom)
            if not logK > bft:
                id_count += 1
                pb_to_freq_dict_t1[seq] = f_c
                pb_to_freq_dict_t1[pb] -= f_c
                seq_to_clust_dict_t1[seq] = id_count
                k_mer_dict = add_seq_to_k_mer_dict(seq, k_mer_dict, q, l, p, eps)
                for s in seqs[i + 1 :]:
                    if s in reassigned_seqs:
                        continue
                    f_c_new = seq_freq_dict[s]
                    d = seq_to_dist_dict_t1[s]
                    d_new = trunc_ham_dist(seq, s, eps, l)
                    logK = get_log_K(f_c_new, f_b, p_no_err, d, l, total_err_rate, logdenom)
                    logK_new = get_log_K(f_c_new, f_c, p_no_err, d_new, l, total_err_rate, logdenom)
                    if logK_new > logK:
                        pb_to_freq_dict_t1[seq] += f_c_new
                        pb_to_freq_dict_t1[pb] -= f_c_new
                        seq_to_clust_dict_t1[s] = id_count
                        reassigned_seqs.add(s)

    return pb_to_freq_dict_t1, seq_to_clust_dict_t1, k_mer_dict


def classify_unassigned(
    unassigned_seq_dict, pb_to_freq_dict, seq_to_clust_dict, k_mer_dict, params
):

    unassigned_seq_dict_copy = unassigned_seq_dict.copy()
    for seq_u, f_u in unassigned_seq_dict_copy.items():
        out = get_closest_pb(
            seq_u, pb_to_freq_dict, k_mer_dict, params[0], params[1], params[2], params[3]
        )
        if out:
            closest_pb, _ = out
            if closest_pb:
                pb_to_freq_dict[closest_pb] += f_u
                seq_to_clust_dict[seq_u] = seq_to_clust_dict[closest_pb]
                del unassigned_seq_dict[seq_u]

    return unassigned_seq_dict, pb_to_freq_dict, seq_to_clust_dict


def cluster_unassigned(seq_list, seq_to_freq_dict, params):
    """Cluster the unassigned reads with the single time point procedure."""

    q, l, p, eps, p_no_err, total_err_rate, bft, logdenom, f, tau = params
    return cluster_reads(
        seq_list, seq_to_freq_dict, q, p, eps, tau, f, l, p_no_err, total_err_rate, logdenom, bft
    )


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
    """Run ``shepherd track`` with arguments parsed by shepherd.cli."""

    f0_prefix = args.f0[:-4]
    multi_filenames = args.fn

    if args.o is None:
        o_fn_prefix = 'multi_freqs'
    else:
        o_fn_prefix = args.o

    with open(f0_prefix + '_pb_freq.csv') as handle:
        pb_to_freq_dict_t0 = {}
        i = 0
        for line in handle:
            if i > 0:
                key, value = line.split(',')
                pb_to_freq_dict_t0[key] = int(value)
            i += 1

    with open(f0_prefix + '_params', 'rb') as f:
        params = pickle.load(f)

    k_mer_dict = build_k_mer_dict(pb_to_freq_dict_t0, params[0], params[1], params[2], params[3])

    print('Starting classification')
    print('\t')

    start = time.time()
    pb_to_freq_dict_list = []
    unassigned_pb_freq_dict = {}
    i = 1
    for filename in multi_filenames:
        print('Classifying time point ' + str(i))
        deletions_dict = {}
        insertions_dict = {}
        seq_freq_dict = {}
        l = params[1]
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

        seq_list = [
            seq for seq, freq in sorted(seq_freq_dict.items(), key=lambda x: x[1], reverse=True)
        ]
        (
            pb_to_freq_dict_t0,
            pb_to_freq_dict_t1,
            seq_to_clust_dict,
            pb_to_seqs_dict,
            seq_to_dist_dict,
            unassigned_seq_dict,
            k_mer_dict,
        ) = classify_reads(
            seq_list, seq_freq_dict, k_mer_dict, pb_to_freq_dict_t0, unassigned_pb_freq_dict, params
        )

        pb_to_freq_dict_t1, seq_to_clust_dict, k_mer_dict = separate_emerging(
            pb_to_freq_dict_t1,
            seq_freq_dict,
            seq_to_clust_dict,
            pb_to_seqs_dict,
            seq_to_dist_dict,
            k_mer_dict,
            params,
        )

        if unassigned_seq_dict:
            unassigned_seq_dict, pb_to_freq_dict_t1, seq_to_clust_dict = classify_unassigned(
                unassigned_seq_dict, pb_to_freq_dict_t1, seq_to_clust_dict, k_mer_dict, params
            )

            unassigned_seq_list = [
                seq
                for seq, freq in sorted(
                    unassigned_seq_dict.items(), key=lambda x: x[1], reverse=True
                )
            ]

            unassigned_seq_to_clust_dict, unassigned_pb_freq_dict = cluster_unassigned(
                unassigned_seq_list, unassigned_seq_dict, params
            )
        else:
            unassigned_pb_freq_dict = {}

        seq_to_clust_dict, pb_to_freq_dict = correct_insertions(
            insertions_dict, pb_to_freq_dict_t1, seq_to_clust_dict, l
        )
        seq_to_clust_dict, pb_to_freq_dict = correct_deletions(
            deletions_dict, pb_to_freq_dict_t1, seq_to_clust_dict, l
        )

        with open(filename[:-4] + '_seq_clust.csv', 'w', newline='') as res:
            writer = csv.writer(res)
            writer.writerow(['sequence', 'cluster'])
            for key, value in seq_to_clust_dict.items():
                writer.writerow([key, value])

        pb_to_freq_dict_list.append(pb_to_freq_dict_t0)
        pb_to_freq_dict_t0 = pb_to_freq_dict_t1.copy()
        print(
            'Time Point '
            + str(i)
            + ' was successfully classified. Results were saved to '
            + filename[:-4]
            + '_seq_clust.csv \n'
        )
        i += 1

    pb_to_freq_dict_list.append(pb_to_freq_dict_t1)
    end = time.time()
    print('\t')
    print('Classification time: ' + str(end - start))

    with open(o_fn_prefix + '.csv', 'w', newline='') as result:
        writer = csv.writer(result, delimiter=',')
        col_names = ['barcode']
        time_points = ['time_point_' + str(i) for i in range(1, len(pb_to_freq_dict_list) + 1)]
        col_names.extend(time_points)
        writer.writerow(col_names)
        pb_set = set([pb for pb_dict in pb_to_freq_dict_list for pb in pb_dict])
        for seq in pb_set:
            line = [seq]
            for pb_to_freq_dict in pb_to_freq_dict_list:
                if seq in pb_to_freq_dict:
                    line.append(pb_to_freq_dict[seq])
                else:
                    line.append(0)

            writer.writerow(line)

    print('Results were saved to ' + o_fn_prefix + '.csv')
