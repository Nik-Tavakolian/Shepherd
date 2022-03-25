from itertools import combinations
import pickle
import json
import time
import math
import sys
import argparse
import csv
from scipy.stats import binom

def trunc_ham_dist(seq_1, seq_2, d, n):

    h = 0
    for n_1, n_2 in zip(seq_1, seq_2):
        if n_1 != n_2:
            h += 1
        if h > d:
            return n

    return h

def find_eps(l, highest_freq, p_no_err, total_err_rate, logdenom):

    f_c = highest_freq
    n_hat = int(f_c / p_no_err)
    for d in range(1, l):
        p_est = (total_err_rate / 3) ** d * (1 - total_err_rate) ** (l - d)
        f_n = int(p_est*(n_hat + 1))
        if f_n == 0:
            f_n = 1
            for dist in range(d, l):
                logK = get_log_K(f_n, f_c, p_no_err, d, l, total_err_rate, logdenom)
                if logK < 0:
                    return d - 1

        logK = get_log_K(f_n, f_c, p_no_err, d, l, total_err_rate, logdenom)
        if logK < 0:
            return d - 1

def estimate_rho(sorted_seq_list, seq_freq_dict_dict, l, N_h):
    correct_freq = 0
    total_SNP_count = 0
    for i, high_freq_seq in enumerate(sorted_seq_list):
        correct_freq += seq_freq_dict_dict[high_freq_seq]
        total_SNP_count += get_SNP_freq_sum(high_freq_seq, seq_freq_dict_dict)
        if i == N_h:
            break
        i += 1
    v = total_SNP_count / (correct_freq * l)
    rho = v / (v + 1)

    return rho


def get_SNP_freq_sum(high_freq_seq, seq_freq_dict_dict):

    SNP_count = 0
    nucleotides = {'A', 'C', 'T', 'G'}
    for i, nuc_current in enumerate(high_freq_seq):
        for nuc_new in nucleotides:
            if nuc_new != nuc_current:
                snp_seq = high_freq_seq[:i] + nuc_new + high_freq_seq[i + 1:]
                if snp_seq in seq_freq_dict_dict:
                    SNP_count += seq_freq_dict_dict[snp_seq]

    return SNP_count

def find_tau(l, p_no_err, total_err_rate, logdenom):

    f_c = 1
    f_n = 1
    for d in range(1, l):
        logK = get_log_K(f_n, f_c, p_no_err, d, l, total_err_rate, logdenom)
        if logK < 0:
            return d - 1

def find_f(l, highest_freq, p_no_err, total_err_rate, logdenom):

    f_c = highest_freq
    d = 1
    n_hat = int(f_c / p_no_err)
    p_est = (total_err_rate / 3) ** d * (1 - total_err_rate) ** (l - d)
    f_n = int(p_est * (n_hat + 1))
    for f in range(f_n, highest_freq):
        logK = get_log_K(f, f_c, p_no_err, d, l, total_err_rate, logdenom)
        if logK < 0:
            return f

def find_q_p(eps, l):

    for p in range(eps + 1, l):
        q = round(l/p)
        res = l % q
        if res == 0:
            if sum(binom.pmf(x, l, 3/4) for x in range(eps + 1, q*eps + 1)) < 0.5:
                    return q, p
        else:
            if sum(binom.pmf(x, l, 3/4) for x in range(eps + 1, l - (q*((p - eps)-1) + res) + 1)) < 0.5:
                    return q, p


def create_k_mer_dict(seq_list, q, p, l, eps):

    k_mer_dict = {}
    non_trivial_keys = set()
    for seq in seq_list:
        for k_mer_comb in combinations(get_k_mers(seq, q, l), p - eps):
            if k_mer_comb in k_mer_dict:
                k_mer_dict[k_mer_comb].add(seq)
                non_trivial_keys.add(k_mer_comb)
            else:
                k_mer_dict[k_mer_comb] = {seq}


    k_mer_dict = {k: k_mer_dict[k] for k in non_trivial_keys}

    return k_mer_dict

def get_k_mers(seq, q, l):

    return [(int(j / q) + 1, seq[j: j + q]) for j in range(0, l, q)]

def get_candidates(seq, k_mer_dict, q, l, p, eps):

    candidates = set()
    for k_mer_comb in combinations(get_k_mers(seq, q, l), p - eps):
        if k_mer_comb in k_mer_dict:
            candidates.update(k_mer_dict[k_mer_comb])

    return candidates

def get_log_K(f_n, f_c, p_no_err, d, l, total_err_rate, logdenom):

    n_hat = max(int(f_c/p_no_err), f_c + f_n)
    p_est = (total_err_rate/3)**d * (1 - total_err_rate)**(l - d)

    return binom.logpmf(f_n, n_hat, p_est) + math.log(p_est) + logdenom


def locate_mins(a):

    smallest = min(a)
    return smallest, [index for index, element in enumerate(a)
                      if smallest == element]


def cluster_reads(seq_list, seq_to_freq_dict, k_mer_dict, q, p,
                  eps, tau, f, l, p_no_err, total_err_rate, logdenom, bft):

    pb_to_freq_dict = {}
    seq_to_clust_dict = {}
    for i, S_c in enumerate(seq_list):
        f_c = seq_to_freq_dict[S_c]
        if f_c < f:
            candidates = get_candidates(S_c, k_mer_dict, q, l, p, eps)
            if candidates:
                pb_neighbors = [cand for cand in candidates if cand in pb_to_freq_dict.keys()]
                if pb_neighbors:
                    min_dist, indices = locate_mins([trunc_ham_dist(S_c, pb_neighbor, eps, l)
                                                     for pb_neighbor in pb_neighbors])
                    if len(indices) == 1:
                        S_b = pb_neighbors[indices[0]]
                    else:
                        S_b = max([pb_neighbors[j] for j in indices], key=lambda x: seq_to_freq_dict[x])
                    if min_dist != l:
                        if (f_c == 1 and min_dist <= tau) or min_dist == 1:
                            seq_to_clust_dict[S_c] = seq_to_clust_dict[S_b]
                            pb_to_freq_dict[S_b] += f_c
                            continue

                        f_b = seq_to_freq_dict[S_b]
                        logK = get_log_K(f_c, f_b, p_no_err, min_dist, l, total_err_rate, logdenom)
                        if logK > bft:
                            seq_to_clust_dict[S_c] = seq_to_clust_dict[S_b]
                            pb_to_freq_dict[S_b] += f_c
                            continue

        seq_to_clust_dict[S_c] = i
        pb_to_freq_dict[S_c] = f_c

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
            corrected_seq = seq[:i] + seq[i + 1:]
            if corrected_seq in pb_to_freq_dict:
                pb_to_freq_dict[corrected_seq] += freq
                seq_to_clust_dict[seq] = seq_to_clust_dict[corrected_seq]
                break

    return seq_to_clust_dict, pb_to_freq_dict

if __name__ == '__main__':

    my_parser = argparse.ArgumentParser(prog='Shepherd Single',
                                        description='Cluster barcode reads at a single time point')
    my_parser.add_argument('-f', action='store', type=str, required=True, help='Input file name')
    my_parser.add_argument('-l', action='store', type=int, required=True, help='Barcode length')
    my_parser.add_argument('-e', action='store', type=float, help='Substitution error rate estimate')
    my_parser.add_argument('-eps', action='store', type=int, help='Hamming distance threshold')
    my_parser.add_argument('-k', action='store', type=int, help='Substring length')
    my_parser.add_argument('-tau', action='store', type=int, help='Distance threshold for frequency 1 sequences')
    my_parser.add_argument('-ft', action='store', type=int, help='Frequency threshold')
    my_parser.add_argument('-bft', action='store', type=float, help='Bayes factor threshold')
    my_parser.add_argument('-Nh', action='store', type=int, help='Number of sequences used for rho estimation')
    args = my_parser.parse_args()

    filename = args.f
    file_prefix = filename[:-4]
    l = args.l

    deletions_dict = {}
    insertions_dict = {}
    seq_freq_dict = {}
    with open(filename, 'r') as a_file:
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

    seq_list = [seq for seq, freq in sorted(seq_freq_dict.items(), key=lambda x: x[1], reverse=True)]

    if args.Nh == None:
        Nh = 500
    else:
        Nh = args.Nh

    if args.e == None:
        total_err_rate = estimate_rho(seq_list, seq_freq_dict, l, Nh)
        if total_err_rate == 0:
            raise ValueError('Error rate could not be estimated from the data. Please provide an error rate estimate.')
    else:
        total_err_rate = args.e

    print('Substitution Error Rate Estimate: ' + str(total_err_rate))

    highest_freq = seq_freq_dict[seq_list[0]]
    logdenom = l * math.log(4) + math.log(highest_freq)
    p_no_err = binom.pmf(0, l, total_err_rate)

    if args.eps == None:
        eps = find_eps(l, highest_freq, p_no_err, total_err_rate, logdenom)
    else:
        eps = args.eps

    if args.tau == None:
        tau = find_tau(l, p_no_err, total_err_rate, logdenom)
    else:
        tau = args.tau

    if args.ft == None:
        f = find_f(l, highest_freq, p_no_err, total_err_rate, logdenom)
    else:
        f = args.ft

    if args.k == None:
        q, p = find_q_p(eps, l)
    else:
        q = args.k
        p = l // q + int(l % q > 0)

    if args.bft == None:
        bft = -4
    else:
        bft = args.bft

    print('Shepherd Single Parameters:')
    print('\t')
    print('Sequence length: ' + str(l))
    print('Substitution Error Rate Estimate: ' + str(total_err_rate))
    print('epsilon: ' + str(eps))
    print('tau: ' + str(tau))
    print('f: ' + str(f))
    print('Bayes Factor Threshold: ' + str(bft))
    print('Substring Length: ' + str(q))
    print('Number of Partitions: ' + str(p))
    print('\t')

    start = time.time()
    k_mer_dict = create_k_mer_dict(seq_freq_dict, q, p, l, eps)
    end = time.time()
    print('k-mer Index creation time: ' + str(end - start))

    start = time.time()
    seq_to_clust_dict, pb_to_freq_dict = cluster_reads(seq_list, seq_freq_dict, k_mer_dict, q, p, eps, tau, f, l,
                                                        p_no_err, total_err_rate, logdenom, bft)
    seq_to_clust_dict, pb_to_freq_dict = correct_insertions(insertions_dict, pb_to_freq_dict, seq_to_clust_dict, l)
    seq_to_clust_dict, pb_to_freq_dict = correct_deletions(deletions_dict, pb_to_freq_dict, seq_to_clust_dict, l)
    end = time.time()
    print('Clustering time: ' + str(end - start))



    with open(file_prefix + '_index', 'wb') as fh:
        pickle.dump(k_mer_dict, fh)

    with open(file_prefix + '_params', 'wb') as fh:
        pickle.dump([q, l, p, eps, p_no_err, total_err_rate, bft, logdenom, f], fh)

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
    print('Total time: ' + str((end_tot - start_tot)/60))