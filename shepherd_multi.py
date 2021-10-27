from itertools import combinations
import pickle
import json
import time
import argparse
import csv

def trunc_ham_dist(seq_1, seq_2, d, n):

    h = 0
    for n_1, n_2 in zip(seq_1, seq_2):
        if n_1 != n_2:
            h += 1
        if h > d:
            return n

    return h

def get_k_mers(seq, q, l):

    return [(int(j / q) + 1, seq[j: j + q]) for j in range(0, l, q)]

def get_candidates(seq, k_mer_dict, q, l, p, eps):

    candidates = set()
    for k_mer_comb in combinations(get_k_mers(seq, q, l), p - eps):
        if k_mer_comb in k_mer_dict:
            candidates.update(k_mer_dict[k_mer_comb])

    return candidates

def locate_mins(a):

    smallest = min(a)
    return smallest, [index for index, element in enumerate(a)
                      if smallest == element]


def get_closest_pb(seq, pb_to_freq_dict_t0, pb_to_freq_dict_t1, k_mer_dict, q, l, p, eps):

    candidates = get_candidates(seq, k_mer_dict, q, l, p, eps)
    if candidates:
        pb_neighbors = [cand for cand in candidates if cand in pb_to_freq_dict_t1]
        if pb_neighbors:
            min_dist, indices = locate_mins([trunc_ham_dist(seq, pb_neighbor, eps, l)
                                             for pb_neighbor in pb_neighbors])
            if len(indices) == 1:
                closest_pb = pb_neighbors[indices[0]]
            else:
                closest_pb = max([pb_neighbors[j] for j in indices],
                                 key=lambda x: pb_to_freq_dict_t0[x])

            return closest_pb


def classify_reads(seq_freq_list_t1, k_mer_dict, pb_to_freq_dict_t0, q, l, p, eps):

    pb_to_freq_dict_t1 = {pb: 0 for pb, freq in pb_to_freq_dict_t0.items() if freq != 0}
    for seq, freq in seq_freq_list_t1:
        if seq in pb_to_freq_dict_t1:
            pb_to_freq_dict_t1[seq] += freq
        else:
            closest_pb = get_closest_pb(seq, pb_to_freq_dict_t0, pb_to_freq_dict_t1, k_mer_dict, q, l, p, eps)
            if closest_pb:
                pb_to_freq_dict_t1[closest_pb] += freq

    return pb_to_freq_dict_t1

if __name__ == '__main__':

    my_parser = argparse.ArgumentParser(prog='Shepherd Multi',
                                        description='Cluster barcode reads at multiple time points')
    my_parser.add_argument('-f0', action='store', type=str, required=True, help='Data file from first time point')
    my_parser.add_argument('-fn', action='store', nargs='+', help='Ordered list of data files from later time points')
    my_parser.add_argument('-o', action='store', type=str, help='Output file name prefix')
    args = my_parser.parse_args()

    f0_prefix = args.f0[:-4]
    multi_filenames = args.fn

    if args.o == None:
        o_fn_prefix = 'multi_freqs'
    else:
        o_fn_prefix = args.o

    with open(f0_prefix + '_index', 'rb') as handle:
        k_mer_dict = pickle.load(handle)

    with open(f0_prefix + '_pb_freq.csv', 'r') as handle:
        pb_to_freq_dict = {}
        i = 0
        for line in handle:
            if i > 0:
                key, value = line.split(',')
                pb_to_freq_dict[key] = int(value)
            i += 1

    with open(f0_prefix + '_params', 'rb') as f:
        q, l, p, eps = pickle.load(f)

    print('Starting Classification')

    start = time.time()
    pb_to_freq_dict_list = [pb_to_freq_dict]
    for filename in multi_filenames:
        seq_freq_list = []
        with open(filename, 'r') as a_file:
            for line in a_file:
                key, value = line.split()
                seq_freq_list.append((key, int(value)))

        pb_to_freq_dict = classify_reads(seq_freq_list, k_mer_dict, pb_to_freq_dict, q, l, p, eps)
        pb_to_freq_dict_list.append(pb_to_freq_dict)
        print(len(pb_to_freq_dict))

    end = time.time()
    print('Classification: ' + str(end - start))

    with open(o_fn_prefix + '.csv', 'w', newline='') as result:
        writer = csv.writer(result, delimiter=",")
        col_names = ['barcode']
        time_points = ['time_point_' + str(i) for i in range(1, len(pb_to_freq_dict_list) + 1)]
        col_names.extend(time_points)
        writer.writerow(col_names)
        for seq in pb_to_freq_dict_list[0].keys():
            line = [seq]
            for pb_to_freq_dict in pb_to_freq_dict_list:
                if seq in pb_to_freq_dict:
                    line.append(pb_to_freq_dict[seq])
                else:
                    line.append(0)

            writer.writerow(line)