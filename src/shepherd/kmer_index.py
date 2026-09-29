"""The k-mer Index (Section 2.1 of the paper).

A sequence of length l is split into p non-overlapping k-mers (the last one
may be shorter). By the pigeonhole principle, two sequences within Hamming
distance eps share at least p - eps of them. The index maps every combination
of p - eps position-tagged k-mers to the set of indexed sequences containing
it, so all sequences within distance eps of a query can be found by looking
up the query's own combinations.
"""

from itertools import combinations


def get_k_mers(seq, q, l):

    return [(int(j / q) + 1, seq[j : j + q]) for j in range(0, l, q)]


def add_seq_to_k_mer_dict(seq, k_mer_dict, q, l, p, eps):
    for k_mer_comb in combinations(get_k_mers(seq, q, l), p - eps):
        if k_mer_comb in k_mer_dict:
            k_mer_dict[k_mer_comb].add(seq)
        else:
            k_mer_dict[k_mer_comb] = {seq}

    return k_mer_dict


def build_k_mer_dict(seq_list, q, l, p, eps):
    k_mer_dict = {}
    for seq in seq_list:
        add_seq_to_k_mer_dict(seq, k_mer_dict, q, l, p, eps)

    return k_mer_dict


def get_candidates(seq, k_mer_dict, q, l, p, eps):

    candidates = set()
    for k_mer_comb in combinations(get_k_mers(seq, q, l), p - eps):
        if k_mer_comb in k_mer_dict:
            candidates.update(k_mer_dict[k_mer_comb])

    return candidates
