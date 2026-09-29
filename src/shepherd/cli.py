"""Command line interface.

    shepherd cluster -f reads_t0.txt -l 20             # single time point
    shepherd track -f0 reads_t0.txt -fn reads_t1.txt   # later time points

The options are the same as those of the original shepherd_t0.py and
shepherd_multi.py scripts.
"""

import argparse
import time

from shepherd import __version__, multi
from shepherd.clustering import cluster
from shepherd.io import output_path, read_counts, write_barcode_counts, write_labels
from shepherd.parameters import DEFAULT_LOG_BF_THRESHOLD, DEFAULT_N_TOP, estimate_parameters


def build_parser():
    parser = argparse.ArgumentParser(
        prog='shepherd',
        description='Shepherd: accurate clustering for correcting DNA barcode errors.',
    )
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    subparsers = parser.add_subparsers(dest='command', required=True, metavar='command')

    cluster_parser = subparsers.add_parser(
        'cluster',
        help='cluster barcode reads at a single time point',
        description='Cluster barcode reads at a single time point',
    )
    cluster_parser.add_argument(
        '-f', action='store', type=str, required=True, help='Input file name'
    )
    cluster_parser.add_argument(
        '-l', action='store', type=int, required=True, help='Barcode length'
    )
    cluster_parser.add_argument(
        '-e', action='store', type=float, help='Substitution error rate estimate'
    )
    cluster_parser.add_argument('-eps', action='store', type=int, help='Hamming distance threshold')
    cluster_parser.add_argument('-k', action='store', type=int, help='Substring length')
    cluster_parser.add_argument(
        '-tau', action='store', type=int, help='Distance threshold for frequency 1 sequences'
    )
    cluster_parser.add_argument('-ft', action='store', type=int, help='Frequency threshold')
    cluster_parser.add_argument('-bft', action='store', type=float, help='Bayes factor threshold')
    cluster_parser.add_argument(
        '-Nh', action='store', type=int, help='Number of sequences used for rho estimation'
    )
    cluster_parser.set_defaults(run=run_cluster)

    track = subparsers.add_parser(
        'track',
        help='track barcodes across later time points using the clustering of the first',
        description='Cluster barcode reads at multiple time points',
    )
    track.add_argument(
        '-f0', action='store', type=str, required=True, help='Data file from first time point'
    )
    track.add_argument(
        '-fn',
        action='store',
        nargs='+',
        required=True,
        help='Ordered list of data files from later time points',
    )
    track.add_argument('-o', action='store', type=str, help='Output file name prefix')
    track.set_defaults(run=multi.run)

    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    args.run(args)


def run_cluster(args: argparse.Namespace) -> None:
    reads = read_counts(args.f, args.l)
    start = time.time()
    params = estimate_parameters(
        reads.barcodes,
        args.l,
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

    clustering = cluster(reads, params)
    print('Clustering time: ' + str(time.time() - start))

    params.save(output_path(args.f, '_params.json'))
    write_labels(output_path(args.f, '_seq_clust.csv'), clustering.labels)
    write_barcode_counts(output_path(args.f, '_pb_freq.csv'), clustering.barcode_counts)
