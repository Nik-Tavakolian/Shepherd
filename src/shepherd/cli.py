"""The ``shepherd`` command.

    shepherd cluster -f reads_t0.txt -l 20                         # first time point
    shepherd track -f0 reads_t0.txt -fn reads_t1.txt reads_t2.txt  # later time points

The short options are those of the original shepherd_t0.py and
shepherd_multi.py scripts; each also has a long name.
"""

import argparse
import logging
from collections.abc import Sequence

from shepherd import __version__
from shepherd.clustering import cluster
from shepherd.io import (
    output_path,
    read_barcode_counts,
    read_counts,
    write_barcode_counts,
    write_count_table,
    write_labels,
)
from shepherd.model import (
    DEFAULT_LOG_BF_THRESHOLD,
    DEFAULT_N_TOP,
    Parameters,
    ShepherdError,
    estimate_parameters,
)
from shepherd.tracking import Tracker

logger = logging.getLogger('shepherd')


def main(argv: Sequence[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(format='%(message)s', level=logging.INFO)
    try:
        args.run(args)
    except (ShepherdError, OSError) as error:
        parser.exit(1, f'shepherd: error: {error}\n')


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog='shepherd', description='Accurate clustering for correcting DNA barcode errors.'
    )
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    commands = parser.add_subparsers(required=True, metavar='command')

    cluster_parser = commands.add_parser(
        'cluster', help='cluster the reads of the first time point'
    )
    cluster_parser.add_argument('-f', '--input', required=True, help='sequences and read counts')
    cluster_parser.add_argument('-l', '--length', type=int, required=True, help='barcode length')
    cluster_parser.add_argument('-e', '--error-rate', type=float, help='substitution error rate')
    cluster_parser.add_argument('-eps', '--epsilon', type=int, help='maximum merge distance')
    cluster_parser.add_argument('-k', '--kmer-length', type=int, help='k-mer length')
    cluster_parser.add_argument(
        '-tau', '--tau', type=int, help='merge distance for one-read sequences'
    )
    cluster_parser.add_argument('-ft', '--count-threshold', type=int, help='count of sure barcodes')
    cluster_parser.add_argument(
        '-bft',
        '--log-bf-threshold',
        type=float,
        default=DEFAULT_LOG_BF_THRESHOLD,
        help='log Bayes factor threshold (default: %(default)s)',
    )
    cluster_parser.add_argument(
        '-Nh',
        '--n-top',
        type=int,
        default=DEFAULT_N_TOP,
        help='sequences used to estimate the error rate (default: %(default)s)',
    )
    cluster_parser.set_defaults(run=run_cluster)

    track_parser = commands.add_parser(
        'track', help='follow the barcodes of the first time point through later ones'
    )
    track_parser.add_argument('-f0', '--first', required=True, help='input of the first time point')
    track_parser.add_argument('-fn', '--later', nargs='+', required=True, help='later inputs')
    track_parser.add_argument(
        '-o', '--output', default='multi_freqs', help='count table name (default: %(default)s)'
    )
    track_parser.set_defaults(run=run_track)

    return parser


def run_cluster(args: argparse.Namespace) -> None:
    reads = read_counts(args.input, args.length)
    params = estimate_parameters(
        reads.barcodes,
        args.length,
        error_rate=args.error_rate,
        epsilon=args.epsilon,
        kmer_length=args.kmer_length,
        tau=args.tau,
        count_threshold=args.count_threshold,
        log_bf_threshold=args.log_bf_threshold,
        n_top=args.n_top,
    )
    logger.info('Parameters: %s', params)

    clustering = cluster(reads, params)
    logger.info('Found %d putative barcodes', len(clustering.barcode_counts))

    write_labels(output_path(args.input, '_seq_clust.csv'), clustering.labels)
    write_barcode_counts(output_path(args.input, '_pb_freq.csv'), clustering.barcode_counts)
    params.save(output_path(args.input, '_params.json'))


def run_track(args: argparse.Namespace) -> None:
    counts_path = output_path(args.first, '_pb_freq.csv')
    if not counts_path.exists():
        raise ShepherdError(
            f'{counts_path} not found; run `shepherd cluster` on {args.first} first'
        )
    params = Parameters.load(output_path(args.first, '_params.json'))
    tracker = Tracker(read_barcode_counts(counts_path), params)

    for filename in args.later:
        clustering = tracker.add_time_point(read_counts(filename, params.barcode_length))
        logger.info('%s: %d putative barcodes', filename, len(clustering.barcode_counts))
        write_labels(output_path(filename, '_seq_clust.csv'), clustering.labels)

    write_count_table(args.output + '.csv', tracker.counts_per_time_point)
