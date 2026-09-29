"""Command line interface.

    shepherd cluster -f reads_t0.txt -l 20             # single time point
    shepherd track -f0 reads_t0.txt -fn reads_t1.txt   # later time points

The options are the same as those of the original shepherd_t0.py and
shepherd_multi.py scripts.
"""
import argparse

from shepherd import __version__, multi, single


def build_parser():
    parser = argparse.ArgumentParser(
        prog='shepherd',
        description='Shepherd: accurate clustering for correcting DNA barcode errors.')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    subparsers = parser.add_subparsers(dest='command', required=True, metavar='command')

    cluster = subparsers.add_parser(
        'cluster', help='cluster barcode reads at a single time point',
        description='Cluster barcode reads at a single time point')
    cluster.add_argument('-f', action='store', type=str, required=True, help='Input file name')
    cluster.add_argument('-l', action='store', type=int, required=True, help='Barcode length')
    cluster.add_argument('-e', action='store', type=float, help='Substitution error rate estimate')
    cluster.add_argument('-eps', action='store', type=int, help='Hamming distance threshold')
    cluster.add_argument('-k', action='store', type=int, help='Substring length')
    cluster.add_argument('-tau', action='store', type=int, help='Distance threshold for frequency 1 sequences')
    cluster.add_argument('-ft', action='store', type=int, help='Frequency threshold')
    cluster.add_argument('-bft', action='store', type=float, help='Bayes factor threshold')
    cluster.add_argument('-Nh', action='store', type=int, help='Number of sequences used for rho estimation')
    cluster.set_defaults(run=single.run)

    track = subparsers.add_parser(
        'track', help='track barcodes across later time points using the clustering of the first',
        description='Cluster barcode reads at multiple time points')
    track.add_argument('-f0', action='store', type=str, required=True, help='Data file from first time point')
    track.add_argument('-fn', action='store', nargs='+', required=True,
                       help='Ordered list of data files from later time points')
    track.add_argument('-o', action='store', type=str, help='Output file name prefix')
    track.set_defaults(run=multi.run)

    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    args.run(args)
