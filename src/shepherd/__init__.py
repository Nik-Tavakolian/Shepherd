"""Shepherd: accurate clustering for correcting DNA barcode errors.

Tavakolian et al., Bioinformatics 38(15), 3710-3716 (2022),
https://doi.org/10.1093/bioinformatics/btac395

Use from Python::

    from shepherd import Tracker, cluster, estimate_parameters, read_counts

    reads = read_counts('reads_t0.txt', barcode_length=20)
    params = estimate_parameters(reads.barcodes, barcode_length=20)
    clustering = cluster(reads, params)

    tracker = Tracker(clustering.barcode_counts, params)
    tracker.add_time_point(read_counts('reads_t1.txt', barcode_length=20))
"""

from importlib.metadata import PackageNotFoundError, version

from shepherd.clustering import Clustering, cluster
from shepherd.errors import ShepherdError
from shepherd.io import ReadCounts, read_counts
from shepherd.parameters import Parameters, estimate_parameters
from shepherd.tracking import Tracker

try:
    __version__ = version('shepherd-barcodes')
except PackageNotFoundError:  # running from a source checkout that is not installed
    __version__ = 'unknown'

__all__ = [
    'Clustering',
    'Parameters',
    'ReadCounts',
    'ShepherdError',
    'Tracker',
    'cluster',
    'estimate_parameters',
    'read_counts',
]
