"""Shepherd: accurate clustering for correcting DNA barcode errors.

Tavakolian et al., Bioinformatics 38(15), 3710-3716 (2022),
https://doi.org/10.1093/bioinformatics/btac395
"""
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version('shepherd-barcodes')
except PackageNotFoundError:  # running from a source checkout that is not installed
    __version__ = 'unknown'
