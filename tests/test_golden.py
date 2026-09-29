"""Check that Shepherd's outputs match the recorded golden outputs.

The golden outputs were recorded with tests/golden/make_golden.py. Rows are
compared as mappings, so a change in row order is not a failure.
"""

import csv
import gzip
import io
from pathlib import Path

import pytest

from golden.make_golden import DATA_DIR, OUTPUT_FILES, run_all_scenarios


def read_csv_rows(text):
    header, *rows = csv.reader(io.StringIO(text))
    return header, {row[0]: row[1:] for row in rows}


@pytest.fixture(scope='module')
def outputs(tmp_path_factory):
    workdir = tmp_path_factory.mktemp('golden')
    run_all_scenarios(workdir)
    return workdir


@pytest.mark.parametrize('filename', OUTPUT_FILES)
def test_output_matches_golden(outputs, filename):
    expected = gzip.decompress((DATA_DIR / (filename + '.gz')).read_bytes()).decode()
    actual = (Path(outputs) / filename).read_text()
    assert read_csv_rows(actual) == read_csv_rows(expected)
