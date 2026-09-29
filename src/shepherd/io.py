"""Reading Shepherd's input files and writing its output files."""

import csv
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from shepherd.model import ShepherdError


@dataclass
class ReadCounts:
    """The unique sequences of one time point and their read counts, split by length."""

    barcodes: dict[str, int] = field(default_factory=dict)
    """Sequences of the barcode length l."""
    insertions: dict[str, int] = field(default_factory=dict)
    """Sequences of length l + 1, possibly with a single insertion error."""
    deletions: dict[str, int] = field(default_factory=dict)
    """Sequences of length l - 1, possibly with a single deletion error."""


def read_counts(path: str | Path, barcode_length: int) -> ReadCounts:
    """Read a file with a sequence and its read count, separated by whitespace, on each line.

    Sequences whose length is not l or l +/- 1 are ignored.
    """
    by_length: dict[int, dict[str, int]] = {
        barcode_length: {},
        barcode_length + 1: {},
        barcode_length - 1: {},
    }
    with open(path) as fh:
        for line_number, line in enumerate(fh, start=1):
            fields = line.split()
            if not fields:
                continue
            try:
                seq, count = fields
                counts = by_length.get(len(seq))
                if counts is not None:
                    counts[seq] = int(count)
            except ValueError:
                raise ShepherdError(
                    f'{path}, line {line_number}: expected a sequence and a read count, '
                    f'got {line.strip()!r}'
                ) from None
    return ReadCounts(
        barcodes=by_length[barcode_length],
        insertions=by_length[barcode_length + 1],
        deletions=by_length[barcode_length - 1],
    )


def output_path(input_path: str | Path, suffix: str) -> Path:
    """The path of an output file: the input path without its extension, plus suffix."""
    input_path = Path(input_path)
    return input_path.with_name(input_path.stem + suffix)


def write_labels(path: str | Path, labels: Mapping[str, int]) -> None:
    _write_csv(path, ['sequence', 'cluster'], labels.items())


def write_barcode_counts(path: str | Path, counts: Mapping[str, int]) -> None:
    _write_csv(path, ['barcode', 'frequency'], counts.items())


def read_barcode_counts(path: str | Path) -> dict[str, int]:
    """Read a file written by :func:`write_barcode_counts`."""
    with open(path, newline='') as fh:
        rows = csv.reader(fh)
        next(rows)
        return {barcode: int(count) for barcode, count in rows}


def write_count_table(path: str | Path, counts_per_time_point: Sequence[Mapping[str, int]]) -> None:
    """Write one row per barcode with its read count at each time point (0 if absent).

    Barcodes appear in the order in which they are first seen.
    """
    barcodes = dict.fromkeys(barcode for counts in counts_per_time_point for barcode in counts)
    header = ['barcode'] + [f'time_point_{i}' for i in range(1, len(counts_per_time_point) + 1)]
    rows = (
        [barcode] + [counts.get(barcode, 0) for counts in counts_per_time_point]
        for barcode in barcodes
    )
    _write_csv(path, header, rows)


def _write_csv(path: str | Path, header: list[str], rows: Iterable[Iterable[object]]) -> None:
    with open(path, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)
