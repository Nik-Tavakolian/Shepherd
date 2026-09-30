# Example data

Small data sets that are easy to check by hand. The barcodes are 20 nt
patterns such as `AAAAACCCCCGGGGGTTTTT`, so errors stand out. Each file has
a sequence and its read count per line, like any Shepherd input.

With this few reads the error rate cannot be estimated reliably, so the
commands give it with `-e 0.01`; the other parameters are determined from the
data (epsilon = 3, tau = 3, count threshold 22-24). The expected results below
are also checked by `tests/test_examples.py`.

Run the commands from this folder.

## Single time point

### single_basic.txt: three barcodes and their errors

```bash
shepherd cluster -f single_basic.txt -l 20 -e 0.01
```

| Barcode | Reads | Error sequences | Expected count |
|---|---|---|---|
| `AAAAACCCCCGGGGGTTTTT` | 1000 | three at distance 1 (4 + 3 + 2 reads), one at distance 2 (1 read) | 1010 |
| `CCCCCAAAAATTTTTGGGGG` | 500 | one at distance 1 (3 reads), one at distance 2 (1 read) | 504 |
| `GGGGGTTTTTAAAAACCCCC` | 50 | none | 50 |

`single_basic_pb_freq.csv` lists these three barcodes, and
`single_basic_seq_clust.csv` gives each error sequence the cluster label of
its barcode.

### single_close_barcodes.txt: two true barcodes two substitutions apart

```bash
shepherd cluster -f single_close_barcodes.txt -l 20 -e 0.01
```

`ACGTACGTACGTACGTAGGA` (15 reads) is two substitutions away from
`ACGTACGTACGTACGTACGT` (800 reads). With a 1% error rate, the barcode with 800
reads is expected to produce about 0.01 reads of any particular sequence at
distance 2, so 15 reads are far too many for an error sequence: the Bayesian
test gives ln K = -75.7 and it stays a separate barcode. A sequence with 2 reads
at the same distance, `ACGTACGTACGTACGTTTGT`, is an error sequence (ln K = 12.7).

| Barcode | Expected count |
|---|---|
| `ACGTACGTACGTACGTACGT` | 805 (800 + 2 + 3) |
| `ACGTACGTACGTACGTAGGA` | 15 |

### single_indels.txt: insertion and deletion errors

```bash
shepherd cluster -f single_indels.txt -l 20 -e 0.01
```

`AAAAACCCCCGGGGGGTTTTT` (21 nt, an extra G) and `AAAAACCCCGGGGGTTTTT` (19 nt,
a missing C) are merged into `AAAAACCCCCGGGGGTTTTT`, which gets 1000 + 7 + 4 =
1011 reads. `TGCATGCATGCATGCATGCAT` (21 nt) is not one insertion away from a
barcode and `ACGTACGTACGTACG` (15 nt) has the wrong length, so neither appears
in the output.

| Barcode | Expected count |
|---|---|
| `AAAAACCCCCGGGGGTTTTT` | 1011 |
| `CCCCCAAAAATTTTTGGGGG` | 200 |

## Multiple time points

```bash
shepherd cluster -f multi_t0.txt -l 20 -e 0.01
```

```bash
shepherd track -f0 multi_t0.txt -fn multi_t1.txt multi_t2.txt
```

Expected `multi_freqs.csv` (time point 1 is t0):

| Barcode | t0 | t1 | t2 | What happens |
|---|---|---|---|---|
| `AAAAACCCCCGGGGGTTTTT` (A) | 1006 | 905 | 802 | keeps its error sequences at every time point |
| `CCCCCAAAAATTTTTGGGGG` (B) | 503 | 603 | 704 | grows |
| `GGGGGTTTTTAAAAACCCCC` (C) | 101 | 0 | 0 | dies out after t0 |
| `AAAAACCCCCGGGGGTTTAA` (E) | 0 | 300 | 350 | emerges at t1, two substitutions from A |
| `TTTTTGGGGGCCCCCAAAAA` (N) | 0 | 204 | 405 | emerges at t1, far from all barcodes |

E and N show the two ways barcodes emerge after the first time point
(Supplementary Section 2):

- E lies within epsilon of A, so at t1 it is first assigned to A. With 300
  reads it fails the Bayesian test as an error sequence of A and is separated
  into its own barcode at the same time point.
- N has no barcode within epsilon, so at t1 it and its two error sequences are
  clustered into a candidate barcode. Candidates are kept only if they are seen
  again at the next time point: N is confirmed at t2, and its count at t1 (200 +
  3 + 1 = 204) is then filled in. Because it is still a candidate after t1, N
  and its error sequences are not listed in `multi_t1_seq_clust.csv`.
