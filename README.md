# ontime-frac

Sample Oxford Nanopore reads by **sequencing-time fraction**: extract the earliest (or
latest) N% of reads by their `st:Z` start time, where N is a fraction of the total read
*count* — not a time span, and not a random subset.

Built on top of [ontime](https://github.com/mbhall88/ontime) for the heavy lifting.

##Why

Dorado basecalls on GPU in batches and writes reads in **completion order, not time
order**. So the first 10% of lines in a basecalled FASTQ is a random temporal sample, not
the first 10% of sequencing. To get "the reads sequenced in the first 10% of the run by
count", you must rank by `st:Z` — that is what this tool does, without ever sorting the
FASTQ itself.

## Install

```bash
conda install -c bioconda seqkit ontime pigz

git clone git@github.com:ilyawarh/ontime-frac.git
cd ontime-frac

bash ontime-frac [options] ...
```

## Quick start

```bash
# earliest 10%, 25% and 50% by sequencing time, in one run
ontime-frac.sh -i reads.fq.gz -o 'sampled_{frac}.fq.gz' -f 0.1,0.25,0.5

# percents work too; latest fraction instead of earliest
ontime-frac.sh -i reads.fq.gz -o 'late_{frac}.fq.gz' -f 10 --from-end

# serial experiments: the second run reuses the cached timestamps and
# decompressed copy, going straight to filtering
ontime-frac.sh -i reads.fq.gz -o 's_{frac}.fq.gz' -f 0.02,0.05
```

Outputs are **cumulative nested subsets**: the 10% output is contained in the 25% output,
etc. Counts are exact (modulo timestamp ties at the cutoff): the cutoff is the timestamp
of the ceil(total × frac)-th read and ontime keeps everything up to it.

## Options

| Option | Default | Meaning |
|---|---|---|
| `-i` | — | input FASTQ (`.fq`/`.fastq`, optionally `.gz`) |
| `-o` | — | output template; `{frac}` is replaced by each fraction |
| `-f` | — | comma-separated fractions (0-1, or 0-100 as percent) |
| `--from-end` | off | sample the latest fraction(s) instead of the earliest |
| `--mode` | `auto` | `auto` / `decompress` / `direct` (see below) |
| `--nested-below` | `0.5` | fractions below this use the nested cascade; at/above run in parallel. `0` = all parallel, `1` = all nested |
| `--workdir` | input's dir | where intermediates (decompressed copy, timestamp cache) live |
| `--no-keep` | off | delete the timestamp cache after the run |
| `--force` | off | overwrite existing outputs |
| `-j`, `--threads` | `nproc` | threads for pigz / seqkit / sort |

## How it works

```
stage 1   pigz -dc in.gz | tee plain.fq | seqkit seq -n - | sed  →  timestamps
stage 2   GNU sort  →  timestamps.sorted.txt   (cached and reused across runs)
stage 3   per fraction: cutoff = sorted[ceil(total × frac)]
          ontime --to <cutoff>  does the actual filtering
```

- **Decompress mode** (auto-chosen for `.gz` when disk allows): ontime reads plain text at
  ~1030 MB/s vs ~67 MB/s through its own single-threaded gunzipper, so the file is
  decompressed once with pigz and every pass reads the plain copy. Stage 1 is fused with
  decompression via `tee`, so the `.gz` is read exactly once.
- **Nested vs parallel**: fractions ≥ `--nested-below` are filtered in parallel from the
  full input (nesting a 0.5 from a 0.75 still reads 75% of the file — little to gain).
  Fractions below it form a cascade: each pass filters the *previous, larger output*
  (0.25 from the 0.5 output, 0.1 from the 0.25 output), so small fractions never touch the
  full file. This is the regime that matters on 100+ GB disk-bound files: total read
  volume is `1 + f_k + … + f₂` instead of `k × 1`.
- Output compression is `ontime | pigz`, ~6.5× faster than ontime's internal gzip.

## Citing

ontime-frac wraps `ontime` — please cite it:

> Hall, M. (2023). mbhall88/ontime: 0.1.3. Zenodo. https://doi.org/10.5281/zenodo.7533053
