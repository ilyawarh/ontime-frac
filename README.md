# ontime-frac

Sample Oxford Nanopore reads by **sequencing-time fraction**: extract the earliest (or
latest) N% of reads by their `st:Z` start time, where N is a fraction of the total read
*count* — not a time span, and not a random subset.

Built on top of [ontime](https://github.com/mbhall88/ontime) for the heavy lifting.

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

# fractions in percents and suffix auto completion, use 16 threads
ontime-frac.sh -i reads.fq -o 'subsample' -f 10,25,50 -j 16

# direct mode (without decompression) with nested subsets below 0.2
ontime-frac.sh -i reads.fq.gz -o 'sampled_{frac}.fq.gz' -f 0.1,0.25,0.5 --mode direct --nested-below 0.2
```

Outputs are **cumulative nested subsets**: the 10% output is contained in the 25% output,
etc. Counts are exact (modulo timestamp ties at the cutoff): the cutoff is the timestamp
of the ceil(total × frac)-th read and ontime keeps everything up to it.

## Options

```
-i PATH              input FASTQ (.fq/.fastq, optionally .gz)
-o TEMPLATE          output template; '{frac}' is replaced by each fraction
                         (without '{frac}', '_<frac>' is inserted before the extension)
-f LIST              comma-separated fractions (0-1, or 0-100 if given as percent)
--from-end           sample the LATEST fraction(s) instead of earliest
--mode M             auto (default) | decompress | direct
                         decompress: gunzip once to plain FASTQ, all passes read plain
                         (ontime is ~15x faster on plain input); direct: read .gz as-is
--nested-below B     fractions < B use the nested cascade, >= B run in
                         parallel from the full input (default 0.5; 0 = all parallel,
                         1 = all nested)
--workdir DIR        directory for intermediates (default: input's directory)
--no-keep            delete the sorted-timestamp cache after the run
--force              overwrite existing outputs
-j, --threads N      threads for pigz/seqkit/sort (default: nproc)
-h                   this help
```

## How it works

- Stage 1: seqkit extracts st:Z timestamps (fused with decompression if .gz)
- Stage 2: GNU sort (lexicographic == chronological for same-timezone RFC3339)
- Stage 3: cutoff = timestamp of the ceil(total*frac)-th read; ontime filters.
            Fractions >= --nested-below run in parallel from the full input;
            smaller fractions run as a nested cascade from the next-larger
            output (each pass reads less data).

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
