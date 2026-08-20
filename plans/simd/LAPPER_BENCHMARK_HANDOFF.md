# rust-lapper 2 benchmark handoff

> This is a reproducibility record for the exact pre-release benchmark commit
> named below, not current installation or API guidance. The release-facing
> summary lives in [`README.md`](../../README.md).

## Objective

Reproduce the rust-lapper 2 release performance check without relying on the
existing `lapper_bakeoff` checkout or its dirty working tree. The required
implementation revisions, workload semantics, data source, expected counts,
and reference results are below.

The release question is: does the portable block/SIMD implementation avoid a
regression against rust-lapper 1.3.0 on the three retained article workloads,
especially the pathological `7-3` case?

## Exact revisions and build

- Candidate query implementation: `sstadick/rust-lapper`, branch
  `worked/portable-simd-index`, commit
  `0eeeafdd0773e10fc1ef7f007e69f471917562d7`.
- Baseline: rust-lapper 1.3.0, commit
  `c545d40` (also tag `v1.3.0`).
- Compiler used for the retained x86 record: Rust 1.95.0.
- Release flags: `RUSTFLAGS="-C target-cpu=native"`.
- Release profile: thin LTO and one code-generation unit.

The later local 2.0.0 release-preparation edits change metadata, documentation,
tests, and development dependencies, not the measured query implementation, so
`0eeeafd` is the correct reproducible performance revision.

Prefer a physical x86-64 host with AVX2. Record `lscpu`, OS/kernel, rustc
version, and the full build flags. Before timing the candidate, verify dispatch:

```sh
cargo +1.95.0 test --all-features --locked \
  simd::dispatch_tests::selected_backend_matches_the_host -- --exact --nocapture
cargo +1.95.0 test --all-features --locked \
  simd::avx2::tests::primitive_masks_match_scalar_when_avx2_is_available \
  -- --exact --nocapture
```

The first test must report/select AVX2 on an AVX2 benchmark host, and the direct
primitive-mask test must pass.

## Data

The source archive is approximately 401 MB:

```text
https://drive.usercontent.google.com/download?id=1lctmude31mSAh9fWjI60K1bDrbeDPGfm&export=download&confirm=t
```

It contains Parquet directories under `databio/`. Convert these three columns
to headerless, tab-separated BED using Polars 1.32.3:

```text
contig, pos_start, pos_end
```

Required directories and integrity counts:

| Article ID | Source directory | BED name | Intervals |
|---|---|---|---:|
| 1 | `fBrain-DS14718` | `fBrain.bed` | 198,621 |
| 2 | `exons` | `exons.bed` | 438,694 |
| 3 | `chainOrnAna1` | `chainOrnAna1.bed` | 1,956,864 |
| 7 | `ex-anno` | `ex-anno.bed` | 1,194,285 |
| 8 | `ex-rna` | `ex-rna.bed` | 9,944,559 |

All coordinates are zero-based, half-open `[start, stop)` ranges.

## Harness contract

A small temporary Rust binary is sufficient; the five-library bakeoff is not
required for the v1-versus-v2 release gate.

1. Read both BED files before starting any timer.
2. Group ranges by chromosome in a `BTreeMap`.
3. For each reported repeat, build a fresh
   `BTreeMap<String, Lapper<u32, ()>>`, one `Lapper` per chromosome, and time
   construction separately.
4. Time queries by iterating every query range and evaluating
   `index.find(start, stop).count()`. Sum and black-box the result.
5. Do **not** substitute `Lapper::count`; the workload measures result
   enumeration through `find`.
6. Emit raw TSV fields:
   `revision, repeat, database_intervals, query_intervals, overlaps, build_seconds, query_seconds`.
7. Run at least three reported repeats per revision and case. Retain every
   sample. Report build and query medians independently, then add them for the
   total.

Use each case in this physical database/query orientation:

| Case | Indexed database | Query input | Expected overlaps |
|---|---|---|---:|
| `1-2` | `exons.bed` | `fBrain.bed` | 54,246 |
| `7-3` | `chainOrnAna1.bed` | `ex-anno.bed` | 4,408,383 |
| `8-7` | `ex-anno.bed` | `ex-rna.bed` | 307,184,634 |

The orientation is intentional. Case `7-3` reproduces rust-lapper 1.3.0's
global-maximum-interval-length worst case. Abort if either revision returns a
different overlap count.

For clean dependency switching, create two temporary worktrees at the exact
commits and build the same harness once against each path. Do not benchmark a
debug build. Alternate revision execution order where practical, and keep the
machine otherwise idle.

## Reference result

The retained native AVX2 run used an AMD Ryzen 9 3950X, Ubuntu 24.04.3, Linux
6.8.0-100-generic, rustc 1.95.0, thin LTO, one codegen unit, and
`-C target-cpu=native`. Independent three-sample medians were:

| Case | v1.3.0 total | Candidate total | Candidate change |
|---|---:|---:|---:|
| `1-2` | 13.715 ms | 8.600 ms | -37.30% |
| `7-3` | 8,756.427 ms | 97.084 ms | -98.89% |
| `8-7` | 1,192.626 ms | 779.564 ms | -34.64% |

The `7-3` total improved by 90.2 times. Absolute timings will vary by host; the
release gate is identical counts and no unexplained regression versus v1.3.0.

## Deliverables

Retain:

- raw TSV samples and stderr/load logs;
- CPU, OS/kernel, compiler, commit IDs, and exact flags;
- interval and overlap-count checks;
- the untrimmed median calculation; and
- a short Markdown summary stating whether any case regressed.

The repository's existing `cargo bench --bench lapper_benchmark` is useful as a
synthetic smoke benchmark, but it does not replace these three release
workloads.
