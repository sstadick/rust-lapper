# Portable SIMD index production plan

The implementation on `worked/portable-simd-index` is feature-complete as an
experimental release candidate. Productionization is limited to the five gates
below; changes should preserve the existing public API, forward borrowed
iteration order, exact overlap semantics, and always-on query algorithm.

## 1. Validate AVX2 on real Intel/AMD hardware

Status: in progress. On 2026-07-28, GitHub's native x86-64 runner reported an
AMD EPYC 7763 with AVX2. Runtime dispatch selected AVX2 and every signed and
unsigned AVX2 mask matched the scalar result. Rosetta separately exercised the
x86-64 scalar selection. Native three-dataset performance measurements and a
physical non-AVX2 x86 host remain outstanding.

- Run the complete test suite on native x86-64 with AVX2 available.
- Verify that the AVX2 backend is selected and executed, rather than only
  inspecting forced-target assembly.
- Run the three retained datasets against rust-lapper 1.3.0, the worked Lapper,
  and the pinned Rust competitors.
- Test an x86-64 build or host without AVX2 and confirm the scalar fallback.
- Record CPU, compiler, build flags, raw samples, medians, and overlap counts.

Pass condition: native AVX2 and non-AVX2 paths are correct, and the native AVX2
performance record contains no unexplained regression.

## 2. Decide and document the effective MSRV

Status: complete. The declared MSRV is Rust 1.59.0. The locked all-feature
library builds on 1.59.0, while 1.58.1 fails because its AArch64 `std::arch`
intrinsics are still unstable. Rust-lapper 1.3.0 had no declared MSRV and its
locked all-feature library still builds on Rust 1.56.1.

- Determine the oldest compiler supported by rust-lapper 1.3.0.
- Identify the oldest compiler accepted by the worked implementation and its
  dependencies.
- Either retain the existing effective MSRV or choose a deliberate new one.
- Add the decision to package metadata, CI, and release documentation.

Pass condition: the declared MSRV builds and tests the supported feature set,
and any increase is intentional and documented.

## 3. Resolve the `I: 'static` API-bound addition

Status: complete and accepted. Stable Rust has no specialization mechanism that
can select primitive SIMD kernels while retaining a blanket scalar path for
every custom `PrimInt`. Safe pointer reinterpretation therefore uses exact
`TypeId` checks, which require `I: 'static`. This excludes only custom coordinate
types carrying non-static borrows; it does not require a `Lapper` value to live
for the program lifetime. Primitive types use SIMD, while `u128`, `i128`, and
other unmatched owned `PrimInt` types use the tested scalar fallback.

- Measure the public API difference from rust-lapper 1.3.0.
- Determine whether safe primitive-type SIMD dispatch can avoid `TypeId` and
  the corresponding `I: 'static` bound without restricting custom `PrimInt`
  implementations.
- If it cannot be removed without a larger compatibility or performance cost,
  document and test the bound as an intentional compatibility decision.

Pass condition: the bound is either eliminated or explicitly accepted with its
actual user impact documented.

## 4. Measure cached backend and type selection

Status: complete and rejected. A construction-time dispatch experiment stored a
typed mask function pointer in each `Lapper`, removing per-iterator backend
selection and per-mixed-block `TypeId` selection. On an Apple M3, paired
alternating measurements showed query regressions of 4.21% on `1-2`, 2.04% on
`7-3`, and 1.49% on `8-7`. The cached variant won only 1/15, 0/15, and 0/10
query pairs respectively. The indirect call costs more than the current
compiler-folded checks, so the existing dispatch is retained.

- Compare the current per-iterator backend detection and per-mask type
  selection with selection cached at `Lapper` construction.
- Measure construction, query, and total time on all three retained datasets.
- Reject function-pointer or cached-dispatch designs that are slower, more
  fragile, or materially more complex without a demonstrated benefit.

Pass condition: retain the fastest defensible dispatch arrangement based on
paired measurements, without adding a workload heuristic or user-visible mode.

## 5. Run the final CI target matrix

- Test AArch64 NEON, x86-64 AVX2, and x86-64 scalar execution.
- Compile and test the scalar fallback on other supported targets.
- Cover default features, `with_serde`, and `sort_unstable`.
- Include the declared MSRV and current stable Rust.

Pass condition: every supported target and feature combination is green, with
native execution for the SIMD backends claimed by the release.
