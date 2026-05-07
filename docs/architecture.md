# phase_tools-rs architecture

`phase_tools-rs` is being organized as a library-first Rust genomics package.
The binaries remain important user-facing entry points, but new algorithms and
shared I/O should live in the `phase_tools` library crate before being exposed
through CLI wrappers.

## Target shape

```text
src/lib.rs
src/assembly/          local assembly and assembly-backed adjudication
src/io/                FASTA/VCF/BAM/CRAM/TSV helpers
src/variant/           alleles, genotypes, phase tags, normalization helpers
src/phase/             read-backed phasing, read selection, MEC/greedy kernels
src/mnv/               MNV/COMPLEX construction and output helpers
src/qc/                BAM/CRAM error, contamination, and ancestry kernels
src/commands/          CLI adapters around library functions
src/bin/               minimal binary entry points
```

The first library boundaries are now in place for fermi-lite assembly and FASTA
reference access:

```text
phase_tools::assembly::fermi_lite
phase_tools::io::fasta
```

`fermi_lite_assemble` and `phase_adjudicate` call the assembly module instead
of path-including assembly code. `phase_adjudicate`, `bam_contamination`, and
`bam_ancestry` share the FASTA/FAI wrapper instead of each owning a separate
htslib `faidx` wrapper.

## Refactor rules

1. Keep behavior-preserving extraction separate from feature changes.
2. Move duplicated primitives into the library before adding new binary options.
3. Keep `rust-htslib` as the only htslib access path in Rust code.
4. Preserve CLI output formats and test fixtures while moving internals.
5. Expose narrow, typed kernels first; keep CLI parsing and printing at the
   command boundary.
6. Avoid adding dependencies unless they materially improve a reusable library
   module.

## Immediate extraction sequence

1. Shared error/result helpers.
2. Complete migration to the shared FASTA/FAI wrapper, including `phase_mnv_rs`
   and `bam_error_model`.
3. BAM record filtering and base/event extraction.
4. VCF/BCF genotype, `PS`, and `HP` parsing helpers.
5. `phase_compare` comparison kernel.
6. BAM phasing read collection, read selection, MEC, and greedy kernels.
7. MNV/COMPLEX observation collection, construction, and writing.
8. BAM/CRAM QC kernels for error, contamination, and ancestry tools.

## Current caveats

The main `phase_mnv_rs` implementation still contains most phasing and MNV logic
in `src/main.rs`. That is intentional during the transition: extraction should
be incremental and validated after each step rather than a large rewrite that
changes behavior and architecture at the same time.
