# CDSeq2: Bioconductor submission candidate

This branch prepares version 0.99.0 for Bioconductor review. It is not yet
accepted by or installable from Bioconductor. It requires R 4.6 or later.

CDSeq2 infers individual-level cell-type expression and RNA fractions from
bulk RNA sequencing with a single-cell reference. Optional gene-type
prioritization helps focus follow-up on supported profiles. Scorer training
requires known component expression and independent validation; applying a
trained scorer does not require truth. Scores are not association p-values.

## Documentation

Run `vignette("CDSeq2", package = "CDSeq2")` after installing with vignettes.
The executable vignette generates synthetic reference, training and validation
data and demonstrates inference, scoring and selection diagnostics.

Full documentation: https://kkang7.github.io/CDSeq2/
Frozen paper release: https://doi.org/10.5281/zenodo.23489824
Paper reproduction workflows remain on the GitHub main branch and in the archive.
The archived version 0.1.0 is unchanged.

## Local installation

From the parent directory, run `R CMD INSTALL CDSeq2_0.99.0.tar.gz` after installing
the dependencies listed in DESCRIPTION. The Bioconductor development check
workflow covers Linux, macOS and Windows when this branch is pushed.

## Changes for review

See NEWS. Public function names are preserved. Explicitly seeded operations
now restore the caller RNG state. Frozen DEG model storage is split into smaller
files without refitting; the loader still accepts the original complete RDS.
Benchmark gene filtering accepts limma and DESeq2-style columns and correctly
labels shared effects only in the affected types.

## Attribution

Kai Kang developed CDSeq2; Kaiyu Xie contributed benchmarking. Supported by the
AMS-Simons Research Enhancement Grant for PUI Faculty and the UNCW CSE Research
Initiative fund. See inst/AI_ASSISTANCE.txt for AI assistance disclosure.
