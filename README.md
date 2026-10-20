# CDSeq2

Archived software: [CDSeq2 0.1.0 (10.5281/zenodo.23489824)](https://doi.org/10.5281/zenodo.23489824).

Reference-based individual-level cell-type expression inference with optional gene-type prioritization.

This repository contains the R package, a Quarto book and reproduction workflows for the accompanying manuscript. The book covers installation, a complete synthetic tutorial, analysis of new data, scorer training and validation, phenotype/eQTL applications, all manuscript figures and underlying tables.

## Installation

Run the following in R:

```r
install.packages("remotes")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("limma", quietly = TRUE)) {
  BiocManager::install("limma", ask = FALSE, update = FALSE)
}
remotes::install_github("kkang7/CDSeq2")
```

R >= 4.1 and a working package compilation toolchain are required. The commands above install the Bioconductor dependency `limma` before installing CDSeq2.

## Documentation

Read the [CDSeq2 Quarto book](https://kkang7.github.io/CDSeq2/) for tutorials, analysis workflows and manuscript reproduction instructions.

## Reproduce manuscript version 16

```bash
python3 reproduction/reproduce.py tables
python3 reproduction/reproduce.py independent-prepare
python3 reproduction/reproduce.py brain-prepare
```

See `book/chapters/reproduce.qmd` for commands to rebuild all 25 figures and run the native independent eQTL benchmark. Statistical figures use the 42 included CSV tables. Figure 1 also requires Python with reportlab and pypdfium2. Verified reconstructed inputs and the source-verification workflow are included in the book downloads. Native reruns and the additional prioritization/null diagnostics are separate commands documented in the book.

The comparison is standard BayesPrism versus CDSeq2 without and with gene-type prioritization. Prioritization selects profiles from the same CDSeq2 fit; it does not change expression or association p values. Report accuracy and coverage together.

## Repository layout

| Directory | Contents |
|---|---|
| `R`, `src`, `man`, `data`, `inst` | R package implementation and help |
| `book` | Quarto sources, figures, table and benchmark downloads |
| `examples` | Four-part runnable synthetic vignette |
| `reproduction` | Paper and brain-study reproduction entry points |
| `tools` | Dependency setup, GitHub configuration and book checks |
| `.github/workflows` | GitHub Pages build and deployment |

`work/`, native fit caches, rendered website output and tutorial output are ignored by Git. Package installation excludes the book and reproduction assets. Original package authorship and the MIT license are preserved. See `RELEASE_NOTES.md` and `VALIDATION.txt` for the delivery scope.
