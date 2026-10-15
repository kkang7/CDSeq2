# CDSeq2

Reference-based individual-level cell-type expression inference with optional gene-type prioritization.

This repository contains the R package, a Quarto book and reproduction workflows for the accompanying manuscript. The book covers installation, a complete synthetic tutorial, analysis of new data, scorer training and validation, phenotype/eQTL applications, all manuscript figures and underlying tables.

## Install and try

```bash
Rscript tools/install_dependencies.R
R CMD INSTALL .
Rscript --vanilla examples/run_tutorial.R
```

The example generates its own data and writes `examples/output/`. Its four scripts can also be read separately. R >= 4.1 and a working package compilation toolchain are required.

## Read and publish the book

```bash
quarto preview book
```

Before publishing, configure your repository owner:

```bash
python3 tools/configure_github.py YOUR_GITHUB_USERNAME --repo CDSeq2
```

Push the repository to GitHub and set **Settings → Pages → Source → GitHub Actions**. The included **Publish Quarto book** workflow renders and deploys the website on pushes to `main`, or when run manually. See `book/chapters/publishing.qmd` for the complete steps. No repository or website has been created by this local delivery.

## Reproduce the manuscript

```bash
python3 reproduction/reproduce.py tables
python3 reproduction/reproduce.py independent-prepare
```

See `book/chapters/reproduce.qmd` for commands to rebuild all 21 figures and run the native independent eQTL benchmark. Statistical figures use the 35 included CSV tables. Figure 1 also requires Python with reportlab and pypdfium2. The original brain-data analyses require prepared study inputs that are listed in `reproduction/brain/DATA_REQUIREMENTS.txt` and are not included.

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
