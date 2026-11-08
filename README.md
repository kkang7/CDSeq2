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

## Quick start: synthetic test dataset

The bundled [test-data generator](examples/01_simulate.R) creates a fully synthetic dataset with a fixed seed (`20261015`). It is included in the [example scripts](examples); no external data download is needed. The example runs deconvolution, trains the GAM6 scorer, evaluates it on separate validation donors, ranks gene–cell-type profiles, and demonstrates phenotype and genotype association analysis.

### Run the example

First install CDSeq2 using the R commands above. Then run the following in a terminal to obtain the example scripts and execute them:

```bash
git clone https://github.com/kkang7/CDSeq2.git
cd CDSeq2
Rscript --vanilla examples/run_tutorial.R
```

If you already have a copy of the repository, run only the final command from its root directory. Alternatively, download and extract the repository using GitHub's **Code → Download ZIP**, open a terminal in the extracted directory, and run that command.

Results are written to `examples/output/`. The script does not overwrite a nonempty output directory. For another run, supply a new directory, for example:

```bash
Rscript --vanilla examples/run_tutorial.R work/tutorial-second-run
```

### Inputs

The script generates these inputs automatically and saves them in the `data` component of `tutorial_objects.rds`:

| Input | Dimensions and meaning |
|---|---|
| Single-cell reference counts | 240 genes × 36 reference cells, with cell-type labels and six reference-donor identifiers |
| Training bulk RNA counts | 240 genes × 20 synthetic training donors |
| Validation bulk RNA counts | 240 genes × 20 separate synthetic validation donors |
| Ground-truth component counts | 240 genes × 3 cell types × 20 donors in each cohort, used for scorer training and evaluation |
| Phenotype and genotype | One continuous phenotype and one simulated genotype (0, 1 or 2) per donor in each cohort |

The three cell types are excitatory neurons, inhibitory neurons and astrocytes. Reference, training and validation donor identifiers are disjoint. Ground truth is available here because the example is simulated; it is not required to apply an already trained scorer to new bulk data.

### Expected outputs

A successful run prints `Tutorial completed:` followed by the output directory and creates these files:

| File in `examples/output/` | Contents |
|---|---|
| `tutorial_objects.rds` | Generated inputs, fitted models, inferred component counts (240 × 3 × 20 per cohort), RNA fractions (20 × 3 per cohort), and prioritization results |
| `toy_scorer.rds` | Fitted GAM6 scorer for this synthetic example |
| `profile_scores.csv` | All 720 validation gene–cell-type profiles, their scores, diagnostics and `retained` flags |
| `recovery_summary.csv` | Recovery correlations and retained counts by cell type; 48 of 240 profiles (20%) are retained in each type, for 144 total |
| `donor_separation.csv` | Checks that training and validation donors do not overlap |
| `phenotype_associations.csv` | All 720 phenotype tests, effect estimates, p values, BH-adjusted q values and calls before/after prioritization |
| `genotype_associations.csv` | Corresponding results for the simulated genotype |
| `sessionInfo.txt` | R version and package environment |
| `run_summary.txt` | Elapsed runtime |

Inspect the results in R, with the repository root as the working directory:

```r
x <- readRDS("examples/output/tutorial_objects.rds")
dim(x$data$validation$bulk)  # 240 20
dim(x$validation$counts)    # 240 3 20
dim(x$validation$fractions) # 20 3
read.csv("examples/output/recovery_summary.csv")
scores <- read.csv("examples/output/profile_scores.csv")
sum(scores$retained)        # 144
```

This small example uses short MCMC chains to demonstrate the workflow. Numerical results can vary across platforms and package versions, and the toy scorer is intended only for the synthetic example. See the [complete tutorial](https://kkang7.github.io/CDSeq2/chapters/quickstart.html) for further explanation.

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
