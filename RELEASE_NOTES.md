# Documentation delivery — 15 October 2026

- Added a Quarto book in the same repository as the existing R package.
- Added a fully synthetic four-part tutorial using the public CDSeq2 variation API.
- Included revision-14 figure assets, complete captions, 35 unchanged supplementary CSVs and explicit figure-to-table links.
- Added one entry point for table extraction, all-figure regeneration and independent-cohort benchmark execution.
- Preserved the frozen independent benchmark and exact package source snapshot in a downloadable archive.
- Added current standard BayesPrism entry points for the original brain-data workflow. Its native input preprocessing is applied and no output selection is used for BayesPrism.
- Added GitHub Pages deployment, local checks and one-command repository-owner configuration.
- Preserved existing package source code, authorship and license. The minimum R version is 4.1 because the existing implementation uses native pipes; documentation/build exclusions were added.

The paper's historical numerical files retain machine labels for reproducibility. Reader-facing text uses full cell-type names and the distinction between CDSeq2 with and without gene-type prioritization. Model development provenance, the post hoc 39-donor continuous sensitivity and accuracy/coverage tradeoffs remain documented.

The shared source directory was not modified. The book does not assign a final manuscript DOI or invent a GitHub owner. Original study input release remains outstanding. Native experiment runs completed earlier are distinguished from this delivery's tutorial and saved-data validation.
