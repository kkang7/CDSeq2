> Historical input version: this clustering audit targets the pre-rerun inputs. Manuscript v16 uses the verified ReconstructedData download and source verifier described in the Quarto book. Do not use historical identity checks to certify the new input version.

# Reconstruct the ROSMAP manuscript inputs

Source data: [AD/Aging Brain Atlas](https://compbio.mit.edu/ad_aging_brain/), accompanying Mathys et al., Cell (2023), doi:10.1016/j.cell.2023.08.039. Download the nine processed Seurat RDS files and the deidentified participant metadata from that page. These inputs are not redistributed here. Follow the source site's access terms; complete raw data and additional participant metadata have separate Synapse/DUC requirements.

## Required downloads

Place these files together in an external DATA_DIR:

- Astrocytes.rds
- Excitatory_neurons_set1.rds
- Excitatory_neurons_set2.rds
- Excitatory_neurons_set3.rds
- Inhibitory_neurons.rds
- Oligodendrocytes.rds
- OPCs.rds
- Immune_cells.rds
- Vasculature_cells.rds
- individual_metadata_deidentified.tsv

The local source RDS files total approximately 174 GiB. Use a high-memory workstation or server and ample space for clustered checkpoints. The wrapper clusters one source component at a time and checkpoints completed components. Primary preparation still combines the components and requires substantial memory. The full replay was launched on a workstation with 128 GiB RAM. Do not place downloaded data or generated objects in Git.

## Dependencies

In R, install CDSeq2 as described in the book, then the preparation dependencies:

```r
install.packages(c("Seurat", "SeuratObject", "harmony", "data.table", "digest"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment", "SummarizedExperiment"),
                     ask = FALSE, update = FALSE)
```

Clustering can change with Seurat, Harmony, numerical libraries and source-data versions. This release includes the original preparation algorithms and frozen donor/gene selections, but does not claim that arbitrary dependency versions reproduce the paper. `original_reference_sessionInfo.txt` records an available historical construction environment; it is provenance, not a complete pinned environment for every stage.

## Run from the repository root

Use distinct DATA_DIR, WORK_DIR and OUTPUT_DIR paths. OUTPUT_DIR must be empty.

```bash
Rscript reproduction/brain/preprocessing/run.R verify-downloads /path/to/ROSMAP_AD /path/to/rosmap_work /path/to/rosmap_inputs
Rscript reproduction/brain/preprocessing/run.R cluster /path/to/ROSMAP_AD /path/to/rosmap_work /path/to/rosmap_inputs
Rscript reproduction/brain/preprocessing/run.R build /path/to/ROSMAP_AD /path/to/rosmap_work /path/to/rosmap_inputs
Rscript reproduction/brain/00_preflight.R /path/to/rosmap_inputs/brain_config.R
Rscript reproduction/brain/01_prepare.R /path/to/rosmap_inputs/brain_config.R
```

The optional `verify-downloads` stage reads all nine input files and compares their MD5 identities with the historical source files. It can take substantial I/O time. The metadata TSV is validated by the preparation engine; it has no frozen source checksum in this release.

Then run parts 02–07 from the brain-benchmark chapter with the generated `brain_config.R` as their argument. The config chooses a separate analysis output directory. The reconstructed request contains precisely the fields consumed by this public direct-aggregation workflow; it is not a replacement for unrelated historical fresh-donor simulation runners.

## What is reconstructed

1. The original Step 1 selects a balanced 200-donor pool, normalizes and clusters each source component using PCA/Harmony and the original settings. The primary training/reference preparation retains the top-ranked donor-shared cluster at resolution 0.01 per component (ranked by qualifying donor coverage, total nuclei, median qualifying donor count and cluster identifier), using the original primary engine.
2. `donors.csv` freezes the manuscript's disjoint 40 training, 20 reference and 40 evaluation donors in their original order; these are the source study's deidentified identifiers, not personal identifiers. `genes.txt` freezes the 2,000-gene panel. This is replay of a published design, not a new random donor split or gene-panel search.
3. Training counts are summed by donor and cell type. Reference cells are sampled with the historical seed 20260919 (the original function applies its internal offset) and historical caps, then restricted to the panel with zero-library cells excluded.
4. Evaluation counts are summed from all nuclei for the frozen evaluation donors in the original downloaded objects, without applying the training/reference cluster filter. The three excitatory-neuron source objects are combined. There is no count resampling or effect injection in this preparation step.
5. The script constructs the four input objects and checks their numerical identity against the manuscript inputs **before exporting them**. Expected identities include counts, sample/donor order, gene/type order, reference cell order and counts, cell counts, pooled reference, and inference settings.

RDS bytes can differ because of compression, serialization, or omitted historical path metadata. The downstream reader accepts either the original byte identities or the frozen numerical-content identities in `expected_content.csv`. It still checks dimensions, count conservation and separation of the 100 donors. A changed count or donor is rejected. Do not edit expected identities to force an altered experiment to pass.

If the numerical check fails, the script does not export a validated input set. Diagnostic candidate objects remain in WORK_DIR/candidate_inputs, and WORK_DIR/numerical_comparison.csv records the expected and observed signatures. These candidates are not verified manuscript inputs. Compare source versions and clustering environments with the historical configuration. Retain the original prepared inputs for exact replay until a matching reconstruction is achieved.

## Validation status

A complete source-data replay was executed on 18 October 2026. All nine source-file MD5 checks passed, all nine components completed clustering, and all four candidate inputs were assembled. **The final numerical identity check failed; this is not exact reproduction of the manuscript inputs.**

The 200-donor starting pool, frozen gene panel, donor identities and fitting settings match. Evaluation counts reproduce exactly (560,000 expression-count values and 280 nucleus counts). Training differs in 132,289 of 560,000 count values. The reconstructed reference contains 39,533 nuclei versus 39,692 archived nuclei, with 29,108 shared; counts for shared nuclei match exactly. Reference-dependent fields also cause the scaffold and request signatures to fail, despite exact evaluation counts and unchanged fitting settings.

The unresolved difference is in training/reference cell selection. The run used R 4.6.1, Seurat 5.5.1, SeuratObject 5.4.0 and Harmony 2.0.5 on a 128-GiB Apple arm64 workstation. Version-dependent clustering is a possible explanation; the original clustering environment is not fully documented, so this explanation is not proven. See [the complete audit](audit_2026-10-18/SUMMARY.txt), its CSV comparisons and session information.

The wrapper updates legacy Seurat object formats while checking exact preservation of raw counts and cell metadata. It saves component checkpoints and retains failed-build candidates for diagnosis, but exports verified manuscript inputs only after the frozen signatures pass. The original expected signatures remain unchanged. A synthetic contract test confirms rejection of altered counts, donor identities, reference data and fit settings.

Until the original cell-selection provenance or a matching environment is recovered and validated, use the archived prepared inputs for exact downstream replay. No manuscript inference results, figures or numerical tables were replaced by this audit.
