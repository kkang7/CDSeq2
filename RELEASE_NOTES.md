# CDSeq2 0.1.0 — manuscript version 16

This release synchronizes the package, Quarto book, 42 supplementary
tables, 25 figures, reconstructed brain inputs and independent LD benchmark.
It adds a directly loadable frozen brain scorer with a model description,
`summarize_variation_selection()`, package tests, Linux/macOS CI, citation
metadata and a benchmark dependency lock. Inference/scoring functions are unchanged.

Additional analyses compare score groups and equally sized ranking/random subsets,
with paired cohort-bootstrap uncertainty, and diagnose continuous-null instability.
They reuse the completed v15 fits. All primary v15 numerical results remain intact.
Contribution ranking performs similarly to GAM6 in the principal eQTL setting.
Donor-38 sparse profiles drive the largest continuous-null call bursts; neither
another donor omission nor a new normalization is applied to primary analyses.

The v14 table archive is replaced by v16 material. Historical clustering audit
files remain identified as historical; current reproduction uses ReconstructedData.
See VALIDATION.txt for completed checks and outstanding external release steps.
