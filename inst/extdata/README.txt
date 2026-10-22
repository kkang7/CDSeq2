Frozen DEG scorer storage

The core, annotations and metadata RDS files together contain the original
deg_scorer_v0_23.rds object from CDSeq2 v0.1.0. Storage is split to keep each
file below 5 MB; fitted models and prediction rules are unchanged.
load_deg_scorer() reconstructs the object and replays its numerical controls.
An explicit path still accepts the original complete model RDS.
The example RDS contains the original teaching inputs.
Training context and original source checksums are in the model metadata.
Source release: https://doi.org/10.5281/zenodo.23489824
