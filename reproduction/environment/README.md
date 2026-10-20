# Exact benchmark environment

`renv.lock` records R 4.6.1 and the installed benchmark dependency versions,
including the pinned BayesPrism commit. It excludes CDSeq2 deliberately: install
the exact archived 0.0.0.9000 source from the independent benchmark for primary
replay. Package 0.1.0 is the current user-facing release.

```bash
Rscript reproduction/environment/restore.R /absolute/path/to/benchmark-library
```

Set `R_LIBS_USER` to that library before invoking a native reproduction command.
The driver adds its archived CDSeq2 library in front of the provided dependency
library. Restore requires network access and suitable compilers/system libraries;
the lock does not supply R or guarantee a cross-platform bitwise binary match.
Some recorded Bioconductor packages were installed from R-universe repositories;
their version and source records are preserved rather than silently replaced.

The lock was created from the validated installed environment and every recorded
benchmark package version was compared with the original environment table.
A complete clean restoration on another operating system is not claimed here.
CI separately installs current compatible package dependencies for Linux/macOS
package checks; it is not a rerun of the full exact benchmark environment.
