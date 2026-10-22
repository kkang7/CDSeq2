


#' Create cell-type-specific pseudo-bulk profiles at the individual level
#'
#' This function aggregates single-cell RNA-seq count data into
#' cell-type-specific pseudo-bulk expression profiles for each individual
#' by summing counts across all cells belonging to the same
#' (cell type, patient) pair.
#'
#' The output is a three-dimensional array with dimensions:
#' genes x cell types x patients, which can be used to define
#' ground-truth cell-type-specific gene expression and differential
#' expression signals for benchmarking bulk RNA-seq deconvolution
#' methods such as CDSeq2.
#'
#' @param counts A numeric matrix or sparse matrix (e.g. \code{dgCMatrix})
#'   of raw counts with genes as rows and cells as columns.
#'   Row names must correspond to gene identifiers.
#'
#' @param cell_type A character or factor vector of length equal to
#'   \code{ncol(counts)}, indicating the annotated cell type for each cell.
#'
#' @param patient_id A character or factor vector of length equal to
#'   \code{ncol(counts)}, indicating the individual (patient / subject)
#'   from which each cell was obtained.
#'
#' @param drop_empty Logical; whether to remove cell types or patients
#'   that have zero total counts across all genes after aggregation.
#'   Default is \code{TRUE}.
#'
#' @param return_long Logical; if \code{FALSE} (default), the function
#'   returns a 3D array of dimensions genes x cell types x patients.
#'   If \code{TRUE}, returns a long-format data frame with columns
#'   \code{gene}, \code{cell_type}, \code{patient}, and \code{count}.
#'
#' @return
#' If \code{return_long = FALSE}, a numeric array with dimensions
#' \code{genes x cell_types x patients}, where each entry represents
#' the summed expression of a gene within a given cell type for a given
#' individual.
#'
#' If \code{return_long = TRUE}, a data frame in long format suitable
#' for downstream differential expression analysis using packages such
#' as \code{limma}, \code{edgeR}, or \code{DESeq2}.
#'
#' @details
#' This function performs simple summation of raw counts and does not
#' apply normalization, scaling, or filtering based on cell counts.
#' These steps should be performed downstream depending on the
#' differential expression framework used.
#'
#' @importFrom Matrix rowSums
#' @export
#' @examples
#' data(cdseq2_toy)
#' x <- cdseq2_toy
#' counts <- x$data$reference
#' ct <- x$data$annotation
#' donors <- x$data$reference_donors
#' pb <- create_ct_pseudobulk(counts, ct, donors)
#' dim(pb)
create_ct_pseudobulk <- function(
    counts,
    cell_type,
    patient_id,
    drop_empty = TRUE,
    return_long = FALSE
) {
  stopifnot(
    ncol(counts) == length(cell_type),
    ncol(counts) == length(patient_id)
  )

  genes <- rownames(counts)
  cell_types <- unique(cell_type)
  patients <- unique(patient_id)

  pb_array <- array(
    0,
    dim = c(length(genes), length(cell_types), length(patients)),
    dimnames = list(
      gene = genes,
      cell_type = cell_types,
      patient = patients
    )
  )

  for (p in patients) {
    idx_p <- which(patient_id == p)
    for (ct in cell_types) {
      idx_ct <- idx_p[cell_type[idx_p] == ct]
      if (length(idx_ct) == 0) next
      pb_array[, ct, p] <- Matrix::rowSums(counts[, idx_ct, drop = FALSE])
    }
  }

  if (drop_empty) {
    keep <- apply(pb_array, c(2, 3), sum) > 0
    keep_ct <- apply(keep, 1, any)
    keep_p  <- apply(keep, 2, any)
    pb_array <- pb_array[, keep_ct, keep_p, drop = FALSE]
  }

  if (!return_long) return(pb_array)

  df <- as.data.frame(as.table(pb_array))
  colnames(df) <- c("gene", "cell_type", "patient", "count")
  df
}


#' Cell-Type-Specific Differential Expression from Pseudo-Bulk Data
#'
#' Perform cell-type-specific differential expression (DE) analysis using
#' a pseudo-bulk strategy and the limma--voom framework. The function supports
#' both unpaired and paired (repeated-measures) experimental designs.
#'
#' For each cell type, single-cell counts are aggregated to the patient level,
#' forming a gene-by-patient pseudo-bulk matrix. Differential expression is
#' then tested across conditions using limma.
#'
#' @param counts A numeric matrix or sparse matrix of raw counts with genes
#'   in rows and cells in columns.
#'
#' @param cell_type A character or factor vector of length \code{ncol(counts)}
#'   specifying the (major) cell type for each cell.
#'
#' @param patient_id A vector of length \code{ncol(counts)} specifying the
#'   patient or subject identifier for each cell. This defines the
#'   pseudo-bulk aggregation unit.
#'
#' @param condition A character or factor vector of length ncol(counts),
#'   giving the condition of each cell in column order. In paired designs,
#'   each patient contributes cells from both conditions.
#'
#' @param min_cells Integer. Minimum number of cells required for a given
#'   cell type within a patient for that patient to be included in DE
#'   analysis. Default is 10.
#'
#' @param paired Logical. If \code{FALSE} (default), an unpaired design is used
#'   with model \code{~ condition}. If \code{TRUE}, a paired design is used with
#'   model \code{~ patient_id + condition}, appropriate for repeated-measures
#'   or longitudinal experiments.
#'
#' @param design_formula Optional model formula. If \code{NULL}, the formula is
#'   chosen automatically based on \code{paired}. Advanced users may override
#'   this to specify custom designs.
#'
#' @param contrast Optional contrast string passed to
#'   \code{limma::makeContrasts} for hypothesis testing. If \code{NULL}, the
#'   default coefficient corresponding to \code{condition} is tested.
#'
#' @param p_adj_method Character string specifying the multiple-testing
#'   adjustment method passed to \code{limma::topTable}. Default is \code{"BH"}.
#'
#' @return A named list of data frames, one per cell type. Each data frame
#'   contains limma differential expression results with at least the following
#'   columns:
#'   \itemize{
#'     \item \code{log2FoldChange}: estimated log-fold change,
#'     \item \code{AveExpr}: average expression,
#'     \item \code{t}: moderated t-statistic,
#'     \item \code{P.Value}: raw p-value,
#'     \item \code{padj}: adjusted p-value,
#'     \item \code{gene}: gene identifier.
#'   }
#'
#' @details
#' \strong{Unpaired design:}
#' Each patient appears in exactly one condition. Differential expression is
#' tested across patients using the model \code{~ condition}.
#'
#' \strong{Paired design:}
#' Each patient appears in multiple conditions (e.g. before/after treatment).
#' Patient-specific baseline effects are accounted for using the model
#' \code{~ patient_id + condition}. This is the recommended mode for datasets
#' such as the Kang et al. IFN-\eqn{\beta} PBMC study.
#'
#' The function automatically skips cell types for which the design matrix
#' is not full rank or for which insufficient patient-level replication is
#' available.
#'
#' @importFrom Matrix rowSums
#' @importFrom limma voom lmFit eBayes topTable makeContrasts contrasts.fit
#' @importFrom stats model.matrix
#' @export
#' @examples
#' data(cdseq2_toy)
#' x <- cdseq2_toy
#' counts <- x$data$reference
#' ct <- x$data$annotation
#' donors <- x$data$reference_donors
#' condition <- setNames(rep(c("control", "treated"), each = 3),
#'     sort(unique(donors)))
#' result <- call_ct_DEGs_pb(counts, ct, donors,
#'     unname(condition[donors]), min_cells = 1, paired = FALSE)
#' names(result)
call_ct_DEGs_pb <- function(
    counts,
    cell_type,
    patient_id,
    condition,
    min_cells = 10,
    paired = FALSE,
    design_formula = NULL,
    contrast = NULL,
    p_adj_method = "BH"
) {

  stopifnot(
    ncol(counts) == length(cell_type),
    ncol(counts) == length(patient_id),
    ncol(counts) == length(condition)
  )

  if (is.null(rownames(counts))) stop("counts must have rownames (gene IDs).")

  genes <- rownames(counts)
  cell_types <- unique(as.character(cell_type))

  ## cell-level metadata
  meta_cell <- data.frame(
    cell_type   = as.character(cell_type),
    patient_id  = as.character(patient_id),
    condition   = as.character(condition),
    stringsAsFactors = FALSE
  )

  ## default formulas
  if (is.null(design_formula)) {
    design_formula <- if (paired) ~ patient_id + condition else ~ condition
  }

  results <- vector("list", length(cell_types))
  names(results) <- cell_types

  for (ct in cell_types) {

    ct_cells <- which(meta_cell$cell_type == ct)
    if (length(ct_cells) == 0) next

    meta_ct <- meta_cell[ct_cells, , drop = FALSE]
    counts_ct <- counts[, ct_cells, drop = FALSE]

    ## ----------------------------
    ## Determine eligible patients
    ## ----------------------------
    if (!paired) {
      ## unpaired: require >= min_cells per patient (within this cell type)
      tab_p <- table(meta_ct$patient_id)
      keep_patients <- names(tab_p)[tab_p >= min_cells]
      if (length(keep_patients) < 2) next

      ## patient-level condition (must be unique per patient in unpaired setting)
      cond_by_patient <- tapply(meta_ct$condition, meta_ct$patient_id, function(x) unique(x))
      bad <- names(cond_by_patient)[lengths(cond_by_patient) != 1]
      if (length(bad) > 0) {
        stop(sprintf(
          "Unpaired mode expects one condition per patient. Patients with multiple conditions in cell type %s: %s",
          ct, paste(bad, collapse = ", ")
        ))
      }

      keep_patients <- intersect(keep_patients, names(cond_by_patient))
      if (length(keep_patients) < 2) next

      ## ----------------------------
      ## Pseudobulk: genes × patient
      ## ----------------------------
      pb_mat <- sapply(keep_patients, function(p) {
        idx <- which(meta_ct$patient_id == p)
        Matrix::rowSums(counts_ct[, idx, drop = FALSE])
      })
      colnames(pb_mat) <- keep_patients
      rownames(pb_mat) <- genes

      meta_pb <- data.frame(
        patient_id = keep_patients,
        condition  = unname(sapply(keep_patients, function(p) cond_by_patient[[p]][1])),
        stringsAsFactors = FALSE
      )

    } else {
      ## paired: require >= min_cells in EACH condition for a patient (within this cell type)
      cond_levels <- sort(unique(meta_ct$condition))
      tab_pc <- table(meta_ct$patient_id, meta_ct$condition)

      keep_patients <- rownames(tab_pc)[apply(tab_pc, 1, function(v) all(v[cond_levels] >= min_cells))]
      if (length(keep_patients) < 2) next

      ## build patient×condition sample keys
      meta_ct2 <- meta_ct[meta_ct$patient_id %in% keep_patients, , drop = FALSE]
      counts_ct2 <- counts_ct[, meta_ct$patient_id %in% keep_patients, drop = FALSE]

      sample_key <- paste(meta_ct2$patient_id, meta_ct2$condition, sep = "__")
      sample_levels <- sort(unique(sample_key))

      ## ----------------------------
      ## Pseudobulk: genes × (patient×condition)
      ## ----------------------------
      pb_mat <- sapply(sample_levels, function(sk) {
        idx <- which(sample_key == sk)
        Matrix::rowSums(counts_ct2[, idx, drop = FALSE])
      })
      colnames(pb_mat) <- sample_levels
      rownames(pb_mat) <- genes

      meta_pb <- data.frame(
        patient_id = sub("__.*$", "", sample_levels),
        condition  = sub("^.*__", "", sample_levels),
        stringsAsFactors = FALSE
      )
    }

    ## factors + order alignment
    meta_pb$patient_id <- factor(meta_pb$patient_id)
    meta_pb$condition  <- factor(meta_pb$condition)

    ## ----------------------------
    ## Design matrix (MUST match pb_mat columns)
    ## ----------------------------
    design <- model.matrix(design_formula, data = meta_pb)

    if (nrow(design) != ncol(pb_mat)) {
      stop(sprintf(
        "Design/pseudobulk mismatch for cell type %s: nrow(design)=%d, ncol(pb_mat)=%d",
        ct, nrow(design), ncol(pb_mat)
      ))
    }

    if (qr(design)$rank < ncol(design)) {
      warning(sprintf("Design matrix not full rank for cell type %s; skipping.", ct))
      next
    }

    ## ----------------------------
    ## limma-voom
    ## ----------------------------
    v <- limma::voom(pb_mat, design, plot = FALSE)
    fit <- limma::lmFit(v, design)

    if (!is.null(contrast)) {
      contr <- limma::makeContrasts(contrasts = contrast, levels = design)
      fit <- limma::contrasts.fit(fit, contr)
      fit <- limma::eBayes(fit)

      tt <- limma::topTable(
        fit,
        coef = 1,
        number = Inf,
        adjust.method = p_adj_method,
        sort.by = "none"
      )
    } else {
      fit <- limma::eBayes(fit)

      ## choose the (first) condition coefficient
      cond_coef <- grep("^condition", colnames(design))
      if (length(cond_coef) == 0) {
        warning(sprintf("No condition coefficient found for cell type %s; skipping.", ct))
        next
      }

      tt <- limma::topTable(
        fit,
        coef = cond_coef[1],
        number = Inf,
        adjust.method = p_adj_method,
        sort.by = "none"
      )
    }

    tt$gene <- rownames(tt)
    results[[ct]] <- tt
  }

  results
}







#' Select Benchmark Genes Based on Cell-Type-Specific Differential Expression
#'
#' This function constructs a gene panel for benchmarking deconvolution or
#' cell-type-specific differential expression methods. It uses pseudo-bulk
#' differential expression (DE) results to identify:
#' \itemize{
#'   \item cell-type-specific DEGs,
#'   \item shared DEGs across multiple cell types,
#'   \item and non-DEGs.
#' }
#'
#' Genes are sampled according to user-specified targets to form a controlled
#' benchmarking gene set, along with a label table indicating the DE status
#' of each selected gene.
#'
#' @param ct_deg_list A named list of data frames, one per cell type.
#'   Each data frame must contain at least the columns:
#'   \code{gene}, \code{log2FoldChange}, and \code{adj.P.Val}.
#'   Names of the list correspond to cell type identifiers.
#'
#' @param target_cell_types Character vector of cell types to include in the
#'   benchmark. Must be a subset of \code{names(ct_deg_list)}.
#'
#' @param all_genes Optional character vector specifying the full gene universe.
#'   If \code{NULL} (default), the gene universe is defined as the union of
#'   genes appearing in \code{ct_deg_list}.
#'
#' @param n_deg_per_ct Integer. Number of cell-type-specific DEGs to sample
#'   for each target cell type.
#'
#' @param n_shared_deg Integer. Number of shared DEGs (significant in two or
#'   more cell types) to sample. Set to 0 to exclude shared DEGs.
#'
#' @param n_non_deg Integer. Number of non-DEGs (not significant in any
#'   target cell type) to sample.
#'
#' @param deg_pval_cutoff Numeric. Adjusted p-value threshold used to define
#'   DEGs. Default is 0.05.
#'
#' @param deg_logfc_cutoff Numeric. Absolute log-fold-change threshold used
#'   to define DEGs. Default is 0.5.
#'
#' @param seed Optional integer. Random seed for reproducible gene sampling.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{genes}{Character vector of selected benchmark genes.}
#'   \item{labels}{Data frame with one row per selected gene, containing:
#'     \itemize{
#'       \item \code{gene}: gene identifier,
#'       \item \code{DE_<cell_type>}: logical indicators of DE status per cell type,
#'       \item \code{DE_class}: categorical label with values
#'         \code{"non-DEG"}, \code{"cell-type-specific DEG"}, or \code{"shared DEG"}.
#'     }}
#'   \item{summary}{List summarizing the number of selected genes in each DE class.}
#' }
#'
#' @details
#' A gene is classified as:
#' \itemize{
#'   \item \strong{cell-type-specific DEG} if it is significant in exactly one
#'   target cell type,
#'   \item \strong{shared DEG} if it is significant in two or more target cell types,
#'   \item \strong{non-DEG} if it is not significant in any target cell type.
#' }
#'
#' This function is typically used after pseudo-bulk DE analysis
#' (e.g. via \code{call_ct_DEGs_pb}) and before synthetic data generation
#' or benchmarking evaluation.
#'
#' @seealso
#' \code{\link{call_ct_DEGs_pb}}
#'
#' @export

#' @examples
#' genes <- paste0("g", seq_len(20))
#' a <- data.frame(gene = genes, log2FoldChange = c(rep(1, 5), rep(0, 15)),
#'     adj.P.Val = c(rep(0.001, 5), rep(0.9, 15)))
#' b <- a
#' b$log2FoldChange <- rev(a$log2FoldChange)
#' b$adj.P.Val <- rev(a$adj.P.Val)
#' selected <- filter_genes_benchmark(list(A = a, B = b), c("A", "B"),
#'     n_deg_per_ct = 2, n_non_deg = 3, seed = 1)
#' selected$genes
filter_genes_benchmark <- function(
    ct_deg_list,
    target_cell_types,
    all_genes = NULL,
    n_deg_per_ct = 200,
    n_shared_deg = 0,
    n_non_deg = 500,
    deg_pval_cutoff = 0.05,
    deg_logfc_cutoff = 0.25,
    seed = NULL
) {

  if (!is.null(seed)) withr::local_seed(seed)

  stopifnot(all(target_cell_types %in% names(ct_deg_list)))

  ## ---- enforce absolute logFC cutoff ----
  deg_logfc_cutoff <- abs(deg_logfc_cutoff)

  ## ---- define gene universe ----
  if (is.null(all_genes)) {
    all_genes <- Reduce(
      union,
      lapply(ct_deg_list, function(df) df$gene)
    )
  }

  ## ---- build DEG indicator matrix ----
  # Assisted-by: OpenAI Codex. Accept the documented limma column names
  # as well as legacy DESeq2-style tables; missing values are not discoveries.
  stopifnot(length(target_cell_types) > 0L, length(all_genes) > 0L)
  deg_mat <- vapply(target_cell_types, function(ct) {
    df <- ct_deg_list[[ct]]
    p_col <- if ("padj" %in% names(df)) "padj" else "adj.P.Val"
    fc_col <- if ("log2FoldChange" %in% names(df)) "log2FoldChange" else "logFC"
    if (!all(c("gene", p_col, fc_col) %in% names(df)))
      stop("Each DEG table needs gene, adjusted p-value and log-fold-change columns.")
    keep <- is.finite(df[[p_col]]) & is.finite(df[[fc_col]]) &
      df[[p_col]] < deg_pval_cutoff & abs(df[[fc_col]]) >= deg_logfc_cutoff
    all_genes %in% df$gene[keep]
  }, logical(length(all_genes)))
  deg_mat <- matrix(deg_mat, nrow = length(all_genes),
    dimnames = list(all_genes, target_cell_types))

  rownames(deg_mat) <- all_genes
  deg_count <- rowSums(deg_mat)

  ## ---- define gene categories ----
  ct_specific <- lapply(seq_along(target_cell_types), function(i) {
    rownames(deg_mat)[deg_mat[, i] & deg_count == 1]
  })
  names(ct_specific) <- target_cell_types

  shared_deg <- rownames(deg_mat)[deg_count >= 2]
  non_deg <- rownames(deg_mat)[deg_count == 0]

  ## ---- sample genes ----
  sel_genes <- character(0)

  for (ct in target_cell_types) {
    g <- ct_specific[[ct]]
    sel_genes <- c(
      sel_genes,
      if (length(g) <= n_deg_per_ct) g else sample(g, n_deg_per_ct)
    )
  }

  if (n_shared_deg > 0 && length(shared_deg) > 0) {
    sel_genes <- c(
      sel_genes,
      if (length(shared_deg) <= n_shared_deg) shared_deg else sample(shared_deg, n_shared_deg)
    )
  }

  if (length(non_deg) > 0) {
    sel_genes <- c(
      sel_genes,
      if (length(non_deg) <= n_non_deg) non_deg else sample(non_deg, n_non_deg)
    )
  }

  sel_genes <- sort(unique(sel_genes))

  ## ---- label table ----
  label_df <- data.frame(gene = sel_genes, stringsAsFactors = FALSE)

  for (ct in target_cell_types) {
    label_df[[paste0("DE_", ct)]] <-
      deg_mat[match(label_df$gene, all_genes), ct]
  }

  label_df$DE_class <- ifelse(
    rowSums(label_df[ , grepl("^DE_", names(label_df)), drop = FALSE]) == 0,
    "non-DEG",
    ifelse(
      rowSums(label_df[ , grepl("^DE_", names(label_df)), drop = FALSE]) == 1,
      "cell-type-specific DEG",
      "shared DEG"
    )
  )

  list(
    genes  = sel_genes,
    labels = label_df,
    summary = list(
      n_total = length(sel_genes),
      n_ct_specific = sum(label_df$DE_class == "cell-type-specific DEG"),
      n_shared = sum(label_df$DE_class == "shared DEG"),
      n_non_deg = sum(label_df$DE_class == "non-DEG")
    )
  )
}



#' Filter single-cell count data by cell-type support and gene expression
#'
#' This function performs **support-aware filtering at the single-cell level**
#' to remove poorly supported cell types and low-information genes *before*
#' pseudobulk aggregation or gene expression profile (GEP) estimation.
#'
#' Filtering at the raw single-cell stage is statistically preferable to
#' filtering pseudobulk data, as it allows exact assessment of cell-type
#' support across samples and prevents unstable GEP estimates caused by
#' sparsely represented cell types (e.g. 1–2 cells in isolated samples).
#'
#' ## Filtering strategy
#'
#' The function applies filtering in two main stages:
#'
#' 1. **Global cell-type filtering**
#'
#' For each cell type \eqn{k}, the number of samples in which it is represented
#' by at least \code{min_cells_per_ct} cells is computed. A cell type is retained
#' **globally** only if:
#'
#' \deqn{
#' \sum_s I\{n_{k,s} \ge \text{min\_cells\_per\_ct}\}
#' \ge \text{min\_samples\_per\_ct}
#' }
#'
#' Cell types failing this criterion are removed from **all samples**.
#'
#' 2. **Gene filtering**
#'
#' Genes are retained if they satisfy both:
#' \itemize{
#'   \item A minimum total read count across all retained cells
#'   \item Expression in at least a minimum number of cells
#' }
#'
#' ## Rationale
#'
#' Cell types with extremely low representation yield non-identifiable and
#' high-variance gene expression profile estimates. Removing such cell types
#' globally enforces minimal data support conditions required for stable
#' deconvolution, reference construction, and benchmarking.
#'
#' This procedure is **not heuristic**; it enforces fundamental statistical
#' identifiability constraints.
#'
#' @param counts A numeric matrix or \code{dgCMatrix} of dimension
#'   genes x cells containing raw single-cell read counts.
#'
#' @param cell_type A character or factor vector of length \code{ncol(counts)}
#'   giving the cell-type annotation for each cell.
#'
#' @param sample_id A character or factor vector of length \code{ncol(counts)}
#'   indicating the sample / patient / individual ID for each cell.
#'
#' @param min_cells_per_ct Integer at least 1.
#'   Minimum number of cells required for a cell type to be considered present
#'   within a sample.
#'
#' @param min_samples_per_ct Integer at least 1.
#'   Minimum number of samples in which a cell type must meet
#'   \code{min_cells_per_ct} to be retained globally.
#'
#' @param min_gene_total Integer at least 0.
#'   Minimum total read count across all retained cells required for a gene
#'   to be kept.
#'
#' @param min_cells_per_gene Integer at least 1.
#'   Minimum number of cells in which a gene must be expressed (nonzero)
#'   to be retained.
#'
#' @param verbose Logical.
#'   If \code{TRUE}, prints a summary of filtering results.
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item \code{counts} Filtered single-cell count matrix (genes x cells)
#'   \item \code{cell_type} Filtered cell-type annotation vector
#'   \item \code{sample_id} Filtered sample ID vector
#'   \item \code{keep_genes} Logical vector indicating retained genes
#'   \item \code{keep_celltypes} Logical vector indicating retained cell types
#'   \item \code{report} A list summarizing filtering statistics
#' }
#'
#' @details
#' The output of this function is guaranteed to be internally consistent:
#' \itemize{
#'   \item \code{counts}, \code{cell_type}, and \code{sample_id} are aligned
#'   \item No missing cell-type or sample annotations remain
#'   \item All retained cell types satisfy global support criteria
#' }
#'
#' The filtered output can be safely used for:
#' \itemize{
#'   \item Pseudobulk construction
#'   \item Reference GEP estimation
#'   \item Deconvolution benchmarking
#'   \item Spatial transcriptomics simulations
#' }
#'
#' @author Kai Kang
#'
#' @export

#' @examples
#' data(cdseq2_toy)
#' x <- cdseq2_toy
#' counts <- x$data$reference
#' ct <- x$data$annotation
#' donors <- x$data$reference_donors
#' filtered <- filter_sc_counts(counts, ct, donors,
#'     min_cells_per_ct = 1, min_samples_per_ct = 1,
#'     min_gene_total = 1, min_cells_per_gene = 1)
#' dim(filtered$counts)
filter_sc_counts <- function(
    counts,
    cell_type,
    sample_id,
    min_cells_per_ct = 5,
    min_samples_per_ct = 3,
    min_gene_total = 10,
    min_cells_per_gene = 3,
    verbose = TRUE
) {

  ## ============================
  ## 0. Basic checks
  ## ============================
  stopifnot(
    is.matrix(counts) || inherits(counts, "Matrix"),
    length(cell_type) == ncol(counts),
    length(sample_id) == ncol(counts)
  )

  if (inherits(counts, "Matrix") && !inherits(counts, "CsparseMatrix")) {
    counts <- as(counts, "CsparseMatrix")
  }


  ## Remove NA cells immediately
  keep_cells <- !is.na(cell_type) & !is.na(sample_id)
  counts <- counts[, keep_cells, drop = FALSE]
  cell_type <- cell_type[keep_cells]
  sample_id <- sample_id[keep_cells]

  ## ============================
  ## 1. Global cell-type support
  ## ============================
  ct_sample_table <- table(cell_type, sample_id)

  ct_supported_samples <- rowSums(ct_sample_table >= min_cells_per_ct)

  keep_celltypes <- ct_supported_samples >= min_samples_per_ct
  names(keep_celltypes) <- rownames(ct_sample_table)

  if (verbose) {
    message("Cell types kept: ",
            sum(keep_celltypes), "/", length(keep_celltypes))
  }

  ## Drop cells from unsupported cell types
  keep_cells <- keep_celltypes[cell_type]
  counts <- counts[, keep_cells, drop = FALSE]
  cell_type <- cell_type[keep_cells]
  sample_id <- sample_id[keep_cells]

  ## ============================
  ## 2. Gene filtering
  ## ============================
  gene_total <- Matrix::rowSums(counts)
  gene_cell_support <- Matrix::rowSums(counts > 0)

  keep_genes <- gene_total >= min_gene_total &
    gene_cell_support >= min_cells_per_gene

  counts <- counts[keep_genes, , drop = FALSE]

  if (verbose) {
    message("Genes kept: ",
            sum(keep_genes), "/", length(keep_genes))
  }

  ## ============================
  ## 3. Final integrity checks
  ## ============================
  stopifnot(
    ncol(counts) == length(cell_type),
    ncol(counts) == length(sample_id),
    !anyNA(cell_type),
    !anyNA(sample_id)
  )

  ## ============================
  ## 4. Report
  ## ============================
  report <- list(
    cells_after = ncol(counts),
    genes_after = nrow(counts),
    celltypes_after = length(unique(cell_type)),
    samples_after = length(unique(sample_id))
  )

  if (verbose) {
    message("filter_sc_counts summary:")
    print(report)
  }

  ## ============================
  ## 5. Return
  ## ============================
  return(list(
    counts = counts,
    cell_type = cell_type,
    sample_id = sample_id,
    keep_genes = keep_genes,
    keep_celltypes = keep_celltypes,
    report = report
  ))
}



#' Generate an aligned single-cell and pseudobulk benchmark dataset
#'
#' @description
#' Constructs a benchmarking dataset by aligning PB single-cell and reference
#' single-cell data, filtering genes and cell types based on PB support,
#' optionally enforcing minimum cell counts per cell type within samples,
#' and generating pseudobulk representations.
#'
#' @param pb_sc_counts Gene-by-cell count matrix for PB single-cell data.
#' @param pb_sc_meta Metadata dataframe for PB single-cell data.
#' @param ref_sc_counts Gene-by-cell count matrix for reference single-cell data.
#' @param ref_sc_meta Metadata dataframe for reference single-cell data.
#' @param gene_set Either "intersection" or "custom".
#' @param custom_genes Gene vector used when gene_set = "custom".
#' @param cell_id_col Column containing cell identifiers.
#' @param subtype_col Column for subtype labels.
#' @param major_type_col Column for major cell-type labels.
#' @param sample_id_col Column for sample identifiers.
#' @param min_cells_per_ct Minimum PB cells required for a cell type.
#' @param min_samples_per_ct Minimum PB samples supporting a cell type.
#' @param min_cells_per_ct_per_sample Minimum cells per cell type per sample.
#' @param min_cells_per_gene Minimum PB cells expressing a gene.
#' @param min_gene_total Minimum PB counts required for a gene.
#' @param verbose Print progress messages.
#'
#' @return A list containing reference data, pseudobulk data, gene universe,
#' and a summary report.
#' @importFrom Matrix sparse.model.matrix rowSums
#' @importFrom methods as
#' @importFrom stats complete.cases setNames
#' @export
#' @examples
#' data(cdseq2_toy)
#' x <- cdseq2_toy
#' counts <- x$data$reference
#' ct <- x$data$annotation
#' donors <- x$data$reference_donors
#' meta <- data.frame(cell_id = colnames(counts), celltype = ct,
#'     sampleID = donors)
#' # A shared input here illustrates formatting only, not independent validation.
#' b <- generate_benchmark_dataset(counts, meta, counts, meta,
#'     min_cells_per_ct = 1, min_samples_per_ct = 1,
#'     min_cells_per_gene = 1, min_gene_total = 1, verbose = FALSE)
#' names(b)
generate_benchmark_dataset <- function(
    pb_sc_counts,
    pb_sc_meta,
    ref_sc_counts,
    ref_sc_meta,
    gene_set = c("intersection","custom"),
    custom_genes = NULL,
    cell_id_col = "cell_id",
    subtype_col = NULL,
    major_type_col = "celltype",
    sample_id_col = "sampleID",
    min_cells_per_ct = 5,
    min_samples_per_ct = 2,
    min_cells_per_ct_per_sample = 0,
    min_cells_per_gene = 50,
    min_gene_total = 50,
    verbose = TRUE
){

  gene_set <- match.arg(gene_set)

  if(verbose){
    message("====================================")
    message("Starting benchmark dataset creation")
    message("====================================")
  }

  ## --------------------------------------------------
  ## Align cells with count matrices
  ## --------------------------------------------------

  pb_cells <- intersect(pb_sc_meta[[cell_id_col]], colnames(pb_sc_counts))
  ref_cells <- intersect(ref_sc_meta[[cell_id_col]], colnames(ref_sc_counts))

  pb_sc_meta <- pb_sc_meta[pb_sc_meta[[cell_id_col]] %in% pb_cells,,drop=FALSE]
  ref_sc_meta <- ref_sc_meta[ref_sc_meta[[cell_id_col]] %in% ref_cells,,drop=FALSE]

  pb_sc_counts <- pb_sc_counts[,pb_sc_meta[[cell_id_col]],drop=FALSE]
  ref_sc_counts <- ref_sc_counts[,ref_sc_meta[[cell_id_col]],drop=FALSE]

  if(verbose){
    message("PB cells retained: ", ncol(pb_sc_counts))
    message("Reference cells retained: ", ncol(ref_sc_counts))
  }

  ## --------------------------------------------------
  ## Gene selection
  ## --------------------------------------------------

  if(gene_set=="intersection"){
    genes_use <- intersect(rownames(pb_sc_counts),rownames(ref_sc_counts))
  }

  if(gene_set=="custom"){
    if(is.null(custom_genes)) stop("custom_genes required when gene_set='custom'")
    genes_use <- intersect(custom_genes,rownames(pb_sc_counts))
    genes_use <- intersect(genes_use,rownames(ref_sc_counts))
  }

  if(verbose){
    message("Genes before filtering: ", length(genes_use))
  }

  pb_sc_counts <- pb_sc_counts[genes_use,,drop=FALSE]
  ref_sc_counts <- ref_sc_counts[genes_use,,drop=FALSE]

  ## --------------------------------------------------
  ## Gene filtering
  ## --------------------------------------------------

  gene_filter <- Matrix::rowSums(pb_sc_counts>0) >= min_cells_per_gene &
    Matrix::rowSums(pb_sc_counts) >= min_gene_total

  genes_keep <- rownames(pb_sc_counts)[gene_filter]

  if(length(genes_keep)==0) stop("No genes passed filters")

  pb_sc_counts <- pb_sc_counts[genes_keep,,drop=FALSE]
  ref_sc_counts <- ref_sc_counts[genes_keep,,drop=FALSE]

  if(verbose){
    message("Genes after filtering: ", length(genes_keep))
  }

  ## --------------------------------------------------
  ## Major cell-type filtering
  ## --------------------------------------------------

  ct_counts <- table(pb_sc_meta[[major_type_col]])
  keep_ct <- names(ct_counts[ct_counts >= min_cells_per_ct])

  pb_sc_meta <- pb_sc_meta[pb_sc_meta[[major_type_col]] %in% keep_ct,,drop=FALSE]
  pb_sc_counts <- pb_sc_counts[,pb_sc_meta[[cell_id_col]],drop=FALSE]

  if(verbose){
    message("Cell types retained: ", length(unique(pb_sc_meta[[major_type_col]])))
  }

  ## --------------------------------------------------
  ## Minimum samples per cell type
  ## --------------------------------------------------

  ct_sample_table <- table(
    pb_sc_meta[[major_type_col]],
    pb_sc_meta[[sample_id_col]]
  )

  valid_ct <- rownames(ct_sample_table)[
    rowSums(ct_sample_table>0) >= min_samples_per_ct
  ]

  pb_sc_meta <- pb_sc_meta[pb_sc_meta[[major_type_col]] %in% valid_ct,,drop=FALSE]
  pb_sc_counts <- pb_sc_counts[,pb_sc_meta[[cell_id_col]],drop=FALSE]

  if(verbose){
    message("Cell types after sample-support filter: ",
            length(unique(pb_sc_meta[[major_type_col]])))
  }

  ## --------------------------------------------------
  ## Enforce per-sample minimum cell counts
  ## --------------------------------------------------

  if(min_cells_per_ct_per_sample>0){

    if(verbose){
      message("Applying min_cells_per_ct_per_sample constraint: ",
              min_cells_per_ct_per_sample)
    }

    ct_sample_table <- table(
      pb_sc_meta[[sample_id_col]],
      pb_sc_meta[[major_type_col]]
    )

    keep_samples <- rownames(ct_sample_table)[
      apply(ct_sample_table >= min_cells_per_ct_per_sample,1,all)
    ]

    pb_sc_meta <- pb_sc_meta[
      pb_sc_meta[[sample_id_col]] %in% keep_samples,
      ,drop=FALSE]

    pb_sc_counts <- pb_sc_counts[,pb_sc_meta[[cell_id_col]],drop=FALSE]

  }

  samples <- unique(pb_sc_meta[[sample_id_col]])

  if(verbose){
    message("Samples retained: ", length(samples))
  }

  ## --------------------------------------------------
  ## Pseudobulk sample counts
  ## --------------------------------------------------

  pb_counts <- Matrix::Matrix(
    0,
    nrow=nrow(pb_sc_counts),
    ncol=length(samples),
    sparse=TRUE
  )

  rownames(pb_counts) <- rownames(pb_sc_counts)
  colnames(pb_counts) <- samples

  for(s in samples){

    cells <- pb_sc_meta[
      pb_sc_meta[[sample_id_col]]==s,
      cell_id_col
    ]

    pb_counts[,s] <- Matrix::rowSums(pb_sc_counts[,cells,drop=FALSE])

  }

  ## --------------------------------------------------
  ## Pseudobulk major cell-type matrices
  ## --------------------------------------------------

  pb_major_list <- list()

  for(s in samples){

    cells <- pb_sc_meta[
      pb_sc_meta[[sample_id_col]]==s,
      ,drop=FALSE]

    cts <- unique(cells[[major_type_col]])

    mat <- Matrix::Matrix(
      0,
      nrow=nrow(pb_sc_counts),
      ncol=length(cts),
      sparse=TRUE
    )

    rownames(mat) <- rownames(pb_sc_counts)
    colnames(mat) <- cts

    for(ct in cts){

      ct_cells <- cells[cells[[major_type_col]]==ct,cell_id_col]

      mat[,ct] <- Matrix::rowSums(pb_sc_counts[,ct_cells,drop=FALSE])

    }

    pb_major_list[[s]] <- mat

  }

  ## --------------------------------------------------
  ## Subtype pseudobulk matrices
  ## --------------------------------------------------

  pb_sub_list <- NULL

  if(!is.null(subtype_col)){

    pb_sub_list <- list()

    for(s in samples){

      cells <- pb_sc_meta[
        pb_sc_meta[[sample_id_col]]==s,
        ,drop=FALSE]

      cells <- cells[!is.na(cells[[subtype_col]]),,drop=FALSE]

      subs <- unique(cells[[subtype_col]])

      mat <- Matrix::Matrix(
        0,
        nrow=nrow(pb_sc_counts),
        ncol=length(subs),
        sparse=TRUE
      )

      rownames(mat) <- rownames(pb_sc_counts)
      colnames(mat) <- subs

      for(st in subs){

        st_cells <- cells[cells[[subtype_col]]==st,cell_id_col]

        mat[,st] <- Matrix::rowSums(pb_sc_counts[,st_cells,drop=FALSE])

      }

      pb_sub_list[[s]] <- mat

    }

  }

  ## --------------------------------------------------
  ## Reference profiles
  ## --------------------------------------------------

  ref_profiles <- sapply(unique(ref_sc_meta[[major_type_col]]),function(ct){

    cells <- ref_sc_meta[
      ref_sc_meta[[major_type_col]]==ct,
      cell_id_col
    ]

    Matrix::rowSums(ref_sc_counts[,cells,drop=FALSE]) / length(cells)

  })

  if(verbose){
    message("Reference cell types: ", ncol(ref_profiles))
  }

  ## --------------------------------------------------
  ## Sample metadata
  ## --------------------------------------------------

  sample_meta <- unique(pb_sc_meta[,sample_id_col,drop=FALSE])
  rownames(sample_meta) <- sample_meta[[sample_id_col]]

  ## --------------------------------------------------
  ## Final report
  ## --------------------------------------------------

  report <- list(
    n_genes = nrow(pb_counts),
    n_pb_cells = ncol(pb_sc_counts),
    n_ref_cells = ncol(ref_sc_counts),
    n_samples = length(samples)
  )

  if(verbose){

    message("====================================")
    message("Final dataset summary")
    message("------------------------------------")
    message("Genes: ", report$n_genes)
    message("PB cells: ", report$n_pb_cells)
    message("Reference cells: ", report$n_ref_cells)
    message("Samples: ", report$n_samples)
    message("====================================")

  }

  ## --------------------------------------------------
  ## Output objects
  ## --------------------------------------------------

  pb_sc_data <- list(
    counts = pb_sc_counts,
    cell_type_major = pb_sc_meta[[major_type_col]],
    cell_type_sub = if(!is.null(subtype_col)) pb_sc_meta[[subtype_col]] else NULL,
    sample_id = pb_sc_meta[[sample_id_col]]
  )

  reference <- list(
    counts = ref_sc_counts,
    cell_type_major = ref_sc_meta[[major_type_col]],
    cell_type_sub = if(!is.null(subtype_col)) ref_sc_meta[[subtype_col]] else NULL,
    sample_id = ref_sc_meta[[sample_id_col]]
  )

  pseudobulk <- list(
    pb_counts = pb_counts,
    pb_major_list = pb_major_list,
    pb_sub_list = pb_sub_list,
    sample_meta = sample_meta,
    pb_sc_data = pb_sc_data
  )

  return(list(
    reference = reference,
    pseudobulk = pseudobulk,
    gene_universe = rownames(pb_counts),
    report = report
  ))

}




