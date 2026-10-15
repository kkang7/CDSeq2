# Internal canonical feature definitions. Original API MD5: 51558d04d14bf1578949d37696d3501f.
# Renamed only to avoid collisions with earlier CDSeq2 scoring APIs.
.cdsdeg_f_base <-
c("contribution_mean", "expected_log_mean_count", "ref_specificity", "bulk_logCPM_sd", "theta_mean", "bulk_max_donor_share", 
"composition_loo_r2", "bulk_theta_r")

.cdsdeg_f_cross <-
c("other_max_abs_effect", "other_mean_abs_effect", "self_abs_share", "other_same_sign_fraction")

.cdsdeg_f_assert <-
function (ok, ...) 
{
    if (length(ok) != 1L || is.na(ok) || !ok) 
        stop(..., call. = FALSE)
}

.cdsdeg_f_cols <-
function (d, n, label) 
.cdsdeg_f_assert(is.data.frame(d) && all(n %in% names(d)), label, " requires columns: ", paste(n, collapse = ", "))

.cdsdeg_f_flag <-
function (x, label) 
{
    .cdsdeg_f_assert(is.logical(x) && !anyNA(x), label, " must be logical with no missing values.")
    x
}

.cdsdeg_f_ids <-
function (d, label) 
{
    .cdsdeg_f_cols(d, c("scenario", "gene", "cell_type"), label)
    for (n in c("scenario", "gene", "cell_type")) .cdsdeg_f_assert((is.character(d[[n]]) || is.factor(d[[n]])) && !anyNA(d[[n]]) && 
        all(nzchar(as.character(d[[n]]))) && !any(grepl("\r", d[[n]], fixed = TRUE)), label, " has invalid ", n)
    .cdsdeg_f_assert(!anyDuplicated(.cdsdeg_f_key(d)), label, " has duplicated scenario/gene/cell-type keys.")
}

.cdsdeg_f_key <-
function (d) 
paste(d$scenario, d$gene, d$cell_type, sep = "\r")

.cdsdeg_f_scenario <-
function (d) 
{
    .cdsdeg_f_assert(is.data.frame(d), "Input must be a data.frame.")
    if (!"scenario" %in% names(d)) 
        d$scenario <- rep("dataset", nrow(d))
    d
}

.cdsdeg_f_context <-
function (context) 
{
    .cdsdeg_f_assert(is.list(context) && !is.null(names(context)) && !anyDuplicated(names(context)) && all(nzchar(names(context))), 
        "context must be an explicitly named list.")
    .cdsdeg_f_assert(identical(context$feature_schema, "pilot_logCPM_r"), "Use canonical log2(CPM+0.5) observables; context$feature_schema must be 'pilot_logCPM_r'.")
    .cdsdeg_f_assert(identical(context$effect_scale, "log2CPM_contrast"), "context$effect_scale must be 'log2CPM_contrast'; simple log2 ratios are a different target.")
    .cdsdeg_f_assert(is.character(context$call_definition) && length(context$call_definition) == 1L && !is.na(context$call_definition) && 
        nzchar(context$call_definition), "Declare context$call_definition, e.g. 'saved_global_BH_0.05'.")
    invisible(context)
}

.cdsdeg_f_rowcor <-
function (x, y) 
{
    xc <- x - rowMeans(x)
    yc <- y - mean(y)
    den <- sqrt(rowSums(xc^2) * sum(yc^2))
    z <- drop(xc %*% yc)/den
    z[!is.finite(z) | den < 1e-14] <- NA_real_
    pmin(1, pmax(-1, z))
}

.cdsdeg_observables <-
function (bulk, reference_mean, rna_proportions) 
{
    .cdsdeg_f_assert(is.matrix(bulk) && is.numeric(bulk) && is.matrix(reference_mean) && is.numeric(reference_mean) && 
        is.matrix(rna_proportions) && is.numeric(rna_proportions), "Feature inputs must be numeric matrices.")
    g <- rownames(bulk)
    s <- colnames(bulk)
    ct <- colnames(reference_mean)
    axes <- list(g, s, ct, rownames(reference_mean), rownames(rna_proportions), colnames(rna_proportions))
    .cdsdeg_f_assert(all(vapply(axes, function(z) !is.null(z) && length(z) > 0 && !anyNA(z) && all(nzchar(z)) && !anyDuplicated(z), 
        TRUE)), "Matrices require unique, nonmissing axis names.")
    .cdsdeg_f_assert(setequal(g, rownames(reference_mean)) && setequal(s, rownames(rna_proportions)) && setequal(ct, 
        colnames(rna_proportions)), "Matrix axes disagree.")
    mu <- reference_mean[g, ct, drop = FALSE]
    P <- rna_proportions[s, ct, drop = FALSE]
    B <- bulk
    .cdsdeg_f_assert(all(is.finite(B) & B >= 0) && all(is.finite(mu) & mu >= 0) && all(is.finite(P) & P >= 0), "Counts, reference profiles and estimated fractions must be finite and nonnegative.")
    .cdsdeg_f_assert(length(ct) >= 2L && all(colSums(B) > 0) && all(rowSums(mu) > 0) && max(abs(colSums(mu) - 1)) < 
        1e-08 && max(abs(rowSums(P) - 1)) < 1e-08, "Invalid panel normalization.")
    logB <- log2(sweep(B, 2, colSums(B), "/") * 1e+06 + 0.5)
    n <- ncol(B)
    X <- cbind(1, P[, -ncol(P), drop = FALSE])
    q <- qr(X)
    .cdsdeg_f_assert(q$rank == ncol(X) && n > ncol(X) + 1, "RNA-composition design is singular or too small.")
    h <- rowSums(qr.Q(q)^2)
    .cdsdeg_f_assert(all(h < 1 - 1e-10), "LOO leverage reaches one.")
    residual <- t(qr.resid(q, t(logB)))
    sst <- rowSums((logB - rowMeans(logB))^2)
    loo <- 1 - rowSums(sweep(residual, 2, 1 - h, "/")^2)/sst
    loo[!is.finite(loo)] <- NA_real_
    share <- apply(B, 1, max)/rowSums(B)
    share[!is.finite(share)] <- NA_real_
    den <- mu %*% t(P)
    .cdsdeg_f_assert(all(den > 0), "Reference mixture contains zero denominators.")
    rows <- lapply(ct, function(k) {
        w <- outer(mu[, k], P[, k])/den
        data.frame(gene = g, cell_type = k, cohort_n = n, contribution_mean = rowMeans(w), expected_log_mean_count = log2(rowMeans(w * 
            B) + 0.5), ref_specificity = mu[, k]/rowSums(mu), bulk_logCPM_sd = sqrt(sst/(n - 1)), theta_mean = mean(P[, 
            k]), bulk_max_donor_share = share, composition_loo_r2 = loo, bulk_theta_r = .cdsdeg_f_rowcor(logB, P[, 
            k]), stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    attr(out, "feature_schema") <- "pilot_logCPM_r"
    out
}

.cdsdeg_prepare_canonical <-
function (observables, de_results, context) 
{
    .cdsdeg_f_context(context)
    schema <- attr(observables, "feature_schema")
    .cdsdeg_f_assert(is.null(schema) || identical(schema, "pilot_logCPM_r"), "Observable schema mismatch.")
    f <- .cdsdeg_f_scenario(observables)
    d <- .cdsdeg_f_scenario(de_results)
    .cdsdeg_f_ids(f, "observables")
    .cdsdeg_f_ids(d, "de_results")
    .cdsdeg_f_cols(f, .cdsdeg_f_base, "observables")
    .cdsdeg_f_cols(d, c("estimated_effect", "p_value", "called"), "de_results")
    .cdsdeg_f_assert(all(vapply(f[.cdsdeg_f_base], is.numeric, TRUE)) && is.numeric(d$estimated_effect) && is.numeric(d$p_value), 
        "Predictors and DE effects/p-values must be numeric.")
    .cdsdeg_f_flag(d$called, "called")
    .cdsdeg_f_assert(!any(is.infinite(as.matrix(f[.cdsdeg_f_base]))) && !any(is.infinite(d$estimated_effect)) && !any(is.infinite(d$p_value)), 
        "Use NA for missing inputs; infinite values are invalid.")
    .cdsdeg_f_assert(all(is.na(d$p_value) | (d$p_value >= 0 & d$p_value <= 1)), "p_value must lie in [0,1] or be NA.")
    j <- match(.cdsdeg_f_key(d), .cdsdeg_f_key(f))
    .cdsdeg_f_assert(!anyNA(j), "Missing observable rows for some DE keys.")
    x <- as.matrix(f[j, .cdsdeg_f_base, drop = FALSE])
    for (nm in c("contribution_mean", "ref_specificity", "theta_mean", "bulk_max_donor_share")) .cdsdeg_f_assert(all(is.na(x[, 
        nm]) | (x[, nm] >= 0 & x[, nm] <= 1 + 1e-10)), nm, " must lie in [0,1].")
    .cdsdeg_f_assert(all(is.na(x[, "bulk_theta_r"]) | abs(x[, "bulk_theta_r"]) <= 1 + 1e-10), "Invalid correlation.")
    for (nm in c("contribution_mean", "ref_specificity", "bulk_logCPM_sd", "theta_mean", "bulk_max_donor_share")) {
        .cdsdeg_f_assert(!any(x[, nm] < 0, na.rm = TRUE), nm, " must be nonnegative.")
        x[, nm] <- log2(x[, nm] + 1e-06)
    }
    a <- matrix(NA_real_, nrow(d), 7L, dimnames = list(NULL, c("estimated_abs_effect", "estimated_sign", "estimated_log_p_support", 
        .cdsdeg_f_cross)))
    a[, "estimated_abs_effect"] <- log1p(abs(d$estimated_effect))
    a[, "estimated_sign"] <- sign(d$estimated_effect)
    a[, "estimated_log_p_support"] <- pmin(50, -log10(pmax(d$p_value, 1e-50)))
    counts <- integer(nrow(d))
    for (ix in split(seq_len(nrow(d)), paste(d$scenario, d$gene, sep = "\r"))) for (i in ix) {
        other <- d$estimated_effect[setdiff(ix, i)]
        other <- other[is.finite(other)]
        counts[i] <- length(other)
        if (!length(other)) 
            next
        a[i, "other_max_abs_effect"] <- log1p(max(abs(other)))
        a[i, "other_mean_abs_effect"] <- log1p(mean(abs(other)))
        a[i, "self_abs_share"] <- abs(d$estimated_effect[i])/(abs(d$estimated_effect[i]) + sum(abs(other)) + 1e-08)
        a[i, "other_same_sign_fraction"] <- mean(sign(other) == sign(d$estimated_effect[i]))
    }
    out <- cbind(d[c("scenario", "gene", "cell_type", "called")], as.data.frame(x), as.data.frame(a))
    out$other_effect_count <- counts
    out$cross_type_available <- counts > 0L
    out$feature_schema <- rep("deg_features_v0.4.0", nrow(out))
    rownames(out) <- NULL
    out
}

