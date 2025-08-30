# DMD pipeline For Identifying Dynamic Biomarkers for Tolerance Vs Syngeneic 
## Adapted from Yeung et al., Nature Communications(2020) 
### Builds per-group 3D arrays (genes x time x replicates)
### Fold-change vs. Control (Tol vs Syn) day-mean baseline (with pseudo-count)
### Standardize per gene across (time x reps)
### Exact DMD with rank reduction (drop NA transitions)
### Finite-horizon observability Gramian in reduced space
### Gene weights & optional minimal-redundancy panel; R^2 check

suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
  library(matrixStats)
  library(MASS)   # for ginv
})

set.seed(42)
options(stringsAsFactors = FALSE)

#  Load Dataset ----
# Cohorts/groups/time grid from your study
groups_keep   <- c("Control Accepted", "Tolerance")     # Control (Syn) vs Tol
days_keep     <- c(7, 14, 28, 42, 56, 70)               # canonical days

# Metadata columns in your DESeq2 object
batch_col     <- "Batch"
group_col     <- "Group"
day_col       <- "Day"
animal_col    <- "Animal"
cohort_col    <- "Cohort"

# Pre-filtering (optional)
top_var_genes <- 5000   # set to NULL to keep all genes

# DMD & panel
rank_target   <- 20     # requested reduced rank; auto-clipped to allowed
horizon_T     <- 15     # horizon for Gramian sum
panel_max_genes <- 20   # target panel size (greedy, de-correlated)
corr_thresh   <- 0.50   # max |corr| allowed between chosen genes

# Fold-change pseudo-count (analogous to tutorial pseudoTPM)
pseudo_eps    <- 20


dds_path <- "/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Robjects/dds_IsletTransplant_master.rds"
dds_all  <- readRDS(dds_path)

sel <- colData(dds_all)[[group_col]] %in% groups_keep &
  colData(dds_all)[[day_col]]   %in% days_keep
dds   <- dds_all[, sel, drop = FALSE]

colData(dds)[[group_col]] <- factor(colData(dds)[[group_col]], levels = groups_keep)
colData(dds)[[day_col]]   <- as.integer(colData(dds)[[day_col]])
stopifnot(batch_col %in% colnames(colData(dds)))
colData(dds)[[batch_col]] <- factor(colData(dds)[[batch_col]])

design(dds) <- as.formula(paste0("~ ", group_col, " + factor(", day_col, ")"))
vsd         <- vst(dds, blind = FALSE)
expr        <- assay(vsd)  # genes x samples

#  Remove batch (For visualization/aggregation ONLY)
design_mat  <- model.matrix(~ colData(vsd)[[group_col]] + factor(colData(vsd)[[day_col]]))
expr_adj    <- removeBatchEffect(expr, batch = colData(vsd)[[batch_col]], design = design_mat)

# Optional HVG filter to speed up later steps-TBD
if (!is.null(top_var_genes) && top_var_genes < nrow(expr_adj)) {
  vars <- rowVars(expr_adj)
  keep <- order(vars, decreasing = TRUE)[seq_len(top_var_genes)]
  expr_adj <- expr_adj[keep, , drop = FALSE]
  cat("\nDONE.\n")
}
dim(expr_adj)

# Build 3D Arrays-arr_syn, arr_tol  (genes x time x replicates) =========================
meta <- as.data.frame(colData(vsd))
meta[[day_col]] <- as.integer(meta[[day_col]])

# Replicate = Cohort.Animal (robust for your data)
meta$RepID <- droplevels(interaction(meta[[cohort_col]], meta[[animal_col]],
                                     drop = TRUE, lex.order = TRUE))
rep_col <- "RepID"

to_array_group <- function(M, meta, group_val, days_keep,
                           rep_col = "RepID", group_col = "Group", day_col = "Day") {
  sel_g  <- meta[[group_col]] == group_val
  M_g    <- M[, sel_g, drop = FALSE]
  meta_g <- droplevels(meta[sel_g, ])
  
  reps <- levels(meta_g[[rep_col]])[levels(meta_g[[rep_col]]) %in% meta_g[[rep_col]]]
  m    <- length(days_keep)
  arr  <- array(NA_real_, dim = c(nrow(M_g), m, length(reps)),
                dimnames = list(rownames(M_g), paste0("Day", days_keep), as.character(reps)))
  for (j in seq_along(reps)) {
    rid   <- reps[j]
    sel_r <- meta_g[[rep_col]] == rid
    M_r   <- M_g[, sel_r, drop = FALSE]
    meta_r<- droplevels(meta_g[sel_r, ])
    for (i in seq_along(days_keep)) {
      d <- days_keep[i]
      sel_rd <- meta_r[[day_col]] == d
      if (any(sel_rd)) {
        vec <- if (sum(sel_rd) > 1) rowMeans(M_r[, sel_rd, drop = FALSE]) else
          as.numeric(M_r[, which(sel_rd)[1], drop = FALSE])
        arr[, i, j] <- vec
      }
    }
  }
  arr
}

arr_syn <- to_array_group(expr_adj, meta, groups_keep[1], days_keep)  # Control Accepted
arr_tol <- to_array_group(expr_adj, meta, groups_keep[2], days_keep)  # Tolerance

# Filter Tol replicates: keep ≥ 4 observed days and ≥1 adjacent observed pair (helps transitions)
filter_reps <- function(arr, days_keep,
                        min_obs_days = length(days_keep) - 2,
                        require_adjacent_pair = TRUE) {
  m <- length(days_keep); r <- dim(arr)[3]
  if (r == 0) stop("No replicates available.")
  keep <- logical(r)
  for (j in seq_len(r)) {
    obs <- colSums(!is.na(arr[ , , j, drop = FALSE])) > 0
    has_pair <- any(obs[-m] & obs[-1])
    keep[j]  <- (sum(obs) >= min_obs_days) && (!require_adjacent_pair || has_pair)
  }
  if (!any(keep)) stop("No replicates passed filtering.")
  arr[ , , keep, drop = FALSE]
}
arr_tol <- filter_reps(arr_tol, days_keep)

# ======= FOLD-CHANGE vs Syn day-mean baseline + STANDARDIZE ==================
# Syn baseline: day-wise mean across Syn replicates
syn_day_mean <- apply(arr_syn, c(1,2),
                      function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE))
# Tile baseline to Tol replicate count
n_tol <- dim(arr_tol)[3]; m <- length(days_keep)
arr_syn_tiled <- array(rep(syn_day_mean, times = n_tol),
                       dim = c(nrow(arr_tol), m, n_tol),
                       dimnames = list(rownames(arr_tol), paste0("Day", days_keep), dimnames(arr_tol)[[3]]))

# Element-wise fold-change with pseudo-count
fc_arr <- (arr_tol + pseudo_eps) / (arr_syn_tiled + pseudo_eps)   # genes x m x n_tol

# Standardize per gene across (time x reps) -> z-score
standardize_time_series <- function(arr) {
  G <- dim(arr)[1]; m <- dim(arr)[2]; R <- dim(arr)[3]
  out <- arr
  for (g in seq_len(G)) {
    v  <- as.numeric(arr[g, , , drop = FALSE])
    mu <- mean(v, na.rm = TRUE)
    sdg<- sd(v,   na.rm = TRUE); if (!is.finite(sdg) || sdg == 0) sdg <- 1
    out[g, , ] <- (arr[g, , ] - mu) / sdg
  }
  out
}
fc_arr_norm <- standardize_time_series(fc_arr)  # genes x m x n_tol

# ======= EXACT DMD with rank reduction ======================================
dmd_array <- function(arr, r = NULL, energy = 0.95, trim = FALSE, trimThresh = 1.5e-3,
                      verbose = TRUE) {
  stopifnot(length(dim(arr)) == 3)
  G <- dim(arr)[1]; m <- dim(arr)[2]; R <- dim(arr)[3]
  
  # Build Xp (past) and Xf (future) as genes x ((m-1)*R)
  X0_list <- X1_list <- vector("list", R)
  for (j in seq_len(R)) {
    X0j <- arr[, 1:(m-1), j, drop = FALSE]; dim(X0j) <- c(G, m-1)
    X1j <- arr[, 2: m   , j, drop = FALSE]; dim(X1j) <- c(G, m-1)
    X0_list[[j]] <- X0j
    X1_list[[j]] <- X1j
  }
  Xp <- do.call(cbind, X0_list)        # "Xp"
  Xf <- do.call(cbind, X1_list)        # "Xf"
  rownames(Xp) <- rownames(arr); rownames(Xf) <- rownames(arr)
  
  # Drop transitions with any NA
  ok <- colSums(is.na(Xp)) == 0 & colSums(is.na(Xf)) == 0
  Xp <- Xp[, ok, drop = FALSE]; Xf <- Xf[, ok, drop = FALSE]
  if (ncol(Xp) < 2) stop("Not enough valid transitions after NA filtering.")
  
  # SVD and rank choice
  sv   <- svd(Xp)
  d    <- sv$d
  rmax <- min(length(d), nrow(Xp), ncol(Xp))
  if (is.null(r)) {
    cumE <- cumsum(d^2) / sum(d^2)
    r <- which(cumE >= energy)[1]; if (is.na(r)) r <- rmax
  }
  r <- max(1, min(rank_target, rmax))  # clip to requested rank_target and rmax
  
  U_r <- sv$u[, 1:r, drop = FALSE]
  s_r <- d[1:r]
  V_r <- sv$v[, 1:r, drop = FALSE]
  rownames(U_r) <- rownames(Xp)
  
  # Reduced operator Atilde = U^T Xf V S^{-1}
  Atilde <- t(U_r) %*% Xf %*% V_r %*% diag(1/s_r, nrow = r, ncol = r)
  
  # Optional sparse trimming (mirrors tutorial "trim" idea)
  if (trim) {
    Atilde[abs(Atilde) < trimThresh] <- 0
  }
  
  if (verbose) {
    message(sprintf("DMD: G=%d, valid pairs=%d, r_used=%d (rmax=%d)",
                    G, ncol(Xp), r, rmax))
  }
  list(U_r = U_r, Atilde = Atilde, r = r, Xp = Xp, Xf = Xf)
}

fit <- dmd_array(fc_arr_norm, r = rank_target, verbose = TRUE)
U_r <- fit$U_r
A_r <- fit$Atilde
r_used <- fit$r

# ======= OBSERVABILITY GRAMIAN (reduced) & GENE WEIGHTS ======================
get_Z0_mat <- function(arr) {
  G <- dim(arr)[1]; R <- dim(arr)[3]
  Z0 <- arr[, 1, , drop = FALSE]; dim(Z0) <- c(G, R)   # genes x reps
  rownames(Z0) <- rownames(arr)
  Z0
}

Z0 <- get_Z0_mat(fc_arr_norm)   # genes x reps (first timepoint across reps)

# Align U_r and Z0 by gene names (bulletproof)
common <- intersect(rownames(U_r), rownames(Z0))
Uuse   <- U_r[common, , drop = FALSE]
Z0use  <- Z0[common,  , drop = FALSE]

Z0_red <- t(Uuse) %*% Z0use     # r x reps

gramian_reduced <- function(Atilde, Z0_red, nT = 15L) {
  r <- nrow(Atilde)
  Gred <- matrix(0, r, r)
  At   <- diag(r)
  for (t in 0:(nT-1)) {
    Xt   <- At %*% Z0_red
    Gred <- Gred + Xt %*% t(Xt)
    At   <- Atilde %*% At
  }
  Gred
}
Gred <- gramian_reduced(A_r, Z0_red, nT = horizon_T)

eG <- eigen(Gred, symmetric = FALSE)
v1 <- eG$vectors[, 1, drop = FALSE]

w_genes <- Uuse %*% v1   # genes x 1 (weights in gene space)
gene_weights <- data.frame(
  gene = rownames(Uuse),
  weight = as.numeric(w_genes),
  abs_weight = abs(as.numeric(w_genes)),
  row.names = NULL
)
gene_weights <- gene_weights[order(gene_weights$abs_weight, decreasing = TRUE), ]

cat("\nTop genes by observability weight:\n")
print(head(gene_weights, 15))

# ======= OPTIONAL: GREEDY, LOW-REDUNDANCY PANEL =============================
# build a flattened matrix of trajectories for correlation screening
flatten_mat <- function(arr) {
  G <- dim(arr)[1]; m <- dim(arr)[2]; R <- dim(arr)[3]
  M <- matrix(NA_real_, nrow = G, ncol = m * R)
  for (g in seq_len(G)) M[g, ] <- as.numeric(arr[g, , , drop = FALSE])
  rownames(M) <- rownames(arr)
  M
}
M_flat <- flatten_mat(fc_arr_norm)  # genes x (m*reps)

greedy_panel <- function(ranked_genes, M_flat, k_max = 20, corr_thr = 0.5) {
  picked <- character(0)
  for (g in ranked_genes) {
    if (!(g %in% rownames(M_flat))) next
    v <- M_flat[g, ]
    if (length(picked) == 0) {
      picked <- c(picked, g)
    } else {
      C <- abs(cor(t(M_flat[picked, , drop = FALSE]), v, use = "pairwise.complete.obs"))
      # C is length(picked) x 1; take max
      if (max(C, na.rm = TRUE) <= corr_thr) picked <- c(picked, g)
    }
    if (length(picked) >= k_max) break
  }
  picked
}

panel_genes <- greedy_panel(gene_weights$gene, M_flat,
                            k_max = panel_max_genes, corr_thr = corr_thresh)
cat("\nProposed minimally-redundant panel (≤ ", panel_max_genes, "):\n", sep = "")
print(panel_genes)

# ======= OPTIONAL: R^2 reconstruction check with chosen sensors ==============
to_indices <- function(genes_or_idx, U_r) {
  if (is.character(genes_or_idx)) {
    idx <- match(genes_or_idx, rownames(U_r))
  } else {
    idx <- as.integer(genes_or_idx)
  }
  idx <- unique(idx[!is.na(idx) & idx >= 1 & idx <= nrow(U_r)])
  idx
}

n_step_pred_r2 <- function(Atilde, U_r, Z0_red, sensors, nT = 6L) {
  G <- nrow(U_r); r <- ncol(U_r)
  if (r == 0L || G == 0L) stop("U_r has zero dimension.")
  idx <- to_indices(sensors, U_r)
  p <- length(idx)
  if (p == 0L) stop("No valid sensors after aligning to U_r.")
  C  <- matrix(0, p, G); C[cbind(seq_len(p), idx)] <- 1
  C_r <- C %*% U_r
  
  O_T <- matrix(0, p * nT, r)
  At  <- diag(r)
  for (t in 0:(nT - 1L)) {
    O_T[(p*t + 1L):(p*(t + 1L)), ] <- C_r %*% At
    At <- Atilde %*% At
  }
  
  z0_true <- rowMeans(Z0_red, na.rm = TRUE)
  y       <- O_T %*% matrix(z0_true, ncol = 1)
  z0_hat  <- MASS::ginv(O_T) %*% y
  
  x0_true <- U_r %*% z0_true
  x0_hat  <- U_r %*% z0_hat
  
  ss_res <- sum((x0_true - x0_hat)^2)
  ss_tot <- sum((x0_true - mean(x0_true))^2)
  1 - ss_res/ss_tot
}

# Try with top 10 by weight
top10_names <- head(gene_weights$gene, 10)
R2_top10    <- n_step_pred_r2(A_r, Uuse, Z0_red, sensors = top10_names, nT = length(days_keep))
cat(sprintf("\nReconstruction R^2 with top 10 genes: %.3f\n", R2_top10))

# And with your de-correlated panel
if (length(panel_genes) >= 5) {
  R2_panel <- n_step_pred_r2(A_r, Uuse, Z0_red, sensors = panel_genes, nT = length(days_keep))
  cat(sprintf("Reconstruction R^2 with panel (%d genes): %.3f\n", length(panel_genes), R2_panel))
}

# ======= PACKAGE OUTPUT ======================================================
dmd_results <- list(
  config = list(groups_keep = groups_keep, days_keep = days_keep,
                rank_used = r_used, horizon_T = horizon_T,
                pseudo_eps = pseudo_eps,
                panel_max_genes = panel_max_genes, corr_thresh = corr_thresh),
  arrays = list(arr_syn = arr_syn, arr_tol = arr_tol, fc_arr_norm = fc_arr_norm),
  dmd    = list(U_r = U_r, Atilde = A_r, Z0_red = Z0_red),
  weights = gene_weights,
  panel   = panel_genes
)

cat("\nDONE.\n")
