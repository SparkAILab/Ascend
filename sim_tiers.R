# ======================================================================
# Two-tier simulator with controllable violations of the tier assumption
# ----------------------------------------------------------------------
# Reviewer comment 4: how does ASCEND degrade when the background /
# foreground split is wrong or incomplete?
#
# One joint linear-Gaussian DAG over all true variables is simulated; the
# analyst then sees the variables through a (possibly wrong) tier labelling.
# Violations, each a fraction in [0, 1):
#   p_x_as_z  true foreground variables labelled as background. These are
#             descendants of other X, so Assumption 1 is violated.
#   p_z_as_x  true background variables labelled as foreground (fewer
#             anchors; Assumption 1 still holds).
#   p_hide_z  background variables not measured (incomplete tier: latent
#             upstream confounders of X).
#   p_feedback background variables that are *caused by* a foreground
#             variable (X -> Z reverse causation / feedback), placed after
#             that X in the causal order. Assumption 1 is violated.
#
# Ground truth is always the ancestral relation among the variables the
# analyst labels foreground, computed on the full joint DAG (so paths
# through hidden or mislabelled variables count).
#
# With all violation rates at 0 this is the same data-generating process as
# sim_dat(z_scale = TRUE): unit-variance background, Z -> X coefficients
# +/-1, X -> X coefficients x_effect / sd(parent), noise set by r2.
# ======================================================================

sim_tiers <- function(n, d_z, d_x,
                      r2 = 0.5, sp = 0.7, p_cross = 0.10, x_effect = 0.9,
                      p_x_as_z = 0, p_z_as_x = 0, p_hide_z = 0, p_feedback = 0,
                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p_edge <- min(max(1 - sp, 0.01), 0.4)

  # --- true graph: order Z (random DAG), then X (random DAG) ----------------
  zn <- paste0("Z", seq_len(d_z)); xn <- paste0("X", seq_len(d_x))
  V  <- c(zn, xn); p <- length(V)
  A  <- matrix(0, p, p, dimnames = list(V, V))              # A[parent, child]
  up_dag <- function(k) { B <- matrix(0, k, k)
  if (k > 1) for (i in 1:(k - 1)) for (j in (i + 1):k) if (runif(1) < p_edge) B[i, j] <- 1
  B }
  A[zn, zn] <- up_dag(d_z)                                  # index order = causal order
  A[xn, xn] <- up_dag(d_x)
  A[zn, xn] <- matrix(rbinom(d_z * d_x, 1, p_cross), d_z, d_x)

  # causal order of all nodes; feedback Z nodes are moved after an X
  order <- c(zn, xn)
  n_fb  <- round(p_feedback * d_z)
  fb    <- if (n_fb > 0) sample(zn, n_fb) else character(0)
  fb_parent <- setNames(character(0), character(0))
  for (z in fb) {
    pos <- sample.int(max(1, d_x - 1), 1)                   # after X_pos, before X_{pos+1}
    xp  <- xn[pos]
    A[z, ] <- 0                                             # drop old out-edges ...
    later  <- xn[seq_len(d_x) > pos]
    A[z, later] <- rbinom(length(later), 1, p_cross * 3)    # ... re-draw them to later X only
    A[xp, z] <- 1                                           # the reverse-causal edge X -> Z
    fb_parent[z] <- xp
    order <- append(setdiff(order, z), z, after = match(xp, setdiff(order, z)))
  }

  # --- sample in causal order -----------------------------------------------
  D <- matrix(0, n, p, dimnames = list(NULL, V))
  is_z <- setNames(V %in% zn, V)
  noise_sd <- function(sig) { v <- var(sig); if (!is.finite(v) || v == 0) 1 else sqrt(v * (1 - r2) / r2) }
  coef <- matrix(0, p, p, dimnames = list(V, V))
  for (v in order) {
    par <- names(which(A[, v] == 1))
    if (is_z[v]) {
      b <- if (length(par)) rnorm(length(par)) else numeric(0)
      if (v %in% fb) b[par == fb_parent[v]] <- x_effect * sample(c(-1, 1), 1)
      sig <- if (length(par)) as.numeric(scale(D[, par, drop = FALSE], scale = TRUE) %*% b) else 0
      D[, v] <- sig + rnorm(n)
      D[, v] <- (D[, v] - mean(D[, v])) / sd(D[, v])        # unit-variance background
    } else {
      pz <- par[is_z[par]]; px <- par[!is_z[par]]
      sig <- rep(0, n)
      if (length(pz)) { b <- sample(c(-1, 1), length(pz), TRUE); sig <- sig + D[, pz, drop = FALSE] %*% b }
      if (length(px)) {
        s <- apply(D[, px, drop = FALSE], 2, sd); s[s == 0] <- 1
        b <- x_effect / s * sample(c(-1, 1), length(px), TRUE)
        sig <- sig + D[, px, drop = FALSE] %*% b
      }
      sig <- as.numeric(sig); if (var(sig) < 1e-6) sig <- rnorm(n)
      D[, v] <- sig + rnorm(n, sd = noise_sd(sig))
    }
  }

  # --- analyst's view --------------------------------------------------------
  n_hide <- round(p_hide_z * d_z)
  hidden <- if (n_hide > 0) sample(setdiff(zn, fb), min(n_hide, length(setdiff(zn, fb)))) else character(0)
  n_x2z  <- round(p_x_as_z * d_x)
  x_as_z <- if (n_x2z > 0) sample(xn, n_x2z) else character(0)
  n_z2x  <- round(p_z_as_x * d_z)
  z_as_x <- if (n_z2x > 0) sample(setdiff(zn, c(hidden, fb)), min(n_z2x, length(setdiff(zn, c(hidden, fb))))) else character(0)

  obs_z <- c(setdiff(zn, c(hidden, z_as_x)), x_as_z)
  obs_x <- c(setdiff(xn, x_as_z), z_as_x)
  dat <- as.data.frame(D[, c(obs_z, obs_x), drop = FALSE])
  lab <- c(paste0("z", seq_along(obs_z)), paste0("x", seq_along(obs_x)))
  true_name <- setNames(c(obs_z, obs_x), lab)
  colnames(dat) <- lab

  # --- truth: ancestral closure on the full joint DAG, restricted to obs X ---
  R <- (A > 0) * 1; P <- R
  for (k in seq_len(p - 1)) { P <- ((P %*% A) > 0) * 1; if (!any(P > 0)) break; R <- ((R + P) > 0) * 1 }
  xl <- paste0("x", seq_along(obs_x))
  truth <- R[obs_x, obs_x, drop = FALSE]; dimnames(truth) <- list(xl, xl); diag(truth) <- 0

  # adj_xx in sim_dat's convention (child, parent) for the observed foreground,
  # direct edges only (paths through other variables are in `truth`)
  adj_xx <- t(A[obs_x, obs_x, drop = FALSE]); dimnames(adj_xx) <- list(xl, xl); diag(adj_xx) <- NA

  list(dat = dat, adj_xx = adj_xx, truth = truth, A_full = A, order = order,
       labels = true_name,
       violations = list(x_as_z = x_as_z, z_as_x = z_as_x, hidden_z = hidden,
                         feedback_z = fb),
       params = list(n = n, d_z = length(obs_z), d_x = length(obs_x), r2 = r2, sp = sp,
                     lin_pr = 1, p_cross = p_cross, x_effect = x_effect,
                     p_x_as_z = p_x_as_z, p_z_as_x = p_z_as_x,
                     p_hide_z = p_hide_z, p_feedback = p_feedback))
}
