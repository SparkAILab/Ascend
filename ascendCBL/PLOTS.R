# ======================================================================
# make_figure.R
#
# Build the headline ASCEND vs CBL figure from results.csv.
#
# Usage:
#   Rscript make_figure.R                       # reads ./results.csv
#   Rscript make_figure.R path/to/results.csv   # explicit path
#
# Outputs (next to the CSV by default):
#   ascend_vs_cbl.pdf
#   ascend_vs_cbl.png
#
# Packages: ggplot2, patchwork, data.table, scales, grid
# ======================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

args     <- commandArgs(trailingOnly = TRUE)
CSV_PATH <- if (length(args) >= 1) args[1] else "results.csv"
OUT_DIR  <- if (length(args) >= 2) args[2] else dirname(normalizePath(CSV_PATH))
stopifnot(file.exists(CSV_PATH))

# ---- Palette / theme ---------------------------------------------------
C_ASCEND  <- "#0072B2"   # blue
C_CBL     <- "#E69F00"   # orange
C_TIMEOUT <- "#C44E52"   # red
GREY      <- "#5A5A5A"

method_levels <- c("CBL (baseline)", "ASCEND (ours)")
method_pal    <- c("CBL (baseline)" = C_CBL, "ASCEND (ours)" = C_ASCEND)

theme_paper <- function(base = 10) {
  theme_bw(base_size = base) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey85"),
      panel.border       = element_blank(),
      axis.line          = element_line(linewidth = 0.45, colour = "grey20"),
      axis.ticks         = element_line(linewidth = 0.45, colour = "grey20"),
      plot.title         = element_text(face = "bold", size = base + 1.5,
                                        margin = margin(b = 4)),
      plot.title.position = "plot",
      strip.background   = element_blank(),
      legend.position    = "none"   # one shared legend on the figure
    )
}

# ---- Load & label ------------------------------------------------------
d <- fread(CSV_PATH)
d[, timed_out := status == "timeout"]
d[, method_lbl := factor(
  ifelse(method == "ascend", "ASCEND (ours)", "CBL (baseline)"),
  levels = method_levels)]

ok <- d[status == "ok"]

# Aggregator: mean and SEM over seeds.
agg <- function(df, group_cols, metric) {
  setDT(df)
  df[, .(mean = mean(get(metric), na.rm = TRUE),
         sd   = sd(get(metric),   na.rm = TRUE),
         n_   = sum(!is.na(get(metric)))),
     by = group_cols][
       , sem := sd / sqrt(pmax(n_, 1))][]
}

# Slice helpers — three 1-D sweeps through the default point.
slice_dx <- d[n == 1024 & d_z == 10]
slice_dz <- d[n == 1024 & d_x == 5]

# ---- Generic sweep panel ----------------------------------------------
# `metric` along `axis`. ASCEND/CBL ribbons = ±1 SEM. If `show_timeouts`
# is TRUE we overlay red ✕ at the mean censored time for cells where
# CBL timed out.
sweep_panel <- function(df, axis, metric, ylab, title,
                        log_y = FALSE, show_timeouts = FALSE,
                        x_breaks = NULL, x_label = NULL,
                        annot = NULL) {
  
  okd <- df[status == "ok"]
  a   <- agg(okd, c("method_lbl", axis), metric)
  setnames(a, axis, "x")
  a <- a[!is.na(mean)]
  
  p <- ggplot(a, aes(x = x, y = mean, colour = method_lbl, fill = method_lbl,
                     group = method_lbl)) +
    geom_ribbon(aes(ymin = pmax(mean - sem, if (log_y) 1e-9 else -Inf),
                    ymax = mean + sem),
                alpha = 0.18, colour = NA) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2.4, shape = 21, colour = "white", stroke = 0.7) +
    scale_colour_manual(values = method_pal) +
    scale_fill_manual  (values = method_pal) +
    labs(title = title, x = x_label, y = ylab) +
    theme_paper()
  
  if (!is.null(x_breaks)) p <- p + scale_x_continuous(breaks = x_breaks)
  if (log_y)              p <- p + scale_y_log10(labels = label_comma())
  
  if (show_timeouts) {
    t <- df[timed_out == TRUE & method == "cbl"]
    if (nrow(t) > 0) {
      t_mean <- t[, .(y = mean(time_sec)), by = axis]
      setnames(t_mean, axis, "x")
      p <- p + geom_point(data = t_mean,
                          aes(x = x, y = y),
                          colour = C_TIMEOUT, shape = 4, size = 4,
                          stroke = 1.4, inherit.aes = FALSE)
    }
  }
  
  if (!is.null(annot)) {
    for (a_ in annot) {
      p <- p + annotate("text",
                        x = a_$x, y = a_$y, label = a_$label,
                        hjust = a_$hjust %||% 0, vjust = a_$vjust %||% 0.5,
                        colour = a_$colour, fontface = "bold", size = 3.4)
    }
  }
  p
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- (a) Runtime vs d_x ----------------------------------------------
p_a <- sweep_panel(
  slice_dx, "d_x", "time_sec",
  ylab = "Wall-clock time (s)",
  title = expression(bold("(a) Runtime vs " * d[x])),
  log_y = TRUE, show_timeouts = TRUE,
  x_breaks = c(5, 10, 15, 20),
  x_label = expression(d[x] * "  (foreground variables)")
) +
  annotate("text", x = 5, y = Inf, label = "n=1024,  d_z=10",
           hjust = -0.05, vjust = 1.6, colour = GREY, size = 3.0)

# ---- (b) Runtime vs d_z ----------------------------------------------
# Compute the slope annotations from the data themselves.
slope_label <- function(method_name, axis_var, slice) {
  s <- slice[status == "ok" & method == method_name,
             .(t = mean(time_sec)), by = axis_var]
  setnames(s, axis_var, "x")
  s <- s[!is.na(t) & t > 0]
  if (nrow(s) < 2) return(NA_real_)
  unname(coef(lm(log(t) ~ log(x), data = s))[2])
}
slope_cbl_dz    <- slope_label("cbl",    "d_z", slice_dz)
slope_ascend_dz <- slope_label("ascend", "d_z", slice_dz)

p_b <- sweep_panel(
  slice_dz, "d_z", "time_sec",
  ylab = "Wall-clock time (s)",
  title = expression(bold("(b) Runtime vs " * d[z])),
  log_y = TRUE, show_timeouts = TRUE,
  x_breaks = c(10, 20, 30, 40, 50),
  x_label = expression(d[z] * "  (background variables)")
) +
  annotate("text", x = 10, y = Inf, label = "n=1024,  d_x=5",
           hjust = -0.05, vjust = 1.6, colour = GREY, size = 3.0) +
  annotate("text", x = 35, y = 200,
           label = sprintf("CBL: t \u221d d_z^{%+0.2f}", slope_cbl_dz),
           colour = C_CBL, fontface = "bold", hjust = 0, size = 3.4) +
  annotate("text", x = 20, y = 0.08,
           label = sprintf("ASCEND: t \u221d d_z^{%+0.2f}", slope_ascend_dz),
           colour = C_ASCEND, fontface = "bold", hjust = 0, size = 3.4)

# ---- (c) Speedup factor bar chart ------------------------------------
# Use mean time per (method, axis-value) including timeouts.
speedup_along <- function(slice, axis) {
  a <- slice[, .(t = mean(time_sec, na.rm = TRUE)), by = c("method", axis)]
  w <- dcast(a, as.formula(paste(axis, "~ method")), value.var = "t")
  w[, speedup := cbl / ascend]
  w
}
sp_dx <- speedup_along(slice_dx, "d_x")
sp_dz <- speedup_along(slice_dz, "d_z")

# Which axis values contained a CBL timeout?
to_dx <- unique(slice_dx[timed_out == TRUE]$d_x)
to_dz <- unique(slice_dz[timed_out == TRUE]$d_z)

bars <- rbind(
  sp_dx[, .(label = paste0("d_x=", d_x),
            value = speedup,
            group = "varying d_x",
            timed_out = d_x %in% to_dx,
            order_key = d_x)],
  sp_dz[, .(label = paste0("d_z=", d_z),
            value = speedup,
            group = "varying d_z",
            timed_out = d_z %in% to_dz,
            order_key = 100 + d_z)]   # offset so d_z bars come after d_x bars
)
bars[, label := factor(label, levels = label[order(order_key)])]

bar_pal <- c("varying d_x" = "#7FB3D5", "varying d_z" = "#1f6f9c")

p_c <- ggplot(bars, aes(x = label, y = value, fill = group)) +
  geom_col(width = 0.78,
           colour = ifelse(bars$timed_out, C_TIMEOUT, "white"),
           linewidth = ifelse(bars$timed_out, 1.1, 0.5)) +
  geom_text(aes(label = paste0(format(round(value), big.mark = ","), "\u00d7")),
            vjust = -0.4, size = 3.3, colour = "grey20") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = GREY,
             linewidth = 0.4) +
  scale_y_log10(labels = label_comma(),
                expand = expansion(mult = c(0.05, 0.22))) +
  scale_fill_manual(values = bar_pal) +
  labs(title = "(c) ASCEND speedup factor",
       x = NULL,
       y = "Speedup  (CBL time / ASCEND time)") +
  theme_paper() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5),
        legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.title = element_blank(),
        legend.background = element_rect(fill = alpha("white", 0.7),
                                         colour = NA),
        legend.key.size = unit(0.4, "cm"),
        legend.text  = element_text(size = 8.5))

# Mark the x-axis tick labels of timeout bars in red, and add a small
# explanatory note inside the panel.
if (any(bars$timed_out)) {
  tick_cols <- ifelse(levels(bars$label) %in%
                        as.character(bars$label[bars$timed_out]),
                      C_TIMEOUT, "grey20")
  p_c <- p_c +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5,
                                     colour = tick_cols)) +
    annotate("text", x = 1, y = max(bars$value) * 2.4,
             label = "red border / red tick = CBL hit timeout\n(speedup is a lower bound)",
             hjust = 0, vjust = 1, size = 2.8, colour = C_TIMEOUT,
             fontface = "italic")
}

# ---- (d) Computational work (CI tests) -------------------------------
p_d <- sweep_panel(
  slice_dz, "d_z", "n_ci_tests",
  ylab = "Conditional independence\ntests performed",
  title = "(d) Computational work",
  log_y = TRUE,
  x_breaks = c(10, 20, 30, 40, 50),
  x_label = expression(d[z])
) +
  annotate("text", x = 10, y = Inf, label = "n=1024,  d_x=5",
           hjust = -0.05, vjust = 1.6, colour = GREY, size = 3.0)

# ---- (e) Decisiveness ------------------------------------------------
p_e <- sweep_panel(
  slice_dz, "d_z", "coverage",
  ylab = "Coverage  (1 \u2212 NA fraction)",
  title = "(e) Decisiveness",
  x_breaks = c(10, 20, 30, 40, 50),
  x_label = expression(d[z])
) +
  coord_cartesian(ylim = c(0.45, 1.03)) +
  annotate("text", x = 10, y = 0.47, label = "n=1024,  d_x=5",
           hjust = -0.05, vjust = 0, colour = GREY, size = 3.0)

# ---- (f) Accuracy vs d_x --------------------------------------------
p_f <- sweep_panel(
  slice_dx, "d_x", "f1",
  ylab = "F1 score",
  title = expression(bold("(f) Accuracy vs " * d[x])),
  x_breaks = c(5, 10, 15, 20),
  x_label = expression(d[x])
) +
  coord_cartesian(ylim = c(-0.02, 1.03)) +
  annotate("text", x = 5, y = -0.02, label = "n=1024,  d_z=10",
           hjust = -0.05, vjust = -0.4, colour = GREY, size = 3.0) +
  # Call out the missing CBL point at d_x=5
  annotate("text", x = 5.3, y = 0.05,
           label = "CBL: 0 TPs across\nall 5 seeds (F1 undef.)",
           hjust = 0, vjust = 0, size = 2.8, colour = C_CBL,
           fontface = "italic")

# ---- Assemble & save -------------------------------------------------
# Custom shared legend strip across the top.
top_legend <- ggplot() +
  annotate("segment", x = 0.06, xend = 0.10, y = 0.5, yend = 0.5,
           colour = C_ASCEND, linewidth = 1.2) +
  annotate("point",   x = 0.08,             y = 0.5,
           colour = C_ASCEND, fill = C_ASCEND, shape = 21, size = 3.0,
           stroke = 0.9) +
  annotate("text",    x = 0.115, y = 0.5, hjust = 0,
           label = "ASCEND (ours)", size = 3.6) +
  annotate("segment", x = 0.30, xend = 0.34, y = 0.5, yend = 0.5,
           colour = C_CBL, linewidth = 1.2) +
  annotate("point",   x = 0.32,             y = 0.5,
           colour = C_CBL, fill = C_CBL, shape = 21, size = 3.0, stroke = 0.9) +
  annotate("text",    x = 0.355, y = 0.5, hjust = 0,
           label = "CBL (baseline)", size = 3.6) +
  annotate("point",   x = 0.555, y = 0.5,
           colour = C_TIMEOUT, shape = 4, size = 4, stroke = 1.5) +
  annotate("text",    x = 0.575, y = 0.5, hjust = 0,
           label = "CBL exceeded 1 h timeout", size = 3.6) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void()

fig <- (p_a | p_b | p_c) / (p_d | p_e | p_f) &
  theme(plot.margin = margin(6, 8, 6, 8))

fig <- top_legend / fig +
  plot_layout(heights = c(0.06, 1)) +
  plot_annotation(
    title = "ASCEND vs CBL  \u2014  scaling, computational cost, and discovery quality  (5 seeds, ribbons = \u00b11 SEM)",
    theme = theme(plot.title = element_text(face = "bold", size = 12.5,
                                            hjust = 0.5,
                                            margin = margin(b = 4)))
  )

pdf_path <- file.path(OUT_DIR, "ascend_vs_cblRR.pdf")
png_path <- file.path(OUT_DIR, "ascend_vs_cblRR.png")
ggsave(pdf_path, fig, width = 13, height = 7.4, device = cairo_pdf)
ggsave(png_path, fig, width = 13, height = 7.4, dpi = 300)

cat("Wrote:\n  ", pdf_path, "\n  ", png_path, "\n", sep = "")