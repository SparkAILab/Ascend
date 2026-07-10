# =============================================================================
# Generate Fig: F1 at matched edge count, n=2000, all (R2, sp) cells.
# Reads the single combined raw CSV (now contains both R2 = 0.5 and 0.7)
# and produces bench_f1_heatmap.pdf
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

d <- fread("benchmark_v2_main_raw.csv")
d[, r2_label := as.character(r2)]

# Filter to n=2000 for the headline figure
dd <- d[n == 2000L]

# Per-cell summary
summ <- dd[, .(F1_m = mean(f1, na.rm = TRUE),
               F1_se = sd(f1, na.rm = TRUE) / sqrt(.N)),
           by = .(method, sp, r2_label)]
summ[, method := factor(method, levels = c("ASCEND","GENIE3","ARACNE","WGCNA"))]
summ[, sp_label := factor(paste0("sp = ", sp),
                          levels = c("sp = 0.5","sp = 0.7","sp = 0.9"))]
summ[, r2_facet := factor(paste0("R^2 == ", r2_label))]

cols <- c(ASCEND = "#C55A11", GENIE3 = "#1F4E79",
          ARACNE = "#538135", WGCNA  = "#7030A0")

p <- ggplot(summ, aes(x = method, y = F1_m, fill = method)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = pmax(0, F1_m - F1_se), ymax = F1_m + F1_se),
                width = 0.25, linewidth = 0.4) +
  facet_grid(r2_facet ~ sp_label, labeller = label_parsed) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_y_continuous(limits = c(0, 0.75), breaks = seq(0, 0.7, 0.1),
                     expand = c(0, 0)) +
  labs(x = NULL, y = "F1 at matched edge count",
       title = "F1 across the simulation sweep (n = 2000, 50 reps/cell)") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 0, size = 9))

ggsave("bench_f1_heatmap.pdf", p, width = 9, height = 5.5, device = cairo_pdf)
ggsave("bench_f1_heatmap.png", p, width = 9, height = 5.5, dpi = 300)
cat("Saved bench_f1_heatmap.pdf and .png\n")