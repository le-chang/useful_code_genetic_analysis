#!/usr/bin/env Rscript
# =============================================================================
#  plot_lava_bivar.R
#  Plots for LAVA bivariate local genetic correlation output (*.bivar.lava)
#
#  Usage (command line):
#      Rscript plot_lava_bivar.R <input.bivar.lava> [output_dir]
#
#  Usage (RStudio / interactive):
#      source("plot_lava_bivar.R")      # loads the functions only
#      res <- run_all("analysis1_local_rg.bivar.lava", "lava_plots")
#      # or work with the pieces:
#      lava <- tidy_lava("analysis1_local_rg.bivar.lava")
#      plot_genome(lava, pairs = c("qtA vs dzX", "qtB vs dzX"))
#
#  Outputs written to output_dir:
#      lava_bivar_tidy.tsv            input + pair labels, adjusted p, significance
#      summary_per_pair.tsv           n tested / n significant (+ / -) per pair
#      significant_cross_category.tsv significant cross-category loci, sorted by p
#      01_pair_counts.{png,pdf}       significant loci per pair, split by sign
#      02_genome_cross.{png,pdf}      genome-wide local rho, cross-category pairs
#      02_genome_all.pdf              genome-wide local rho, every pair
#      03_heatmap_cross.{png,pdf}     locus x pair heatmap, cross-category pairs
#      03_heatmap_all.pdf             locus x pair heatmap, every pair
#      04_forest_cross.{png,pdf}      rho with 95% CI per locus, cross-category pairs
#      05_volcano.{png,pdf}           rho vs -log10(p) for all tests
#
#  Requires: dplyr, tidyr, ggplot2, scales, ggrepel
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(ggrepel)
})

## ---- 1. Settings ------------------------------------------------------------

alpha        <- 0.05
p_adjust     <- "bonferroni"  # "bonferroni" = alpha / n bivariate tests (LAVA default); or "BH"
genome_build <- "hg19"        # "hg19" or "hg38"; LAVA's default 2,495 blocks are hg19
neglog10_cap <- 50            # LAVA reports p = 0 for some loci; cap -log10(p) here for plotting

# How to turn a phenotype name into a trait group, and how groups map to
# categories. This is only used to split pairs into "within-category" vs
# "cross-category", so the cross-trait pairs can be plotted on their own.
# EDIT these two lines to match your phenotype names (see README).
trait_group    <- function(x) sub("_.*$", "", x)          # strip a suffix: "qtA_male" -> "qtA"
group_category <- c(qtA = "quantitative", qtB = "quantitative", qtC = "quantitative",
                    dzX = "disease",      dzY = "disease")  # groups not listed get their own category
# group_category <- NULL                                   # <- use this if every group is its own category

# Colours
col_pos <- "#C2403C"   # significant positive rho
col_neg <- "#2C6FB7"   # significant negative rho
col_ns  <- "grey70"    # not significant
cls_cols <- c("positive" = col_pos, "negative" = col_neg, "n.s." = col_ns)

## ---- 2. Helpers ---------------------------------------------------------------

chr_lengths <- list(
  hg19 = c(`1` = 249250621, `2` = 243199373, `3` = 198022430, `4` = 191154276,
           `5` = 180915260, `6` = 171115067, `7` = 159138663, `8` = 146364022,
           `9` = 141213431, `10` = 135534747, `11` = 135006516, `12` = 133851895,
           `13` = 115169878, `14` = 107349540, `15` = 102531392, `16` = 90354753,
           `17` = 81195210, `18` = 78077248, `19` = 59128983, `20` = 63025520,
           `21` = 48129895, `22` = 51304566, X = 155270560),
  hg38 = c(`1` = 248956422, `2` = 242193529, `3` = 198295559, `4` = 190214555,
           `5` = 181538259, `6` = 170805979, `7` = 159345973, `8` = 145138636,
           `9` = 138394717, `10` = 133797422, `11` = 135086622, `12` = 133275309,
           `13` = 114364328, `14` = 107043718, `15` = 101991189, `16` = 90338345,
           `17` = 83257441, `18` = 80373285, `19` = 58617616, `20` = 64444167,
           `21` = 46709983, `22` = 50818468, X = 156040895)
)

chr_table <- function(build = genome_build) {
  len <- chr_lengths[[build]]
  tibble(chr = names(len), len = unname(len)) %>%
    mutate(chr_idx = row_number(),
           offset  = cumsum(c(0, len[-n()])),
           center  = offset + len / 2)
}

theme_lava <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")
}

save_plot <- function(p, path_stem, width, height, png = TRUE, pdf = TRUE, dpi = 300) {
  height <- min(height, 48); width <- min(width, 48)     # keep within device limits
  if (png) ggsave(paste0(path_stem, ".png"), p, width = width, height = height, dpi = dpi, limitsize = FALSE)
  if (pdf) ggsave(paste0(path_stem, ".pdf"), p, width = width, height = height, limitsize = FALSE)
  invisible(p)
}

## ---- 3. Read and tidy ----------------------------------------------------------

tidy_lava <- function(infile, build = genome_build) {
  d <- read.table(infile, header = TRUE, stringsAsFactors = FALSE,
                  na.strings = c("NA", ""), check.names = FALSE)
  need <- c("locus", "chr", "start", "stop", "phen1", "phen2", "rho", "rho.lower", "rho.upper", "p")
  miss <- setdiff(need, names(d))
  if (length(miss)) stop("Missing columns in ", infile, ": ", paste(miss, collapse = ", "))

  ct <- chr_table(build)
  n_tests  <- sum(!is.na(d$p))
  p_thresh <- if (p_adjust == "bonferroni") alpha / n_tests else NA_real_

  g1 <- trait_group(d$phen1); g2 <- trait_group(d$phen2)
  cat_of <- function(g) {
    if (is.null(group_category)) return(g)
    cc <- unname(group_category[g]); cc[is.na(cc)] <- g[is.na(cc)]; cc
  }

  d <- d %>%
    mutate(chr   = sub("^chr", "", as.character(chr)),
           pair  = paste(phen1, phen2, sep = " vs "),
           pair_type = ifelse(cat_of(g1) == cat_of(g2), "within-category", "cross-category"),
           locus_label = sprintf("chr%s:%.2f-%.2f Mb (locus %s)", chr, start / 1e6, stop / 1e6, locus),
           p_adj = if (p_adjust == "bonferroni") pmin(p * n_tests, 1) else p.adjust(p, method = "BH"),
           sig   = !is.na(p_adj) & p_adj < alpha,
           class = factor(ifelse(!sig, "n.s.", ifelse(rho > 0, "positive", "negative")),
                          levels = c("positive", "negative", "n.s.")),
           neglog10p = pmin(-log10(p), neglog10_cap)) %>%
    left_join(ct %>% select(chr, chr_idx, offset), by = "chr") %>%
    mutate(pos_mid = offset + (start + stop) / 2)

  if (any(is.na(d$chr_idx)))
    warning("Chromosomes not in the ", build, " table were dropped from genome plots: ",
            paste(unique(d$chr[is.na(d$chr_idx)]), collapse = ", "))

  # pair ordering: cross-category first, then by number of tests
  lev <- d %>% count(pair, pair_type) %>%
    arrange(pair_type != "cross-category", desc(n)) %>% pull(pair)
  d$pair <- factor(d$pair, levels = lev)

  attr(d, "n_tests")  <- n_tests
  attr(d, "p_thresh") <- p_thresh
  attr(d, "build")    <- build
  message(sprintf("Read %d bivariate tests over %d loci and %d trait pairs; %d significant (%s, alpha = %g%s).",
                  nrow(d), n_distinct(d$locus), n_distinct(d$pair), sum(d$sig), p_adjust, alpha,
                  if (!is.na(p_thresh)) sprintf(", p < %.3g", p_thresh) else ""))
  d
}

summarise_pairs <- function(dat) {
  dat %>%
    group_by(pair, pair_type) %>%
    summarise(n_tested = n(),
              n_sig    = sum(sig),
              n_pos    = sum(sig & rho > 0),
              n_neg    = sum(sig & rho < 0),
              min_p    = min(p, na.rm = TRUE),
              .groups  = "drop") %>%
    arrange(pair)
}

## ---- 4. Plots -------------------------------------------------------------------

# 4a. Significant loci per pair, by sign ---------------------------------------
plot_pair_counts <- function(dat) {
  s <- summarise_pairs(dat)
  long <- s %>%
    select(pair, pair_type, n_tested, n_sig, positive = n_pos, negative = n_neg) %>%
    pivot_longer(c(positive, negative), names_to = "class", values_to = "n") %>%
    mutate(class = factor(class, levels = c("positive", "negative")))
  ggplot(long, aes(x = n, y = pair, fill = class)) +
    geom_col(width = 0.75) +
    geom_text(data = s, aes(x = n_sig, y = pair, label = sprintf("%d / %d", n_sig, n_tested)),
              inherit.aes = FALSE, hjust = -0.15, size = 3, colour = "grey30") +
    scale_fill_manual(values = cls_cols[1:2], name = "Significant local rho") +
    scale_y_discrete(limits = rev) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
    facet_grid(pair_type ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Number of significant loci (label: significant / tested)", y = NULL) +
    theme_lava() + theme(panel.grid.major.x = element_line(colour = "grey92"),
                         panel.grid.major.y = element_blank())
}

# 4b. Genome-wide local rho, one panel per pair ----------------------------------
plot_genome <- function(dat, pairs = NULL, ncol = 3, show_ci = TRUE, point_size = c(1.2, 2.2)) {
  if (!is.null(pairs)) dat <- dat %>% filter(pair %in% pairs)
  dat <- dat %>% filter(!is.na(chr_idx)) %>% droplevels()
  ct  <- chr_table(attr(dat, "build") %||% genome_build) %>%
    filter(chr_idx <= max(dat$chr_idx)) %>%
    mutate(band  = rep(c("a", "b"), length.out = n()),
           label = ifelse(chr_idx > 16 & chr_idx %% 2 == 0, "", chr))   # thin crowded labels
  ggplot(dat) +
    geom_rect(data = ct, aes(xmin = offset, xmax = offset + len, ymin = -Inf, ymax = Inf, fill = band),
              inherit.aes = FALSE) +
    scale_fill_manual(values = c(a = "white", b = "grey93"), guide = "none") +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey45") +
    { if (show_ci) geom_linerange(aes(x = pos_mid, ymin = rho.lower, ymax = rho.upper, colour = class),
                                  alpha = 0.45, linewidth = 0.4) } +
    geom_point(aes(x = pos_mid, y = rho, colour = class, size = sig)) +
    scale_colour_manual(values = cls_cols, name = NULL, drop = FALSE) +
    scale_size_manual(values = c(`FALSE` = point_size[1], `TRUE` = point_size[2]), guide = "none") +
    scale_x_continuous(breaks = ct$center, labels = ct$label, expand = c(0.005, 0)) +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
    facet_wrap(~ pair, ncol = ncol) +
    labs(x = "Chromosome", y = expression("Local " * r[g] * " (" * rho * ", 95% CI)")) +
    theme_lava() +
    theme(axis.text.x = element_text(size = 7), panel.spacing.x = unit(0.8, "lines"))
}

# 4c. Locus x pair heatmap ---------------------------------------------------------
plot_heatmap <- function(dat, pairs = NULL, loci_sig_only = TRUE, star_size = 3.5) {
  if (!is.null(pairs)) dat <- dat %>% filter(pair %in% pairs)
  if (loci_sig_only) dat <- dat %>% filter(locus %in% locus[sig])
  dat <- dat %>% droplevels() %>%
    arrange(chr_idx, start) %>%
    mutate(locus_label = factor(locus_label, levels = rev(unique(locus_label))))
  ggplot(dat, aes(x = pair, y = locus_label, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(sig, "*", "")), size = star_size, vjust = 0.8) +
    scale_fill_gradient2(low = col_neg, mid = "white", high = col_pos, limits = c(-1, 1),
                         name = expression("Local " * rho)) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL,
         caption = "* significant after multiple-testing correction\nblank = not tested (no univariate signal in one of the traits)") +
    theme_lava() +
    theme(axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0),
          axis.text.y = element_text(size = 7),
          panel.grid = element_blank(), panel.border = element_blank(),
          legend.position = "right", plot.caption = element_text(hjust = 0))
}

# 4d. Forest plot: rho with 95% CI per locus ---------------------------------------
plot_forest <- function(dat, pairs = NULL, sig_only = FALSE, ncol = 3) {
  if (!is.null(pairs)) dat <- dat %>% filter(pair %in% pairs)
  if (sig_only) dat <- dat %>% filter(sig)
  dat <- dat %>% droplevels() %>%
    arrange(chr_idx, start) %>%            # genomic order, shared across panels
    mutate(locus_label = factor(locus_label, levels = rev(unique(locus_label))))
  ggplot(dat, aes(x = rho, y = locus_label, colour = class)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey45") +
    geom_errorbar(aes(xmin = rho.lower, xmax = rho.upper), width = 0.35, linewidth = 0.5) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = cls_cols, name = NULL, drop = FALSE) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
    facet_wrap(~ pair, scales = "free_y", ncol = ncol) +
    labs(x = expression("Local " * r[g] * " (" * rho * ", 95% CI)"), y = NULL) +
    theme_lava() + theme(axis.text.y = element_text(size = 7))
}

# 4e. Volcano: rho vs -log10(p) -----------------------------------------------------
plot_volcano <- function(dat, label_n = 12, label_cross_only = TRUE) {
  p_thresh <- attr(dat, "p_thresh")
  lab <- dat %>% filter(sig, if (label_cross_only) pair_type == "cross-category" else TRUE) %>%
    arrange(p) %>% head(label_n) %>%
    mutate(lab = sprintf("%s\nchr%s:%.1f Mb", pair, chr, start / 1e6))
  dat <- dat %>% arrange(pair_type == "cross-category")   # draw cross-category points on top
  ggplot(dat, aes(x = rho, y = neglog10p)) +
    { if (!is.na(p_thresh)) geom_hline(yintercept = -log10(p_thresh), linetype = 2, colour = "grey45") } +
    geom_vline(xintercept = 0, colour = "grey85") +
    geom_point(aes(colour = class, shape = pair_type), alpha = 0.75, size = 1.8) +
    geom_label_repel(data = lab, aes(label = lab), size = 2.4, lineheight = 0.9,
                     fill = alpha("white", 0.9), label.size = 0.15,
                     nudge_x = -0.35 * sign(lab$rho), xlim = c(-1, 1),   # push labels toward the centre
                     box.padding = 0.6, max.overlaps = Inf, segment.colour = "grey40",
                     min.segment.length = 0, seed = 1) +
    scale_colour_manual(values = cls_cols, name = NULL) +
    scale_shape_manual(values = c(`cross-category` = 17, `within-category` = 16), name = NULL) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
    labs(x = expression("Local " * r[g] * " (" * rho * ")"),
         y = bquote(-log[10] * "(p)" ~ .(if (any(dat$neglog10p >= neglog10_cap)) sprintf("(capped at %d)", neglog10_cap) else ""))) +
    theme_lava()
}

## ---- 5. Run everything ----------------------------------------------------------

run_all <- function(infile, outdir = "lava_plots", build = genome_build) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  out <- function(f) file.path(outdir, f)

  lava  <- tidy_lava(infile, build)
  summ  <- summarise_pairs(lava)
  cross <- levels(droplevels(lava$pair[lava$pair_type == "cross-category"]))
  n_cross_loci <- n_distinct(lava$locus[lava$pair_type == "cross-category"])

  write.table(lava %>% select(-chr_idx, -offset, -pos_mid), out("lava_bivar_tidy.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(summ, out("summary_per_pair.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(lava %>% filter(sig, pair_type == "cross-category") %>% arrange(p) %>%
                select(locus, chr, start, stop, n_snps, phen1, phen2, rho, rho.lower, rho.upper, r2, p, p_adj),
              out("significant_cross_category.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  # 01 counts
  save_plot(plot_pair_counts(lava), out("01_pair_counts"),
            width = 8, height = 2 + 0.26 * nlevels(lava$pair))

  # 02 genome
  if (length(cross)) {
    nc <- min(3, length(cross))
    save_plot(plot_genome(lava, pairs = cross, ncol = nc), out("02_genome_cross"),
              width = 4.5 * nc + 1, height = 1.6 + 1.8 * ceiling(length(cross) / nc))
  }
  save_plot(plot_genome(lava, ncol = 3), out("02_genome_all"), png = FALSE,
            width = 14.5, height = 1.6 + 1.8 * ceiling(nlevels(lava$pair) / 3))

  # 03 heatmap
  if (length(cross)) {
    save_plot(plot_heatmap(lava, pairs = cross, loci_sig_only = FALSE), out("03_heatmap_cross"),
              width = max(6.5, 3.5 + 0.32 * length(cross)), height = 2.5 + 0.22 * n_cross_loci)
  }
  n_loci_sig <- n_distinct(lava$locus[lava$sig])
  save_plot(plot_heatmap(lava, loci_sig_only = TRUE), out("03_heatmap_all"), png = FALSE,
            width = max(6.5, 3.5 + 0.32 * nlevels(lava$pair)), height = 2.5 + 0.18 * n_loci_sig)

  # 04 forest (panel height scales with the largest number of loci in any pair)
  if (length(cross)) {
    nc <- min(3, length(cross))
    max_rows <- max(table(droplevels(lava$pair[lava$pair_type == "cross-category"])))
    save_plot(plot_forest(lava, pairs = cross, ncol = nc), out("04_forest_cross"),
              width = 5 * nc + 1,
              height = 1.6 + (0.7 + 0.18 * max_rows) * ceiling(length(cross) / nc))
  }

  # 05 volcano
  save_plot(plot_volcano(lava), out("05_volcano"), width = 8, height = 6.5)

  message("Written to ", normalizePath(outdir))
  invisible(list(data = lava, summary = summ))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

if (!interactive()) {
  args   <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) stop("Usage: Rscript plot_lava_bivar.R <input.bivar.lava> [output_dir]")
  infile <- args[1]
  outdir <- if (length(args) >= 2) args[2] else "lava_plots"
  run_all(infile, outdir)
}
