# locuszoomr for colocalization results — cheatsheet

*How to draw a stacked regional plot (trait 1, trait 2, coloc posterior, genes) at one locus with `locuszoomr`. Written for lab members; the example uses generic trait names.*

---

## The one idea

**One `locus()` object per data frame, then stack the panels.**

`locus()` is the only function that turns *your* data into a track. The `link_*()` functions (`link_LD()`, `link_eqtl()`, `link_recomb()`) are API fetchers that pull *external* data (1000G LD, GTEx eQTLs, recombination rates) into an existing locus object. Coloc output is not external data — it's a data frame you already have — so it goes through `locus()`, not `link_eqtl()`.

| You have | Use | Not |
|---|---|---|
| GWAS summary stats (either trait) | `locus(data = ...)` | — |
| eQTL summary stats you ran yourself | `locus(data = ...)` | `link_eqtl()` |
| `coloc.abf()` per-SNP results | `locus(data = ..., yvar = "SNP.PP.H4")` | `link_eqtl()` |
| Want LD colouring | `link_LD(loc, token = <LDlink token>)` | — |
| Want GTEx eQTLs from the web | `link_eqtl(loc, token = <LDlink token>)` → `eqtl_plot()` | — |

The `token` in `link_eqtl()` / `link_LD()` is an **LDlink personal access token** (register at ldlink.nih.gov). It has nothing to do with coloc.

---

## What coloc output looks like

`res <- coloc.abf(dataset1, dataset2)` returns:

- `res$summary` — PP.H0 … PP.H4 for the region (one number each; goes in a table or plot title, not a track)
- `res$results` — one row per SNP **shared by both datasets**, with `snp`, `position`, per-dataset z / lABF, and **`SNP.PP.H4`** = posterior that this SNP is the shared causal variant, conditional on H4

`SNP.PP.H4` is the only per-SNP quantity worth a track. The figure most readers expect is the **two association tracks stacked at the same window** (that is the visual argument for colocalization), with `SNP.PP.H4` as an optional third panel.

---

## Full example (base graphics)

```r
library(locuszoomr)
library(EnsDb.Hsapiens.v75)   # GRCh37; use v86 / v105 for GRCh38

# 1. Same window on every panel — the EXACT region you passed to coloc,
#    not a gene-centred flank.
chr <- "7"
xr  <- c(start, end)

# 2. Track 1: trait-1 GWAS
loc1 <- locus(data = gwas1_sub, seqname = chr, xrange = xr,
              chrom = "chr", pos = "pos", p = "p", labs = "rsid",
              ens_db = "EnsDb.Hsapiens.v75")

# 3. Track 2: the other dataset you gave coloc (second GWAS, or eQTL)
loc2 <- locus(data = gwas2_sub, seqname = chr, xrange = xr,
              chrom = "chr", pos = "pos", p = "p", labs = "rsid",
              ens_db = "EnsDb.Hsapiens.v75",
              index_snp = loc1$index_snp)      # same lead SNP on every panel

# 4. Track 3: coloc per-SNP posterior. yvar= plots the raw value (0-1), no -log10.
pp <- res$results[, c("snp", "position", "SNP.PP.H4")]
pp$chr <- chr
loc3 <- locus(data = pp, seqname = chr, xrange = xr,
              chrom = "chr", pos = "position", yvar = "SNP.PP.H4", labs = "snp",
              ens_db = "EnsDb.Hsapiens.v75", index_snp = loc1$index_snp)

# 5. Optional: LD colouring (needs an LDlink token, internet access)
# loc1 <- link_LD(loc1, token = "LDLINK_TOKEN")
# loc2 <- link_LD(loc2, token = "LDLINK_TOKEN")

# 6. Stack: n scatter panels + gene track at the bottom
pdf("locus_coloc.pdf", width = 6, height = 9)
set_layers(3, heights = c(3, 3, 2, 2))
scatter_plot(loc1, xticks = FALSE, ylab = expression("Trait 1  -log"[10] ~ "P"))
scatter_plot(loc2, xticks = FALSE, ylab = expression("Trait 2  -log"[10] ~ "P"))
scatter_plot(loc3, xticks = FALSE, ylab = "SNP PP.H4")
genetracks(loc1)
dev.off()
```

Add `main = sprintf("PP.H4 = %.2f", res$summary["PP.H4.abf"])` to the top panel if you want the regional posterior on the figure.

---

## ggplot2 variant

```r
library(patchwork)

(gg_scatter(loc1, xticks = FALSE, ylab = expression("Trait 1  -log"[10] ~ "P")) /
 gg_scatter(loc2, xticks = FALSE, ylab = expression("Trait 2  -log"[10] ~ "P")) /
 gg_scatter(loc3, xticks = FALSE, ylab = "SNP PP.H4") /
 gg_genetracks(loc1)) +
  plot_layout(heights = c(3, 3, 2, 2))
```

Same locus objects, same rules; only the drawing layer changes.

---

## Checklist before you trust the figure

- [ ] `xrange` is the exact coloc window for all three panels.
- [ ] Every panel uses the same `index_snp` (`loc1$index_snp`), otherwise each panel highlights its own lead and the plot argues against you.
- [ ] Positions are in the same build as the `EnsDb` (v75 → GRCh37). Coloc's `position` column is whatever you fed in — check it.
- [ ] The PP.H4 panel has fewer points than the GWAS panels. Expected: `res$results` only contains SNPs present in both datasets.
- [ ] Chromosome names match between your data (`"7"` vs `"chr7"`) and the EnsDb (no `chr` prefix).
- [ ] If you ran `coloc.susie()` instead of `coloc.abf()`, plot the credible sets (colour/shape by `cs`), not a single PP.H4 track.

---

## Common mistakes

- Calling `link_eqtl()` on coloc output — it will try to hit the LDlink API and store nothing useful.
- Using a gene-centred `flank` for the plot but a different window for coloc, so the SNP with max PP.H4 sits outside the panel.
- Letting each panel pick its own index SNP.
- Forgetting `yvar =` for PP.H4, so `locus()` tries to -log10 a posterior probability.
- Mixing GRCh37 summary stats with a GRCh38 EnsDb (genes land in the wrong place, no error is raised).

---

## References

- Lewis MJ, Wang S. locuszoomr: an R package for visualizing publication-ready regional gene locus plots. *Bioinformatics Advances* 2025;5(1):vbaf006.
- Giambartolomei C, et al. Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. *PLoS Genet* 2014;10:e1004383.
- locuszoomr reference manual: https://cran.r-universe.dev/locuszoomr/doc/manual.html
- LDlink token registration: https://ldlink.nih.gov/?tab=apiaccess
