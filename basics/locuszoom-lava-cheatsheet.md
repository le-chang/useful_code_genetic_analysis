# Cheatsheet: LocusZoom plots & visualizing LAVA results

For students working with the cross-trait pipeline (LDSC → conjFDR / LAVA → coloc).
Trait names are kept generic (**Trait A / Trait B**).

---

## TL;DR

| Question | Short answer |
|---|---|
| Which region should I plot first in LocusZoom? | Start from **coloc**, not LAVA. Rank regions by **PP.H4** and plot the top one first. |
| What's the best plot for LAVA results? | A **forest plot** of significant loci (local rho ± 95% CI). Add a **genome-wide overview** (Manhattan/Miami-style) for context. |

---

## 1. LAVA vs coloc: what each one tells you

| | LAVA (bivariate) | coloc (ABF) |
|---|---|---|
| Unit of analysis | Pre-defined LD block (~2,500 blocks genome-wide, ~1 Mb on average) | Region around a lead SNP |
| Question asked | Are SNP effects **correlated** across the whole block? | Do both traits share **one causal variant** here? |
| Key output | `rho` (with 95% CI), `p` | PP.H0–H4 (**PP.H4** = shared causal variant) |
| Needs a strong SNP-level peak? | No. It can be driven by many small effects. | Yes, in **both** traits. |

**Why this matters for LocusZoom:** LocusZoom is a SNP-level plot. A LAVA block can be highly significant without a clear peak, so ranking by the LAVA p-value often gives a flat, uninformative plot. Coloc results map directly onto what LocusZoom shows.

---

## 2. LocusZoom: which regions, in what order

1. **Filter** coloc results to PP.H4 ≥ 0.8 (strong evidence). Treat 0.5–0.8 as suggestive and label it that way.
2. **Rank** by PP.H4 and plot the **top region first** as a test run.
3. **Check robustness to the p12 prior.** The coloc default is p12 = 1e-5; also check a stricter value (e.g. 5e-6 or 1e-6), for example with `coloc::sensitivity()`. If a call only passes at one p12, report it with that caveat.
4. Once the first plot looks right, **plot the remaining strong regions**.
5. *Optional:* note whether each colocalized locus falls in a LAVA block with significant local rg. Use this as **supporting evidence**, not as the ranking criterion.

---

## 3. LocusZoom: making the plot

- **Index SNP:** centre on the SNP with the highest `SNP.PP.H4` (the candidate shared variant). If that SNP is missing in one trait, use the lead SNP of the trait with the stronger signal.
- **Window:** use the same boundaries as the coloc region, or about ±500 kb.
- **Stack both traits** on the same x-axis, and colour LD relative to the **same index SNP** in both panels.
- **Genome build:** both GWAS and the LD reference must use the same build. The default LAVA locus file is GRCh37/hg19.
- **LD reference ancestry** must match the GWAS (e.g. 1000 Genomes EUR for European-ancestry GWAS).
- **Add a LocusCompare scatter** (−log10 p Trait A vs −log10 p Trait B). This is often the clearest visual evidence of colocalization.

**Tools:** [my.locuszoom.org](https://my.locuszoom.org) (web upload), `locuszoomr` (R), `LocusCompareR` (R).

**A convincing plot:** both peaks sit on the same SNP with consistent LD colouring, and LocusCompare points rise along the diagonal toward the top-right.

**Red flags:**
- Peaks offset from each other: likely distinct causal variants (closer to H3 than H4).
- Weak signal in one trait: PP.H4 may be driven by the prior, so interpret with caution.
- Multiple independent peaks: this violates coloc's single-causal-variant assumption. Consider `coloc.susie`.

---

## 4. Visualizing LAVA results

### A. Forest plot (recommended main figure)
- One row per significant locus (or locus × trait pair).
- Point = local `rho`, bars = 95% CI (`rho.lower` / `rho.upper`), vertical line at 0.
- Order rows by chromosome and position. Label with locus ID and a nearby gene.
- Colour by trait pair if there are several.
- This shows **direction, effect size and uncertainty** together, which a p-value plot cannot.

### B. Genome-wide overview (supporting)
- Manhattan-style: LD blocks on x, −log10(p) on y, points coloured by the **sign of rho**.
- Draw the Bonferroni line for the **number of bivariate tests actually run**.
- Miami variant: positive correlations above the axis, negative below.

### C. Heatmap (when testing several trait pairs)
- Rows = loci, columns = trait pairs, fill = rho on a diverging scale centred at 0.
- Mark significant cells with `*`. Leave untested cells grey (loci that failed the univariate filter).

### Skip
- **Circos plots.** They look nice on posters, but direction and CIs are hard to read.

### Rules for all LAVA plots
- Plot only loci that **passed the univariate heritability filter**, and state the threshold in the legend. Otherwise a missing point looks like a null.
- If you have the LDSC global rg, add a **dashed line at that value** on the forest plot to show whether local signals agree with the genome-wide direction.
- `rho` near ±1 with very wide CIs usually reflects **low local heritability**, not a near-perfect correlation.

---

## 5. Minimal forest plot code (R)

```r
library(dplyr)
library(ggplot2)

# Assumes LAVA bivariate output with columns:
#   locus, chr, start, stop, phen1, phen2, rho, rho.lower, rho.upper, p
# If chr/start/stop are missing, join them from the LAVA locus file.
# Adjust column names to your pipeline's output.
bivar <- read.table("lava_bivar.tsv", header = TRUE)

n_tests <- nrow(bivar)   # number of bivariate tests actually run

sig <- bivar %>%
  filter(p < 0.05 / n_tests) %>%
  arrange(chr, start) %>%
  mutate(
    locus_label = paste0("chr", chr, ":",
                         round(start / 1e6, 1), "-", round(stop / 1e6, 1), " Mb"),
    locus_label = factor(locus_label, levels = rev(unique(locus_label))),
    pair = paste(phen1, "x", phen2)
  )

ggplot(sig, aes(x = rho, y = locus_label, colour = pair)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = rho.lower, xmax = rho.upper),
                  position = position_dodge(width = 0.5)) +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(x = "Local genetic correlation (rho, 95% CI)", y = NULL, colour = "Trait pair") +
  theme_bw()
```

---

## 6. Checklist before sharing figures

- [ ] LocusZoom regions come from coloc (PP.H4 ≥ 0.8), not from ranking by LAVA p-value
- [ ] Coloc calls checked for sensitivity to p12
- [ ] Same genome build across both GWAS, the LD reference and the locus file
- [ ] LD reference ancestry matches the GWAS
- [ ] Both LocusZoom panels coloured by LD to the same index SNP
- [ ] LAVA plots show only loci that passed the univariate filter, with the threshold stated
- [ ] Bonferroni correction uses the number of bivariate tests actually run
- [ ] No unpublished trait names in figures or files pushed to the repo
