# Multiple-testing correction for LDSC genetic correlations — cheatsheet

*For cross-trait and sex-stratified designs. Written for lab members; adapt the worked example to your own traits.*

---

## The one rule

**Correct for the number of hypotheses you pre-specified, not the number of rows the software gives you.**

LDSC (and most wrappers) computes r_g for every pair of input files. With *n* summary-statistic files you get *n(n−1)/2* rows. Most of those pairs are not questions you are asking.

---

## Why you get more rows than you expect

| Input files (*n*) | Pairwise rows |
|---|---|
| 5 | 10 |
| 8 | 28 |
| 10 | 45 |
| 12 | 66 |

The row count is a property of the software, not of your study design.

---

## Worked example: sex-stratified cross-trait analysis

Design: 3 exposure traits (X1, X2, X3) × 2 disease outcomes (D1, D2), each sex-stratified (F/M) → 10 files → **45 rows**.

Hypothesis: *within each sex, is each exposure trait genetically correlated with each outcome?*

| Pair type | Count | Hypothesis? | Why / why not |
|---|---|---|---|
| Exposure–outcome, **same sex** (X1-F/D1-F, X1-M/D1-M, …) | **12** | **Yes** | This is the question |
| Exposure–outcome, **cross-sex** (X1-F/D1-M, …) | 12 | No | Mixes sex differences in two traits at once; not interpretable without a specific hypothesis |
| Exposure–exposure (X1/X2/X3 × 4 sex combinations) | 12 | No | Nuisance output |
| Outcome–outcome (D1/D2 × 4 sex combinations) | 4 | No | Nuisance output |
| Same trait, cross-sex (X1-F/X1-M, …, D2-F/D2-M) | 5 | Different question | QC check; the interesting null is r_g = 1, not r_g = 0 |
| **Total** | **45** | | |

**Family size m = 12 → α = 0.05 / 12 ≈ 4.2 × 10⁻³.**

---

## Checklist: defining the family

1. **Pre-specify.** Write the list of pairs down *before* looking at results. Filtering to the "interesting" pairs after seeing all rows is post hoc selection — a reviewer can then fairly insist you correct for every row.
2. **Count every analysis you will make an inferential claim from.**
   - Sex-combined LDSC reported alongside sex-stratified? Either add it to the family (12 + 6 = 18) or declare sex-combined primary and sex-stratified secondary, each with its own threshold. Decide up front.
   - Testing whether r_g differs between sexes? That is a separate family (one test per exposure–outcome pair, e.g. 6).
3. **One family per analysis layer.** LDSC global r_g, LAVA local r_g, conjFDR, colocalization each get their own correction, sized by what is tested in that layer (LAVA: number of loci × number of pairs, or per-pair, as pre-specified).
4. **Don't mix nulls.** Tests of r_g ≠ 0 and tests of r_g ≠ 1 answer different questions and should not share a family.

---

## Useful tests

**Genetic correlation differs from zero** (LDSC default): *z* = r_g / SE, *p* from the standard normal.

**Same trait across sexes, is r_g < 1?**
*z* = (1 − r_g) / SE. Reported as a sex-similarity check; use a one-sided or two-sided *p* as pre-specified.

**Does r_g differ between sexes** (independent samples)?
*z* = (r_g,F − r_g,M) / √(SE_F² + SE_M²)

---

## Bonferroni vs alternatives

- **Bonferroni** (α / m): standard in cross-trait genetics papers; reviewers accept it. Conservative when tests are correlated (exposure traits share genetics; F and M GWAS of the same disease share most of theirs).
- **Benjamini–Hochberg FDR**: reasonable companion when you want to be less conservative; report as a sensitivity analysis, not a replacement.
- **Effective number of tests** (Li & Ji 2005): estimates *m_eff* from the eigenvalues of the correlation matrix among tests; defensible but adds a step reviewers may question.

Recommended reporting: Bonferroni as primary; call *p* < 0.05 "nominal" or "suggestive"; state the threshold explicitly in Methods and in every table legend.

---

## Methods sentence template

> Genetic correlations were estimated with LDSC for the *N* pre-specified exposure–outcome pairs within each sex (*m* = 12). Statistical significance was defined as Bonferroni-corrected *P* < 0.05/12 = 4.2 × 10⁻³; *P* < 0.05 was considered nominally significant. Other pairwise correlations produced by the software (cross-sex and within-category pairs) were not hypotheses of this study and are reported in Supplementary Data *X* for completeness only.

---

## Minimal code pattern (R)

```r
# rg: LDSC output with columns p1, p2, rg, se, p
# Keep only the pre-specified pairs; do not compute alpha from nrow(rg)
keep <- with(rg,
  trait(p1) %in% exposures & trait(p2) %in% outcomes & sex(p1) == sex(p2))
tests <- rg[keep, ]

m      <- nrow(tests)             # should equal the number you pre-specified
alpha  <- 0.05 / m
tests$bonf_sig <- tests$p < alpha
tests$fdr_q    <- p.adjust(tests$p, method = "BH")   # sensitivity only
```

`trait()` and `sex()` are placeholders for however your file names encode trait and sex — check that `m` matches your pre-registered count before looking at any *p*-values.

---

## Common mistakes

- Dividing α by the number of output rows instead of the number of hypotheses.
- Choosing the "relevant" pairs after seeing the results.
- Forgetting to count the sex-combined analysis when it is also reported.
- Putting r_g ≠ 0 and r_g ≠ 1 tests in the same family.
- Using one α for LDSC and silently reusing it for LAVA or colocalization.

---

## References

- Bulik-Sullivan B, et al. An atlas of genetic correlations across human diseases and traits. *Nat Genet* 2015;47:1236–1241.
- Benjamini Y, Hochberg Y. Controlling the false discovery rate. *J R Stat Soc B* 1995;57:289–300.
- Li J, Ji L. Adjusting multiple testing in multilocus analyses using the eigenvalues of a correlation matrix. *Heredity* 2005;95:221–227.
