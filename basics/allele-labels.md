These are **two different ways of labeling alleles**: one describes the reference genome; the other describes a statistical association.

| Term                         | Meaning                                                                                       | Does it indicate the direction of an effect?                                                   |
| ---------------------------- | --------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------- |
| **Reference allele (REF)**   | The allele present at that position in a specified reference genome assembly, such as GRCh38. | No                                                                                             |
| **Alternative allele (ALT)** | An allele that differs from REF at that position. A site can have multiple ALT alleles.       | No                                                                                             |
| **Effect allele (EA)**       | The allele whose association is reported by the effect estimate, such as β or an odds ratio.  | Yes—the estimate describes the effect per additional copy of this allele in an additive model. |
| **Other allele (OA)**        | The allele against which EA is compared, usually in a biallelic analysis.                     | It provides the comparison.                                                                    |

**The effect allele can be either REF or ALT.** You need the dataset’s documentation to know which was used.

For example, suppose a SNP has **REF = A** and **ALT = G**:

| Reporting convention            | Effect allele | Other allele | β     |
| ------------------------------- | ------------- | ------------ | ----- |
| Association reported for G      | G             | A            | +0.20 |
| Same association reported for A | A             | G            | −0.20 |

These describe **the same association**. Switching EA and OA changes the sign of β. For odds ratios, it takes the reciprocal: **OR = 1.25 becomes OR = 0.80**. The standard error and two-sided P value remain unchanged.

A few important distinctions:

* **REF does not necessarily mean common, ancestral, healthy, or protective.**
* **ALT does not necessarily mean rare, harmful, or disease-causing.**
* **EA does not necessarily mean risk-increasing**: a negative β or OR < 1 indicates a decreasing effect on the analyzed outcome.
* Labels such as **A1/A2 are software-dependent**; do not assume A1 always means ALT or EA.

The key is to **align effect estimates to the same effect allele across datasets**, while also checking genome build and strand orientation. Matching REF/ALT alone is insufficient unless you know which allele each reported β refers to.

