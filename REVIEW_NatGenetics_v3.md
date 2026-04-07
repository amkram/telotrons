# Nature Genetics Critical Review v3

**Verdict: Minor revisions**

The manuscript has improved substantially from v2. The metatranscriptomic P-value is now properly framed, the ortholog limitation is acknowledged, the within-Haptista Pif1 control is reported, strand denominators are unified, the CCCTAA splice-site problem is stated, and ED figures are sequentially numbered. The Data Availability section is now present. Most v2 issues are resolved. The remaining issues are primarily internal consistency problems and a few missing standard sections.

---

## Moderate Issues (3)

### 1. Abstract still overstates splice competence for ~44% of telotrons

The abstract states: "Splice competence is intrinsic: tandem TTAGGG arrays inherently encode GT donor and AG acceptor dinucleotides, enabling immediate intronization." The body now explicitly acknowledges (strand paragraph) that template-only telotrons (44.1% of all high-purity telotrons) cannot derive GT from the CCCTAA repeat. The abstract sentence remains unqualified and effectively claims all telotrons have intrinsic splice competence. This was flagged in v2 and addressed in the body text but not in the abstract. Suggest qualifying: "In converging telotrons, tandem TTAGGG arrays inherently encode GT donor and AG acceptor dinucleotides..." or "The TTAGGG repeat encodes splice signals that enable immediate intronization in the majority of configurations."

### 2. Multiple shifting MAG denominators without explanation

The manuscript uses several different MAG subsets without clearly defining transitions between them:
- **257 MAGs** — full Tara Oceans dataset
- **30 MAGs with ≥50 telomeric-repeat introns** — for boundary analysis
- **20 telotron-bearing MAGs** — for TERT detection ("5 of 20")
- **8 most affected MAGs** — for the 2,371 high-purity count
- **53 MAGs** — for burden–youth correlation
- **48 MAGs** — for Ku70/Pif1 correlations

What distinguishes "30 MAGs with ≥50 telomeric-repeat introns" from "20 telotron-bearing MAGs"? Why are there 53 MAGs in the burden–youth analysis but only 48 for repair gene correlations? A single sentence in Methods or main text establishing these nested subsets (e.g., "Of the 257 MAGs, 53 contained ≥N telomeric introns and were used for burden analyses; of these, 48 had sufficient gene model counts for repair gene density normalization; 30 had ≥50 telomeric-repeat introns for boundary signature analysis; and 8 harboured ≥100 high-purity telotrons") would resolve this.

### 3. Pif1 within-Haptista result (n = 7) is fragile

The within-Haptista Pif1 correlation (ρ = −0.82, P = 0.023, n = 7) is technically significant but extremely underpowered. With n = 7, a single outlier can drive the entire result. The manuscript should acknowledge this limitation explicitly — e.g., "though the small sample limits confidence" — rather than presenting the Ku70 and Pif1 within-Haptista results with equal weight. Consider whether the Pif1 result should be framed as "consistent with" rather than "significant."

---

## Minor Issues (10)

### 4. "strict" still appears in closing paragraph

The closing paragraph reads "844 strict telotrons and over 32,000 telo-containing introns." The user previously requested that the word "strict" not be used. Replace with "high-purity" for consistency.

### 5. 721 vs 844 denominator unexplained

The intronization paragraph gives "93.5% of telotrons in the primary MAG use canonical GT-AG splice sites (674 of 721)." The census paragraph establishes 844 high-purity telotrons. Why is the splice-site denominator 721 rather than 844? Presumably 721 are the telotrons with annotated splice dinucleotides at both boundaries, and 123 have ambiguous or non-standard sites. This should be stated.

### 6. "transcription-associated DSB repair²¹" — ref 21 is broad

Ref 21 (Aguilera & García-Muse, 2013) is a general review on genome instability. The specific claim about template-strand bias reflecting "transcription-associated DSB repair" would benefit from a more targeted reference on transcription-coupled DNA damage and R-loop-mediated breaks (e.g., Aguilera & Gómez-González, *Nat. Rev. Genet.* 2008, or Crossley et al., *Mol. Cell* 2019).

### 7. "BLASTp-like sequence comparison" — vague

The Methods state orthologs were identified by "BLASTp-like sequence comparison across MAGs using a custom pipeline." Either state it is BLASTp, or name the actual algorithm/tool. "BLASTp-like" is imprecise.

### 8. "improving sensitivity by ~35%" — unanchored claim

The repair gene Methods state motif-based approaches improved sensitivity "by ~35%." Relative to what baseline? Simple regex? HMM? BLAST? This claim needs a comparator.

### 9. Coding-strand-only splice competence unaddressed

The strand paragraph addresses converging (GT from TTAGGG array) and template-only (GT from flanking) splice mechanisms, but does not explain how coding-strand-only telotrons (TTAGGG on coding strand, n = 10) achieve splice competence. Presumably the TTAGGG repeat on the coding strand directly provides the GT donor. A brief statement would close the gap.

### 10. Data Availability BioProject accession may be incorrect

"BioProject PRJEB402" looks unusually short for a Tara Oceans eukaryotic MAG accession. The Delmont et al. (2022) Cell Genomics paper lists PRJEB402 as the umbrella project, but the specific MAG assemblies may be under a sub-accession (e.g., PRJEB52452 or similar). Verify and consider listing the specific sub-project for reproducibility.

### 11. Missing standard Nature Genetics sections

The manuscript lacks: (a) Acknowledgments, (b) Author Contributions, (c) Competing Interests declaration, and (d) Code Availability (currently merged into Data Availability — Nature Genetics typically separates these). These are all required for submission.

### 12. Metatranscriptomic length explanation is untested

The text states the higher telotron validation rate "likely reflects their shorter median length, since shorter introns are more efficiently captured by junction-spanning reads." This is presented as fact but is untested. Consider softening to "may reflect" or, better, testing the hypothesis by binning introns by length and comparing validation rates within length-matched subsets.

### 13. Tetrahymena parenthetical is clunky

"(which uses TTGGGG; the TBE mechanism is structurally conserved but boundary positions differ by species)" disrupts the flow of an important paragraph. Consider moving to a footnote-style aside or integrating more smoothly: "as shown by crystallography of the TERT–TBE complex in Tetrahymena⁸ — a species with a different repeat (TTGGGG) but a structurally conserved boundary mechanism."

---

## Summary of Required Changes

| # | Issue | Severity | Effort |
|---|-------|----------|--------|
| 1 | Abstract splice-competence claim overstated | Moderate | Low |
| 2 | MAG denominator nesting unexplained | Moderate | Low |
| 3 | Pif1 n=7 needs caveat | Moderate | Low |
| 4 | "strict" in closing paragraph | Minor | Trivial |
| 5 | 721 vs 844 denominator | Minor | Low |
| 6 | Ref 21 too broad | Minor | Low |
| 7 | "BLASTp-like" vague | Minor | Trivial |
| 8 | "~35% sensitivity" unanchored | Minor | Trivial |
| 9 | Coding-strand splice mechanism | Minor | Low |
| 10 | BioProject accession | Minor | Low |
| 11 | Missing standard sections | Minor | Low |
| 12 | Length explanation untested | Minor | Low |
| 13 | Tetrahymena parenthetical | Minor | Trivial |

**Overall assessment:** The manuscript is close to submission-ready. The science is strong and the evidence for telomerase-mediated intron creation is compelling across multiple independent lines. The remaining issues are mostly about precision and completeness rather than substance. One round of targeted revisions should bring it to final form.
