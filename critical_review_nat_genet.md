# Critical Review: Nature Genetics Manuscript Assessment

**Manuscript:** "Widespread de novo intron creation by telomerase-mediated double-strand break repair in marine eukaryotes"

**Review type:** Pre-submission internal review simulating Nature Genetics referee standards

**Reviewer approach:** Each point backed by specific literature citations where applicable

---

## 1. FORMAT AND SCOPE: Wrong submission category

**Severity: Must resolve before submission**

The manuscript in its current form does not fit the Nature Genetics Brief Communication format. Brief Communications require 1,000–1,500 words of main text (including abstract, references, and figure legends), a 3-sentence abstract of ≤70 words, a title of ≤10 words, ≤2 display items, and ~20 references. The current manuscript has:

- **Main text:** ~2,132 words (exceeds the 1,500-word ceiling by >40%)
- **Abstract:** ~142 words (double the 70-word limit)
- **Title:** 13 words ("Widespread de novo intron creation by telomerase-mediated double-strand break repair in marine eukaryotes"), exceeding the 10-word limit
- **Methods:** ~2,232 words (Brief Communications typically have short online methods, not standalone methods sections this long)

**Recommendation:** This is an Article, not a Brief Communication. The scope (new mechanism of intron gain, multi-supergroup survey, mechanistic model with G4 feed-forward, repair gene depletion, Introner exclusion) is far too broad for a Brief Communication. Submit as an Article (main text typically 3,000–5,000 words) or drastically cut. If staying as a Brief Communication, the repair gene correlation analysis, the G4 analysis, and the Introner exclusion argument would all need to move entirely to supplementary, and the main text would need to be halved.

---

## 2. ASSEMBLY ARTIFACTS: The central vulnerability

**Severity: High — reviewers will raise this first**

The entire study relies on metagenome-assembled genomes (MAGs), which are known to harbour chimeric contigs, especially for eukaryotes. Eukaryotic MAGs are far less mature than prokaryotic MAGs. The EukCC tool (Saary et al., *Genome Biology* 2020) was specifically developed because existing quality metrics systematically underperform for eukaryotic bins, and contamination rates can be substantial. The Tara Oceans eukaryotic MAG catalogue itself acknowledges incompleteness and chimeric assembly as ongoing challenges.

**Specific concerns:**

(a) Telomeric repeats are among the most difficult sequences for short-read assemblers. TTAGGG arrays can cause misassembly by collapsing non-homologous loci that share telomeric sequence, creating artifactual chimeric contigs where telomeric sequence appears internal. The manuscript claims telotrons are "uniformly distributed along contigs rather than clustering near boundaries," but this test is necessary but not sufficient — chimeric joins can occur anywhere along contigs, not only at boundaries.

(b) The manuscript states assembly validation includes five predictions (Extended Data Fig. 3), but does not report BUSCO completeness scores, EukCC contamination estimates, or read-depth consistency across telotron-containing contigs. A reviewer will ask: what is the BUSCO score for TARA_PSW_86_MAG_00284? If it's <50% (common for eukaryotic MAGs), the assembly quality concern becomes acute.

(c) The metatranscriptomic validation (46.7% vs 41.8% junction-read support) is the strongest rebuttal, but the manuscript should quantify how many telotrons have ≥2 independent junction reads (not just ≥1), since single reads could reflect mapping artifacts at repetitive loci.

**Recommendation:** Report BUSCO/EukCC scores for all 30 telotron-bearing MAGs in a supplementary table. Show read-depth consistency across telotron-containing regions. Strengthen metatranscriptomic validation by reporting the distribution of junction-read counts per telotron rather than just the binary validated/unvalidated rate.

---

## 3. PRIOR WORK ON ITS AT DSBs: Critical missing citations

**Severity: High — appears to be a significant oversight**

The manuscript frames telotrons as "a mechanism of de novo intron creation at genomic scale not previously observed" and states ITS "have been widely documented in other contexts." However, the specific mechanism proposed — telomerase inserting TTAGGG repeats at double-strand breaks — has direct precedent in mammalian genomics that the manuscript does not cite:

- **Nergadze et al., *Genome Research* 14:1704–1710 (2004):** Compared 10 human ITS with orthologs in 12 primate species. Showed that 9 ITS insertion events (dated 40–6 Ma) arose suddenly during primate evolution via DSB repair, not by gradual TTAGGG expansion. The ancestral sequences were interrupted precisely by telomeric repeats with typical DSB repair modifications (short deletions, random additions, duplications). The authors explicitly proposed "capture of telomeric DNA at break sites OR telomerase-mediated insertion."

- **Nergadze et al., *Genome Biology* 8:R260 (2007):** Extended the analysis to rodent genomes and provided direct evidence that telomerase RNA retrotranscription contributed to DSB repair during mammalian genome evolution.

These papers describe essentially the same mechanism the manuscript proposes — telomerase acting at DSBs — but in a non-intronic context. The novelty of telotrons is that the insertions land in coding regions and become spliced introns. This is an important distinction, but the manuscript needs to properly contextualize the mechanism as an extension of known ITS biology rather than claiming it as entirely unprecedented. Reviewers familiar with the ITS literature will see this omission immediately.

**Recommendation:** Cite Nergadze et al. 2004 and 2007 explicitly. Frame telotrons as demonstrating that this previously documented ITS-at-DSB mechanism can also create functional spliceosomal introns at genomic scale — that is the true novelty, not the mechanism itself.

---

## 4. MICROSATELLITE EXPANSION: The alternative hypothesis needs stronger exclusion

**Severity: High — this is the primary competing explanation**

The manuscript devotes one sentence to excluding microsatellite expansion: "telotron arrays maintain constant purity independent of array length, whereas slippage-expanded microsatellites accumulate mutations proportional to expansion time." This is a critical claim that needs more rigorous treatment.

**Issues:**

(a) The purity–length independence argument assumes a specific mutation model (clock-like point mutations during expansion). However, if expansion is rapid and recent (as the "burden–youth" correlation suggests), there may not have been time for mutations to accumulate regardless of mechanism. A young microsatellite expansion would also show high purity.

(b) Telomeric repeats (TTAGGG)n are among the most unstable microsatellite sequences known, highly prone to slippage-induced rearrangements (reviewed in Aksenova & Mirkin, *Genes* 2019). The ITS instability literature shows these sequences are hotspots for chromosomal fragility. This makes slippage a biologically plausible alternative that cannot be dismissed with a single sentence.

(c) The "terminal partial repeats show TTAGGG-specific enrichment expected from template-boundary dissociation" argument is compelling but is buried in the text. This should be a prominent, quantified finding — it is the strongest discriminator between telomerase and slippage, since slippage would produce random truncation points rather than template-boundary-specific termini.

**Recommendation:** Expand the microsatellite exclusion into a dedicated paragraph or supplementary analysis. Quantify the expected vs. observed terminal repeat distributions under both models. Consider whether long-read sequencing data from any Tara Oceans samples could independently confirm the arrays.

---

## 5. REPAIR GENE CORRELATIONS: Phylogenetic confounding remains a concern

**Severity: Medium-High**

The repair gene analysis is one of the manuscript's most provocative claims, but it has a structural weakness: the global correlation between Ku70 and telotron burden is flat (ρ ≈ −0.001 across all 228 MAGs), and the signal comes entirely from within Haptista (n = 65), where 60 of 65 MAGs are Prymnesiophyceae. The manuscript reports permutation tests (P < 0.0001) and partial correlations, but there is a deeper issue: the signal may reflect a Prymnesiophyceae-specific trait (e.g., a founder effect where the ancestor of this clade lost repair capacity for reasons unrelated to telotrons) rather than a causal relationship between repair depletion and telotron accumulation.

**Supporting literature:** Rijal et al. (*PLOS ONE* 2025; 20(3):e0308593) analyzed the evolutionary history of Ku proteins across eukaryotes and found that the Ku70 SAP domain is absent in Haptophyta specifically (also Discoba, Amoebozoa, and some Fungi and Alveolata). This is actually supportive of the manuscript's claim, but it also means the Ku70 depletion may be an ancient Haptophyta trait predating telotron accumulation — consistent with the manuscript's own statement that "depletion predates telotron accumulation" but undermining the causal narrative.

**Additional concerns:**

(a) The motif-based repair gene search (six-frame translation, dual-motif approach) is creative but non-standard. A reviewer will question sensitivity and specificity. How many known Ku70/Lig4/Rad51 proteins from model organisms (human, yeast, Arabidopsis) are recovered by these motifs? A benchmarking analysis against curated orthologs would strengthen confidence.

(b) The Pif1 null result (ρ = −0.006) is interesting but the interpretation ("Pif1 regulation operating post-translationally via Mec1 phosphorylation rather than at the copy-number level") is speculative. An equally valid interpretation is that the Pif1 motif search has poor sensitivity and the null result reflects detection failure.

**Recommendation:** Cite Rijal et al. 2025 explicitly — it independently confirms Ku70 structural divergence in Haptophyta. Add benchmarking of the motif search against known orthologs from reference genomes. Temper the causal language around repair gene depletion; "associated with" rather than "implicates" would be more appropriate given the correlational nature of the evidence.

---

## 6. INTRONER EXCLUSION: Rigorous but potentially over-argued

**Severity: Low-Medium**

The four-test Introner exclusion (TIR, TSD, length CV, pairwise identity) is thorough and convincing, but the manuscript spends substantial text on it. A Nature Genetics reviewer will note that introners and telotrons are so structurally different (one is a TE, the other is a tandem repeat array) that extensive exclusion feels like arguing against a strawman. The 99.1% hexamer repeat content alone makes introner origin implausible at a glance.

**However:** There is a subtle issue with the TSD analysis. The manuscript reports that 71% of telotrons with intron-free orthologs have TSDs of 3–21 bp, but the Introner exclusion states only 10.6% show TSDs ≥2 bp. These appear contradictory. The difference is that the ortholog-based TSDs are detected by alignment comparison (looking at flanking exonic sequence), while the Introner diagnostic looks at immediately flanking sequence identity. This should be explicitly reconciled, as a reviewer will flag the apparent inconsistency.

**Recommendation:** Condense the Introner exclusion to 2–3 sentences in the main text and move the full four-test analysis to Extended Data. Reconcile the TSD discrepancy between the two analyses.

---

## 7. STATISTICAL REPORTING AND MULTIPLE TESTING

**Severity: Medium**

(a) The manuscript reports many P-values from different analyses (Fisher's exact, Spearman, χ², Wilcoxon, permutation, KS, logistic regression, mixed-effects model) but does not apply a global multiple testing correction. The Methods state "Multiple comparisons in the G4 analyses were controlled using Benjamini–Hochberg FDR at q = 0.05" but no correction is applied to the repair gene correlations (4 genes × multiple supergroups × partial correlations). With Lig4 ρ = −0.58, Ku70 ρ = −0.44, and Rad51 ρ = −0.48 all tested separately within the same dataset, a Bonferroni or FDR correction should be applied.

(b) The χ² test for hexamer phase bias (χ² = 324.4, P < 0.001) has 5 degrees of freedom (6 phases), and the manuscript reports this is based on n = 218 observations. The expected count per cell under uniformity is ~36, which is adequate, but the test assumes independence of observations. If telotrons from the same genomic locus or the same gene tend to share phase, the effective n is lower and the P-value is inflated.

(c) Effect sizes are inconsistently reported. Some analyses give ρ, some give fold-enrichment, some give odds ratios, and some give only P-values. Nature Genetics strongly prefers consistent effect-size reporting throughout.

**Recommendation:** Apply FDR correction to all within-Haptista repair gene tests. Verify independence assumption for the hexamer phase χ² test. Standardize effect-size reporting.

---

## 8. THE "FEED-FORWARD" MODEL: Interesting but unsupported

**Severity: Medium**

The G4 feed-forward mechanism (G4 → DSB → telomerase → TTAGGG → more G4 → more DSB) is an appealing model, but the manuscript provides only correlative evidence: G4 enrichment at telotron flanks and telotron clustering in G4-dense regions. The model predicts that telotron density should increase over time in G4-rich regions, but the "burden–youth" correlation actually shows the opposite direction (most-affected genomes have the youngest telotrons), which could equally reflect recent activation rather than progressive accumulation via feed-forward.

**Recommendation:** Label the feed-forward mechanism explicitly as speculative/hypothetical. The correlative G4 data support preferential insertion at G4-rich loci but do not demonstrate the feed-forward loop.

---

## 9. SINGLE-AUTHOR MANUSCRIPT

**Severity: Low but notable**

Single-author papers are unusual at Nature Genetics (though not unprecedented). A reviewer may be concerned about:

(a) Independent verification of computational results. All analyses were performed by one person. Some journals request independent code review or reproduction by a second analyst for papers with complex computational pipelines.

(b) The Acknowledgments thank "members of the Carvunis laboratory for discussions" — if this reflects genuine intellectual contribution, consider whether co-authorship is appropriate. If not, the acknowledgment should be more specific about what was discussed.

**Recommendation:** Consider adding a collaborator who can independently verify key analyses. At minimum, ensure all code is well-documented and reproducible from raw data.

---

## 10. MISSING EXPERIMENTAL VALIDATION

**Severity: Medium — expected criticism but not fatal for a computational paper**

The manuscript acknowledges this: "Validation in cultured haptophyte and cryptophyte genomes, where telomerase knockdown should suppress new telotron formation, will be an important next step." This is honest and appropriate. However, a Nature Genetics reviewer may additionally ask:

(a) Can you detect telotron insertion in real-time using long-read sequencing of a telotron-rich organism in culture?

(b) Are there any published haptophyte or Phaeocystis genomes assembled from long reads that could confirm telotron structure without MAG assembly uncertainty?

(c) The Phaeocystis globosa genome was published (Liu et al., *Molecular Plant* 2022), and Emiliania huxleyi (Read et al., *Nature* 2013) is a well-studied haptophyte. Do these reference genomes show telotrons?

**Recommendation:** Search for Phaeocystis or related haptophyte reference genomes and check whether telotrons are present. If found, this would substantially strengthen the manuscript. If absent, explain why (e.g., different species, cultured vs. wild populations).

---

## 11. SPECIFIC TEXT ISSUES

(a) The abstract uses "recently" ("has recently generated tens of thousands of new introns") — on what timescale? The burden–youth correlation suggests recent activity but "recently" is vague for an abstract.

(b) "Telotrons demonstrate that an endogenous enzyme can drive intron gain at a scale rivalling transposable elements" — this is a strong concluding claim. The comparison with P. glacialis Introner-derived introns (~12,000) is apt, but the 844 high-purity telotrons vs. ~12,000 Introner introns suggests telotrons are actually an order of magnitude less abundant. The comparison is with the ~32,000 telo-containing introns (including degraded ones), which is less clean.

(c) The strand configuration analysis is elegant but very long for a Brief Communication. This level of detail belongs in an Article format.

(d) "Depletion of Ku70 in affected lineages implicates shifted break-repair balance as a predisposing factor" in the abstract — the abstract should mention Lig4 and Rad51 as well, since the updated analysis includes all three. Or, since the abstract is already too long, simplify to "depletion of NHEJ and HR repair genes."

---

## 12. KEY LITERATURE TO ADD

The manuscript currently has ~21 references. These should be added:

- **Nergadze et al., *Genome Research* (2004):** ITS arise from DSB repair in primates — direct mechanistic precedent
- **Nergadze et al., *Genome Biology* (2007):** Telomerase RNA retrotranscription in mammalian DSB repair
- **Rijal et al., *PLOS ONE* (2025):** Ku70 SAP domain absent in Haptophyta — independent support for repair gene depletion
- **Saary et al. / EukCC, *Genome Biology* (2020):** For citing eukaryotic MAG quality assessment methodology

---

## OVERALL ASSESSMENT

This is a genuinely exciting discovery with strong mechanistic evidence. The telomerase fingerprints (template boundary signatures, strand bias, hexamer phase bias) are compelling and difficult to explain by any alternative mechanism. The breadth of evidence — splicing validation, ortholog comparison, lineage-specific repeat units, G4 enrichment, repair gene depletion — builds a convincing cumulative case.

**However, the manuscript needs:**

1. **Format decision:** Submit as an Article (recommended) or drastically cut for Brief Communication
2. **Critical citations:** Nergadze 2004/2007 and Rijal 2025 must be added
3. **Assembly quality data:** BUSCO/EukCC scores for all telotron-bearing MAGs
4. **Stronger microsatellite exclusion:** Dedicated analysis, not a single sentence
5. **Tempered causal language** around repair gene associations
6. **Reference genome check:** Search for telotrons in published haptophyte genomes
7. **Multiple testing correction** for repair gene analyses

**Verdict:** With these revisions, this is a strong candidate for Nature Genetics as an Article. The core finding — telomerase creating spliceosomal introns — is novel and significant. The mechanistic evidence is substantially stronger than a typical correlative genomics paper. The main risk is that a reviewer dismisses the finding as an assembly artifact, so the assembly validation needs to be bulletproof.
