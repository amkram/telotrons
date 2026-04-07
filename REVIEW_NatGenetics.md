# Peer Review: Nature Genetics

**Manuscript:** Widespread de novo intron creation by telomerase-mediated double-strand break repair in marine eukaryotes

**Author:** Alexander M. Kramer

**Reviewer assessment:** This manuscript presents a striking and well-supported claim — that telomerase acting at DNA double-strand breaks creates de novo spliceosomal introns (termed "telotrons") across diverse marine eukaryotic lineages. The central finding, if validated, would represent a fundamentally new mechanism of intron gain, orthogonal to Introner elements, and would be of broad interest to the genomics, molecular evolution, and chromatin biology communities. The data are extensive, the statistical analyses are generally appropriate, and the multiple independent lines of evidence (telomerase boundary signatures, strand asymmetry, lineage-specific repeat types, ortholog comparisons, repair gene depletion) are individually compelling and mutually reinforcing.

I recommend this manuscript for publication in Nature Genetics after revision to address formatting constraints, several evidentiary gaps, and a number of points where claims outrun the data.

---

## 1. MAJOR ISSUES

### 1.1 Abstract exceeds Nature Genetics limit (Critical)

The abstract is 241 words; Nature Genetics requires ≤150 words. This must be cut by ~90 words. The abstract currently reads more like a mini-review than a structured summary. Suggestion: remove the mechanistic detail about GT/AG dinucleotides within the repeat (this belongs in the main text, where it already appears), compress the cryptophyte/arthropod comparison, and tighten the opening framing.

### 1.2 Extended Data figures far exceed allowance (Critical)

There are 23 Extended Data figures. Nature Genetics typically allows ~10. Many of the current ED figures appear to cover overlapping aspects of the G4 cascade analysis (ED Figs 7, 10, 11, 13 all address cascade/G4 relationships). Recommendations:

- Consolidate ED Figs 7, 10, 11, 13 into 1–2 multi-panel figures
- ED Figs 18, 19, 20, 21 (purity spectrum, architecture, multi-MAG comparison, cutoff determination) could be merged into 2 figures
- Consider moving some panels to Supplementary Information (which has no figure limit at Nature Genetics)
- The new ED Fig 23 (gapped alignments) is strong and should be retained — it provides the most direct evidence for de novo insertion

### 1.3 Single-author MAG-based study: reproducibility and assembly concerns

The entire study rests on metagenome-assembled genomes, which are known to be susceptible to chimeric contigs, contamination, and mis-assembly — particularly in repetitive regions. While the manuscript addresses this (paragraph 6, ED Fig 8), the controls could be strengthened:

- **Missing control:** A synthetic benchmark where known telomeric arrays of varying lengths are spiked into simulated metagenomic reads, assembled, and tested for recovery fidelity. This would directly address the concern that assemblers might preferentially collapse or expand tandem repeats.
- **Missing control:** Independent long-read validation. Even a single MAG re-sequenced with PacBio HiFi or Oxford Nanopore would substantially strengthen the assembly validation. The manuscript acknowledges this implicitly by suggesting future work in cultured organisms, but a reviewer at Nature Genetics will likely request at least one long-read dataset.
- The statement that telotrons are "uniformly distributed along contigs (not boundary-clustered)" is reassuring but not shown quantitatively. Consider adding a formal test (e.g., comparing observed telotron distance-from-contig-end distribution to uniform expectation).

### 1.4 Ortholog evidence is limited to a single MAG's gene set

The new paragraph on ortholog confirmation (paragraph 10) is an important addition, but 52 orthologs out of 844 strict telotrons (6.2%) is a modest confirmation rate. This limitation should be stated explicitly. The low rate likely reflects the difficulty of finding intron-free orthologs in related MAGs rather than a problem with the telotrons themselves, but this should be discussed. Can the authors quantify how many telotrons have orthologs that *also* contain the intron (shared telotrons), which would indicate insertion before speciation?

### 1.5 Target site duplications need mechanistic interpretation

The discovery that 71% of ortholog-confirmed telotrons carry TSDs (3–21 bp) is significant and new. However, the manuscript does not discuss the mechanistic implications adequately:

- TSDs are typically associated with staggered-cut mechanisms (e.g., transposase or endonuclease). What enzyme would create the staggered cut at a DSB? If the break is blunt (as most DSBs are), TSDs should not arise. If they do, this suggests either: (a) the breaks are processed by exonucleases creating staggered ends before telomerase extension, (b) fill-in synthesis after telomerase creates the duplication, or (c) the insertion mechanism involves a target-primed reverse transcription-like step.
- The TSD size distribution (3–21 bp, mean 6 bp) is unusually variable compared to transposon-generated TSDs (which are typically fixed-length for a given element). This variability needs explanation.
- The 29% of telotrons *without* TSDs also needs discussion — are these consistent with blunt-ended DSB repair?

### 1.6 The Methods section is too brief

At 420 words, the Methods section is underdeveloped for a study of this complexity. Key missing details:

- **Intron annotation:** How were introns called? MetaEuk gene models are mentioned, but what parameters? How were splice sites validated beyond the GT-AG census?
- **Ortholog identification:** How were orthologs identified across MAGs? What alignment tool, identity threshold, and synteny requirements were used? The `homolog_results.pkl` pipeline is not described.
- **TSD detection:** The algorithm for finding TSDs should be described (search window, minimum length, handling of ambiguity).
- **Statistical tests:** Multiple comparisons are reported without correction. The G4 analyses involve many tests across overlapping subsets — what FDR control was applied?
- **Code availability:** "Code is available at [repository URL upon publication]" — reviewers will want to see the code during review. Consider providing a temporary repository link.

---

## 2. MODERATE ISSUES

### 2.1 Telomerase identification is indirect

The manuscript claims telomerase origin based on four fingerprints: (1) template-boundary enrichment, (2) strand asymmetry, (3) lineage-specific repeat type, and (4) repeat purity independent of length. While collectively compelling, no direct evidence of telomerase protein (TERT) activity is provided. The manuscript mentions searching for TERT domains (Methods), but does not report whether TERT was found in the primary MAG or other affected MAGs. If TERT is absent, this would be a significant problem; if present, this should be stated prominently.

### 2.2 The Ku70/Pif1 correlation is observational

The anti-correlation between Ku70/Pif1 density and telotron burden (Fig. 2a,b) is presented as evidence for a mechanistic link, but this is purely correlational. The manuscript acknowledges that repair genes are not preferentially disrupted by telotrons, but alternative explanations should be discussed more thoroughly:

- MAG completeness could confound both measurements (less complete MAGs have fewer genes AND possibly more repetitive content)
- The negative correlation with total gene count (ρ = −0.79) is itself concerning — it suggests systematic differences in MAG quality
- Protein domain detection by motif searching has inherent false-negative rates that could vary with MAG quality/completeness

### 2.3 The "intrinsic splice competence" argument needs qualification

The elegant observation that TTAGGG contains GT at repeat junctions and AG internally is central to the narrative. However:

- The manuscript states "89% of telotrons use canonical GT donors and 95% use canonical AG acceptors" (paragraph 12), while the corrected census shows 93.5% GT-AG. These numbers should be reconciled — the 89%/95% may be from an earlier analysis across all MAGs. State clearly which population each number refers to.
- The GT and AG dinucleotides are necessary but not sufficient for splicing. Branch point sequences, polypyrimidine tracts, and ESE/ESS elements all matter. The manuscript mentions "non-telomeric leader and trailer sequences" providing splice context, but does not analyse these computationally. A branch-point prediction analysis on telotron 3' boundaries would strengthen this section.

### 2.4 Metatranscriptomic validation: effect size vs. statistical significance

The splicing validation (46.7% vs 41.8%, P = 3.0 × 10⁻¹³) has high statistical significance but modest effect size. The 4.9 percentage-point difference could reflect many things — the key question is whether telotrons splice at all, which they clearly do. But the framing ("telotrons splice at 46.7% vs 41.8%") could be misread as suggesting telotrons splice *better* than normal introns, which is probably not the intended claim. Clarify that the comparison demonstrates comparable splicing efficiency, not superiority.

### 2.5 Species count uncertainty is large

"32–49 distinct species-level lineages" is a wide range (53% uncertainty). For a Nature Genetics paper, the taxonomic framework should be tightened. Consider using ANI (average nucleotide identity) clustering at 95% to define species-level OTUs, which would give a more precise count.

---

## 3. MINOR ISSUES

### 3.1 Reference list is sparse

19 references for a paper of this scope is unusually low. Missing citations that reviewers will expect:

- Kramer & Carvunis (2018) on de novo gene birth and intron evolution — relevant context
- Huff et al. (2016) on telomere-to-telomere assembly and ITS characterization
- Nandakumar & Cech (2013) on telomerase mechanism — more authoritative than the Tetrahymena structure alone
- Myler et al. (2017) on DSB end processing and the competition between repair pathways
- Any reference on ITS biology beyond Bolzán (2017), which focuses on vertebrates

### 3.2 Terminology

- "Strict telotron" vs "very strict telotron" vs "telo-containing intron" — three tiers plus "ambiguous" and "non-telo" makes five categories. Consider simplifying the naming or providing a clear table/schematic in the main text.
- "Telotron" is a good coinage, but it could be confused with "telomeron" (a different concept). Briefly note this in the introduction.

### 3.3 Figures

- **Fig. 1** has 5 panels but only 2 main figures. Consider whether Fig. 1e (boundary sharpness by age) belongs in ED — it's interesting but tangential to the core claim.
- **Fig. 2** panel descriptions are minimal. The figure legend for Fig. 2d ("Exonic ITS are shorter than telotrons and intergenic ITS") is a conclusion, not a description.
- **ED Fig 23** (gapped alignments) is excellent and could arguably be promoted to a main figure panel — it provides the most intuitive visual evidence of de novo insertion.

### 3.4 Narrative structure

The manuscript flows well but is front-loaded with characterization (paragraphs 5–8) before reaching the mechanistic model (paragraph 9) and the direct ortholog evidence (paragraph 10). Consider moving the ortholog comparison earlier — perhaps right after the splicing validation — since it provides the most accessible proof of de novo insertion for a general reader.

### 3.5 Small numerical issues

- Paragraph 5: "844 strict de novo telotrons (698 strict, 146 very strict)" — 698 + 146 = 844, check. But then "2,371 strict de novo telotrons" across 8 MAGs, "of which 248 are very strict" — is 248 across all MAGs while 146 is in the primary MAG alone? Clarify.
- Paragraph 9: "Among 462 converging strict and very strict telotrons" but earlier "844 strict de novo telotrons" — the remainder are presumably single-strand. State this explicitly.
- Paragraph 12: "89% of telotrons use canonical GT donors" conflicts with the corrected value of 93.5% GT-AG in the primary MAG. If 89% is across all MAGs, say so.

### 3.6 Writing

The prose is generally strong and clear. A few suggestions:

- The opening sentence of the abstract ("Intron gain mechanisms remain poorly understood despite evidence for episodic bursts throughout eukaryotic evolution") is somewhat generic. Consider leading with the finding.
- "Orthogonal to Introner elements" (last sentence of abstract) — "orthogonal" is imprecise here. "Mechanistically distinct from" is clearer.
- Avoid "near-perfect" (paragraph 9) — either quantify or say "complete."

---

## 4. SUMMARY ASSESSMENT

**Novelty:** Exceptional. If the core claim holds, this is a genuinely new mechanism of intron creation with broad implications for genome evolution.

**Evidence quality:** Strong overall, with some gaps. The multi-pronged fingerprint evidence (boundary signatures, strand asymmetry, lineage-specific repeats) is the strongest part. The ortholog confirmation with TSDs is a valuable recent addition. The main weaknesses are the lack of long-read validation, incomplete methods, and the correlational nature of the repair gene analysis.

**Significance:** High. Telotrons would join Introner elements as only the second well-characterized mechanism of de novo intron gain at scale, and the first driven by an endogenous enzyme rather than a mobile element.

**Suitability for Nature Genetics:** Yes, after revision. The core finding is of sufficient novelty and significance for this journal. The main barriers to acceptance are formatting (abstract length, ED figure count), methodological detail, and the need for at least one orthogonal validation (ideally long-read sequencing of a single MAG).

### Recommended verdict: **Revise and resubmit (major revisions)**

**Priority revisions:**
1. Cut abstract to ≤150 words
2. Consolidate ED figures to ≤10 (move remainder to Supplementary)
3. Expand Methods to ≥1,500 words with full algorithmic detail
4. Add TERT presence/absence data for affected MAGs
5. Reconcile the 89%/93.5% GT-donor discrepancy
6. Discuss TSD mechanistic implications
7. Address assembly validation with at least one quantitative control
