# Nature Genetics Referee Report (v2)

**Manuscript:** Widespread de novo intron creation by telomerase-mediated double-strand break repair in marine eukaryotes

**Verdict: Revise and resubmit (minor–moderate revisions)**

This manuscript has improved substantially since the prior round. The abstract is tighter, the narrative is better ordered, the Methods are greatly expanded, and weak evidence has been trimmed. The central claim—that telomerase creates de novo introns at DSBs—remains extraordinary and is now supported by a more focused evidence chain. Below I identify remaining issues organized by severity.

---

## Major issues

### 1. The 50,000 number in the abstract is misleading without qualification

The abstract states "over 2,000 de novo intron gains mediated by telomerase ('telotrons'), plus 50,000 additional telomere-containing introns." The 50,000 includes introns with as little as 10% telomeric content—many of which may be pre-existing introns that acquired telomeric insertions rather than de novo telomerase-generated introns. The degenerate analysis (3–4× near-telo enrichment) supports a telomerase origin for many of these, but not all. The abstract should temper this number, e.g., "plus ~50,000 partially degraded telomere-containing introns representing a continuous decay series." Without qualification, a reader could interpret this as 52,000 confirmed telomerase-generated introns, which would overstate the evidence.

### 2. The metatranscriptomic comparison is statistically confusing

The text reports that telotrons splice at 46.7% vs 41.8% for non-telotron introns with P = 3.0 × 10⁻¹³, and then says this demonstrates telotrons are "efficiently spliced." But the significant P-value actually shows telotrons splice at a *higher* rate than controls—which is the opposite of what a skeptic would expect for a parasitic insertion. This result is never explained. Is the higher rate due to length bias (shorter introns splice more efficiently), expression bias (telotron-bearing genes are more expressed), or ascertainment bias (MGT reads preferentially capture shorter intron junctions)? This needs either an explanation or a reframing. Simply calling the rates "comparable" while reporting a P-value showing they are significantly different is internally contradictory.

### 3. Ortholog confirmation covers only 6% of strict telotrons

Only 52 of 844 strict telotrons (6.2%) have ortholog-confirmed de novo insertion. While the text previously explained this was removed in editing, some acknowledgment of this limitation is needed in the body text. Without it, a reviewer will flag it. Even one sentence ("The modest confirmation rate reflects the difficulty of identifying intron-free orthologs in incomplete MAG assemblies; most orthologous genes either share the telotron or lack sufficient coverage for alignment") would suffice.

### 4. Ku70/Pif1 correlations may reflect phylogenetic confounding

The anti-correlation between Ku70 density and telotron burden (ρ = −0.67) across 48 MAGs is reported, and the within-Haptista control (ρ = −0.68, n = 18) is good. However, the Pif1 correlation (ρ = −0.41, P = 4.3 × 10⁻³) is modest and the within-Haptista P-value for Pif1 is not reported. Is the Pif1 correlation significant within Haptista alone? If not, it may be entirely explained by phylogenetic structure (haptophytes have both high telotron loads and low Pif1 by lineage). This needs to be tested and reported.

---

## Moderate issues

### 5. Reference 21 (Myler et al.) is cited for transcription-coupled DSB repair but doesn't support the specific claim

Myler et al. (2017) describes Mre11-Rad50-Nbs1 initiation of break repair via single-molecule imaging. It does not specifically address transcription-coupled repair or template-strand 3′-end exposure during transcription. The claim that telomerase preferentially engages "the 3′ end exposed on the template strand during transcription-coupled DSB repair" is mechanistically reasonable but needs a more appropriate citation—perhaps Aguilera & García-Muse (2012) "Causes of genome instability" or similar work on transcription-associated DSBs.

### 6. The 844 vs 2,371 numbers create confusion

The primary MAG has 844 strict telotrons. The eight most affected MAGs have 2,371. These numbers appear in different paragraphs without clear context linking them. A reader encountering "844" in paragraph 3 and "2,371" later may wonder how they relate. Recommend stating "844 of 2,371 total across eight MAGs" or similar bridging language when the 844 is first introduced.

### 7. Strand configuration paragraph uses two different denominators

"462 of 844 high-purity telotrons (55%)" uses 844 as the denominator for converging fraction. Then "328 of 674 GT-AG telotrons (49%)" switches to 674 (GT-AG only) as the denominator for template-only fraction. This makes it impossible to compare the three categories against a single total. Present all three against a single denominator—either 844 (all strict) or 674 (GT-AG only)—to allow the reader to verify the percentages sum correctly.

### 8. "Intronization" model needs the CCCTAA problem stated explicitly

The ITS/intronization paragraph explains that TTAGGG arrays encode GT and AG dinucleotides. But CCCTAA (the reverse complement) does NOT contain a GT dinucleotide anywhere—meaning template-only telotrons, which are nearly half of all telotrons, cannot get their GT donor from the repeat itself. This is a fundamental asymmetry that should be stated explicitly: "In template-only telotrons (CCCTAA on the coding strand), the GT donor must be supplied by flanking genomic sequence or the non-telomeric leader, since CCCTAA lacks a GT dinucleotide in any reading frame." The current text's general claim about TTAGGG encoding splice sites obscures this important caveat.

### 9. No data availability or code availability section

The Methods state code will be at a GitHub URL "to be made public upon publication," but Nature Genetics requires a Data Availability statement describing where the raw data (MAG sequences, intron annotations, MGT reads) can be accessed. Since these are from the Tara Oceans catalogue, the appropriate accession numbers (BioProject, ENA, etc.) should be listed.

---

## Minor issues

### 10. Abstract is 152 words; Nature Genetics limit is 150

Two words over. Suggest cutting "de novo" from "de novo intron gains mediated by telomerase" (since "mediated by telomerase" already implies they are new), or "additional" from "50,000 additional."

### 11. Figure 1 legend lists 4 panels but the text references Fig. 1c and 1d

The figure legend describes panels (a)–(d) but the text only references 1a and 1b in the body. Panels 1c (telomeric base coverage) and 1d (template-strand bias) were referenced in the previous version but those text references were removed during trimming. Either restore the in-text references or move panels c and d to Extended Data.

### 12. The "controlling for expression level" parenthetical is imprecise

The metatranscriptomic comparison states it controls for expression by restricting to genes with both intron types. This is a valid design, but it controls for *gene identity*, not expression level per se. Genes with both types could still have varying expression. Rephrase to "controlling for gene identity" or explain the expression-matching procedure.

### 13. "32–49 distinct species-level lineages" is vague

The Methods explain this range comes from counting distinct classes (32) vs genera (49). This should be stated in the main text rather than requiring the reader to find the explanation in Methods. E.g., "representing 32 distinct classes and up to 49 genera."

### 14. Missing ED Figure reference: Fig. 13

The text references "Extended Data Fig. 13" for the within-gene G4 test but the ED figure numbering appears non-sequential (the referenced figures are 1, 5, 8, 9, 13, 15, 16, 17, 18, 22, 23). If the ED has been consolidated, renumber sequentially (ED Figs. 1–11). If it hasn't been consolidated yet, do so.

### 15. "Telotrons constitute a previously undescribed mechanism" — hedge

The discussion opens with "previously undescribed mechanism," which will be evaluated by the reviewers and editors. Consider whether any prior report of telomeric sequence in introns exists, even incidentally. Bolzán (ref 11) discusses ITS broadly. If no one has previously noted telomeric introns, the claim stands; if there are incidental observations, acknowledge them.

### 16. Ref 20 (Nandakumar & Cech) is cited in the Methods for TERT RT motif searching

This reference is about telomerase recruitment, not RT motif identification. Consider citing Lingner et al. (1997) "Reverse transcriptase motifs in the catalytic subunit of telomerase" (Science 276:561–567) or Nakamura et al. (1997) instead, which defined the TERT RT motifs.

### 17. The Tetrahymena TERT–TBE structure (ref 8) may not generalize

The boundary element structure is from Tetrahymena, which uses a different repeat (TTGGGG) than haptophytes (TTAGGG). The text should note that while the TBE mechanism is structurally conserved, the specific boundary position differs by species, making the haptophyte CTAA enrichment a prediction that matches expectation rather than a direct structural confirmation. Alternatively, cite vertebrate TERC boundary work if available.

---

## Summary

The manuscript presents a compelling discovery with strong computational evidence from multiple independent lines. The improvements since the last draft are significant—the abstract is sharper, the Methods are adequate, redundancy has been removed, and the strand-configuration taxonomy is now clearly presented. The major remaining issues are: (1) the 50,000 number needs qualification, (2) the metatranscriptomic P-value needs explanation, (3) the ortholog limitation needs acknowledgment, and (4) the Pif1 within-lineage control is missing. These are all addressable without new analyses. With these revisions, the manuscript would be appropriate for Nature Genetics.

**Recommendation: Minor–moderate revision.**
