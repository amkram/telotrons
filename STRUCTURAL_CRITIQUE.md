# Structural Critique — Telotron Manuscript
## Weaknesses, Unsupported Claims, Incorrect Inferences

Generated: 2026-03-28

---

## A. CRITICAL STRUCTURAL WEAKNESSES

### 1. Near-total dependence on a single low-completeness MAG

Almost every quantitative analysis in the paper — strand orientation, phase enrichment, splice sites, metatranscriptomic validation, orthologs, G4 enrichment — comes from TARA_PSW_86_MAG_00284. This MAG has **BUSCO completeness of only 29.0%**. The top 10 telotron-bearing MAGs average BUSCO ~32%. A reviewer will immediately question whether patterns observed in a 29%-complete MAG are representative or assembly artifacts.

**The paper already addresses this** with consistency across 8 independent MAGs and assembly validation tests, but these controls are buried in parentheticals. Given the low BUSCO, these controls deserve more prominent treatment.

**FIX**: Add a sentence in para 1 or para 3 acknowledging the low BUSCO and explicitly noting the multi-MAG consistency check as partial mitigation. Consider stating the primary MAG BUSCO value.

### 2. "30 MAGs across four supergroups" — misleading framing

The per_mag_telotron_statistics.csv shows 257 MAGs total, spanning **7 supergroups** (Archaeplastida, Haptista, Opisthokonta, Stramenopiles, Alveolata, Cryptista, Rhizaria). The "30 MAGs" claim presumably refers to MAGs with HIGH-PURITY telotrons specifically (≥85% purity), not just any telomeric introns. This filtering criterion should be stated explicitly. Similarly, "four supergroups" should specify which four.

**FIX**: Change "30 MAGs across four supergroups" to specify the HP-telotron threshold and name the supergroups.

### 3. "691 tested" for introner exclusion — unexplained count

The main text says "0.1% of 691 tested" for terminal inverted repeats, and the Methods says "691 high-purity telotrons." But HP telotrons number 422 unique (844 doubled). 691 is neither 422 nor 844. This count is unexplained and likely comes from a different analysis pipeline or filtering. **This will confuse reviewers.**

**FIX**: Either verify and explain the 691 count, or update to match the 422 unique HP count.

---

## B. UNSUPPORTED OR WEAKLY SUPPORTED CLAIMS

### 4. "CAGG as the top motif in the stramenopile MAG analysed"

Based on a SINGLE stramenopile MAG (TARA_MED_95_MAG_00464). The boundary data for this MAG appears to have very low total array count. A single-MAG observation with potentially few arrays is weak evidence for a "lineage-specific" pattern. CTAA in haptophytes is strong (n=1,812 arrays, 74×); CAGG in stramenopiles may be noise.

**FIX**: Either (a) provide the n for stramenopile arrays and the enrichment fold, or (b) remove the CAGG claim and instead say "with other boundary motifs in non-haptophyte lineages (data not shown)."

### 5. "40–57% near-telomeric hexamers... indicating degradation"

Near-telomeric hexamers (≤1 substitution from canonical) comprise **5.6% of all possible hexamers** (228/4096). The 13.2% rate in non-telo introns is already 2.4× the random expectation, suggesting either compositional bias in intron sequences or some residual telomeric signal in the "non-telo" class. The 40-57% in telo-containing is clearly elevated, but the argument would be stronger with the proper null expectation stated.

**FIX**: Add the 5.6% random expectation to contextualize the 13.2% non-telo background.

### 6. Phase analysis n=248 out of 422 — what about the missing 174?

Only 248 of 422 HP telotrons had identifiable 5' boundary phase (a canonical hexamer at the first position after GT). The remaining 174 (41.2%) did not match any canonical rotation. This is never mentioned. If phase is a diagnostic "catalytic fingerprint," why does it fail for 41% of HP telotrons?

Possible explanations: boundary degeneration, non-standard telomerase initiation, or mixed-strand arrays. But the paper should acknowledge this.

**FIX**: Add a clause noting that phase was determinable for 248 of 422 HP telotrons, with the remainder carrying degenerate boundary sequences.

### 7. "Phaeocystis globosa" nuclear genome not published

The manuscript cites ref 13 for P. globosa genome, but the only P. globosa genome paper found on PubMed (Song et al. 2021) describes the **mitochondrial** genome, not the nuclear genome. If no nuclear genome assembly exists for P. globosa, this claim is unverifiable. E. huxleyi (Read et al. 2013 Nature) does have a nuclear genome.

**FIX**: Verify ref 13 actually describes a nuclear genome. If not, remove P. globosa from this claim.

### 8. "3 bp of template matching—the minimum for stable annealing"

The claim that 3 bp is the minimum for stable annealing of telomerase at a DSB is presented without citation. Is this from telomerase biochemistry literature? The in vitro studies of telomerase template annealing typically show the template alignment region is longer. This specific number needs a citation or should be qualified.

**FIX**: Add a citation for the 3bp minimum, or soften to "a short microhomology."

---

## C. INCORRECT OR OVERREACHING INFERENCES

### 9. Strand bias as "selection filter" — conflates mechanism with selection

The 11.5:1 template:coding ratio (172:15) is framed as "a selection filter favouring the splice-competent class." But this analysis EXCLUDES the 231 converging introns (54.7% of all HP telotrons). Converging introns have both TTAGGG and CCCTAA — they ARE the predicted DSB product (telomerase extending both overhangs). The single-strand classes (template-only or coding-only) represent cases where only one overhang was extended, which could reflect mechanistic bias (e.g., asymmetric resection) rather than selective removal of coding-strand introns.

The 96.5% vs 40.0% GT-AG difference IS strong evidence for the splice-site prediction. But the 11.5:1 count ratio could be partly mechanistic.

**FIX**: Acknowledge that the strand ratio could reflect both mechanistic bias and purifying selection. The splice-site quality difference (96.5% vs 40%) is the stronger evidence.

### 10. "Self-reinforcing system" — unsupported positive feedback

The final paragraph claims "once repair erodes past a threshold, telomerase captures an increasing fraction of DSBs, and the resulting introns become permanent genomic features." This implies positive feedback: repair erosion → more telotrons → more repair erosion. But there's no evidence that telotron insertion impairs repair capacity. The repair gene depletion could be independent of telotron accumulation (e.g., driven by drift in small populations). The correlation between Ku70/Lig4/Rad51 depletion and telotron burden could be a one-way causal arrow, not a feedback loop.

**FIX**: Remove the "self-reinforcing" framing, or explicitly state this is speculative and note the lack of evidence for the feedback direction.

### 11. ITS in yeast extrapolation

The paper cites Aksenova et al. 2013 to argue that "unspliced telomeric insertions in coding regions are highly deleterious." But Aksenova placed ITS **within an existing intron** of a reporter gene. The rearrangements arose from replication fork stalling at the telomeric repeats, not from failure to splice. The manuscript's inference — that the yeast data demonstrates the lethality of unspliced telotrons specifically — overstates the parallel.

**FIX**: Reframe: "interstitial telomeric sequences engineered into yeast chromosomes cause replication-dependent gross chromosomal rearrangements²⁰, indicating that telomeric tracts within genes are inherently destabilizing."

### 12. "Unlike all previously described intron-gain mechanisms"

This claim that telotrons are the first intron-gain mechanism requiring only an endogenous enzyme is strong, but could a reviewer argue that reverse-splicing intron gain (Roy & Irimia 2009) also requires only endogenous machinery (the spliceosome acting in reverse)? The distinction is that telotrons don't require an RNA intermediate, but the manuscript doesn't make this distinction explicitly.

**FIX**: Clarify: "unlike all previously described intron-gain mechanisms, which require either an exogenous transposable element or a spliced RNA intermediate."

---

## D. MISSED OPPORTUNITIES FOR STRENGTHENING

### 13. Metatranscriptomic validation is buried

The MGT validation (213/422 = 50.5%, Fisher P = 4.4×10⁻⁸) is arguably the strongest evidence in the paper — direct physical proof that telotrons are actively spliced. It's currently in paragraph 3, sandwiched between the ortholog comparison and the microsatellite exclusion. This deserves more prominence. The fact that HP telotrons have a HIGHER confirmation rate than non-telo introns (50.5% > 37.0%) is remarkable and underutilized.

### 14. Converging introns as strongest DSB evidence

The 231 converging introns (54.7%) — with head-to-head TTAGGG and CCCTAA arrays — are the most direct prediction of the DSB model. They're barely mentioned in the main text. A converging intron is the unmistakable signature of telomerase extending BOTH 3' overhangs at a single DSB. This is a powerful observation that could be highlighted.

### 15. Abstract word count

The abstract appears to exceed 70 words (Nature Genetics Brief Communication limit). Current abstract is ~68 words — verify and trim if needed.
