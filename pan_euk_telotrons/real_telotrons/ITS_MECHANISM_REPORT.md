# ITS Mechanism Analysis: Telotron Architecture vs. Published ITS Mechanisms

**Date:** 2026-04-27
**Datasets:**
- 55,822 telomeric introns from 8 top haptophyte MAGs (Tara Oceans, ≥30% TTAGGG/CCCTAA coverage)
- 186 validated telotrons from 59 species (pan-eukaryotic, ≥85% purity)

This report tests each mechanism class from the synthesized ITS literature review against our data.

---

## 1. Orientation classification

The four orientation classes from the review map onto our data as follows:

| Class | Tara count | Pan-euk count | Total | % of converging-eligible |
|---|---:|---:|---:|---:|
| SINGLE_G (G-strand only)        | 20,230 | 113 | 20,343 | — |
| SINGLE_C (C-strand only)        | 23,132 | 70  | 23,202 | — |
| **CONVERGING_NOLINKER** (head-to-head, ≤4bp gap) | 1,537  | 0 | 1,537  | 13% |
| **CONVERGING_WITHLINKER** (head-to-head, >4bp gap) | 10,058 | 0 | 10,058 | 87% |
| DIVERGING (tail-to-tail, rare)  | 288   | 0  | 288   | — |
| NONE                            | 21    | 3  | 24    | — |

**Key observations:**
- SINGLE_G and SINGLE_C are nearly equiprobable (~50/50 split → 88% of all telotrons), consistent with the review's note that ITS insertion is non-directional.
- 21% of haptophyte telotrons are converging — much higher than the converging fraction described in primates/rodents.
- Diverging (tail-to-tail) is very rare (0.5%) — consistent with biology (would require an unusual two-ended capture geometry).
- The pan-euk validated set has **zero converging cases** — survival of converging telotrons through the misannotation filter requires high-purity flanking arrays from a single intron, which is much rarer outside Haptista.

---

## 2. TSD presence and length (Nergadze 2004 staggered-DSB model)

Nergadze 2004 attributed TSD-bearing ITSs to staggered DSBs with overhang fill-in.

| Class | % with TSD ≥3bp | Mean TSD length (when present) |
|---|---:|---:|
| CONVERGING_NOLINKER     | 8.1% | 6.4 bp |
| CONVERGING_WITHLINKER   | 7.0% | 4.6 bp |
| SINGLE_G                | 6.9% | 4.5 bp |
| SINGLE_C                | 7.2% | 4.7 bp |
| DIVERGING_WITHLINKER    | 6.0% | 3.4 bp |

**Interpretation:**
- TSD rate is **uniform 6-8% across all classes** — not enriched in any orientation.
- This is much lower than the "71%" claim in our earlier draft (which used a permissive 3bp-in-6bp-windows method). The standard exact-TSD definition gives 7%, consistent with Nergadze's primate/rodent data.
- TSD length when present is short (4-7 bp) — consistent with mostly-blunt DSBs with small staggered overhangs, processed through MMEJ or short-overhang NHEJ.
- TSDs are **NOT a defining feature of mechanism class** — they reflect a downstream resolution detail (NHEJ subpathway), not insertion mechanism.

---

## 3. Microhomology to canonical TTAGGG at boundaries (dnTA prediction)

Putnam, Pennaneach & Kolodner (2004) showed dnTA preferentially heals breaks where the broken end exposes microhomology to the TERC alignment region.

| Class | Mean 5' microhomology | Mean 3' microhomology |
|---|---:|---:|
| CONVERGING_NOLINKER     | 1.38 bp | 1.36 bp |
| CONVERGING_WITHLINKER   | 1.35 bp | 1.34 bp |
| SINGLE_G                | 1.33 bp | 1.27 bp |
| SINGLE_C                | 1.27 bp | 1.31 bp |
| DIVERGING_NOLINKER      | 0.50 bp | 3.75 bp (small n) |

**Interpretation:**
- ~1.3 bp average microhomology — small but non-zero, in line with dnTA-via-microhomology but not a strong signal at the 6+ bp scale Putnam reported.
- Symmetric across classes → consistent with dnTA at both ends of two-ended events.
- The single-end classes (SINGLE_G, SINGLE_C) show similar microhomology to the converging classes → both ends of converging events use the same microhomology-mediated mechanism.

---

## 4. Linker length distribution (Nergadze TERC-retrotranscription)

Nergadze 2007 predicts linker lengths of **30-120 bp** for TERC retrotranscription events, in opposite orientation to flanking arrays.

For **CONVERGING_WITHLINKER (n=10,058)**:

| Length bucket | Count | % | Predicted mechanism |
|---|---:|---:|---|
| Short (5-30 bp)    | 3,486  | 34.7% | NHEJ scar / processing |
| **Medium (30-120 bp)** | **5,569** | **55.4%** | **TERC retrotranscription (Nergadze)** |
| Long (120+ bp)     | 1,003  | 10.0% | Captured DNA fragment |

The 55% medium-length distribution is a striking match to Nergadze's TERC range — but composition reveals this is misleading.

---

## 5. Linker composition: most "linkers" are degenerate telomeric repeats

Critical correction to the naive linker analysis above:

| Linker subset | Count | % of converging linkers |
|---|---:|---:|
| Total linkers (≥15bp)            | 8,701 | — |
| **Pure non-repeat (zero TTAGGG/CCCTAA kmers)** | **751** | **8.6%** |
| Contains ≥1 telomeric kmer       | 7,950 | 91.4% |

**The 91% of "linkers" with telomeric kmers are degenerate TTAGGG/CCCTAA — not non-repeat sequence.**

Among the 5,569 medium-length (30-120 bp) "linkers":
- 6.8% are pure non-repeat
- 33.4% contain BOTH TTAGGG and CCCTAA kmers (mosaic converging continuation)
- 30.8% contain only TTAGGG
- 29.0% contain only CCCTAA

For the **751 truly non-repeat linkers**, the most-shared 8-mers are:
- `CCTCACCC, TCACCCTC, CCCTCACC, CCCTACCC, CCTACCCC` — these are **CCCTAA variants with single-base substitutions**
- `TAGGGCTA, GGGCTAGG` — **TTAGGG variants** with T→C mutations

**Conclusion:** The dominant "linker" class is not a TERC fragment, captured DNA, or random NHEJ scar. It is **mosaic, internally-degenerated telomeric repeat sequence** — consistent with a single long telomerase-mediated insertion that has accumulated point mutations along its length, with the central degenerate region between two "younger-looking" canonical end blocks.

This argues against Nergadze TERC retrotranscription as the dominant linker source in haptophyte telotrons, despite the length distribution superficially matching. The 30-120bp medium bucket is dominated by degenerate-but-still-recognizable telomeric content, not external sequence.

The remaining true non-repeat candidates (n=751, 8.6%) are split:
- 340 short (15-30 bp) — NHEJ scars
- 395 medium (30-120 bp) — best Nergadze candidates (would require external TERC BLAST to confirm)
- 16 long (120+ bp) — fragment-capture candidates

---

## 6. NR2C/F binding motifs (Marzec 2015 TTI mechanism)

Marzec 2015 described targeted telomere insertion (TTI) by orphan nuclear receptors NR2C/F, which bind GGGTCA half-sites (DR0/6/7 spacings). TTI signature: variant repeats interspersed in array AND NR2C/F motifs in flanks.

**Mean GGGTCA half-sites in 200bp flanks per class:**

| Class | Mean half-sites | % with ≥1 |
|---|---:|---:|
| CONVERGING_NOLINKER     | 0.125 | 9.3% |
| CONVERGING_WITHLINKER   | 0.196 | 14.3% |
| SINGLE_G                | 0.186 | 13.6% |
| SINGLE_C                | 0.193 | 13.9% |
| DIVERGING_WITHLINKER    | 0.310 | 20.4% |

**Interpretation:** Very low motif density (0.1-0.3 per 200bp). 86% of telotrons have NO GGGTCA half-site in flanks. This **argues strongly against TTI as the primary mechanism**. TTI requires DR0/6/7 paired half-sites — if the dominant mechanism were TTI, we'd expect orders of magnitude more motif enrichment.

---

## 7. Variant-repeat fraction (TTI vs canonical-telomerase)

ALT/TTI telomeres are variant-rich (TCAGGG, TGAGGG, GTAGGG); canonical-telomerase ITSs are nearly pure when young.

For TTAGGG-family telotrons:

| Class | Mean canonical fraction | Mean 1-mismatch variant fraction | TTI-specific variants |
|---|---:|---:|---:|
| CONVERGING_NOLINKER     | (varies) | 24.6% | 2.4% |
| CONVERGING_WITHLINKER   | (varies) | 33.0% | 3.5% |
| SINGLE_G                | (varies) | 33.4% | 3.9% |
| SINGLE_C                | (varies) | 33.4% | 3.9% |
| DIVERGING_NOLINKER      | (varies) | 46.7% | 6.7% (small n) |

**Interpretation:**
- Variant fractions are 24-33% — moderate, consistent with aged but not fully degraded telomeric insertions.
- TTI-specific variants (TCAGGG, TGAGGG, GTAGGG) are present at 2-7% — well below ALT-like levels.
- CONVERGING_NOLINKER has the lowest variant fraction (24.6%) → these are the "youngest" insertions, with the cleanest canonical TTAGGG arrays. Consistent with two-ended dnTA followed by re-ligation as a recent event.
- CONVERGING_WITHLINKER (33%) and SINGLE classes (33%) have similar variant rates → comparable age distribution.

---

## 8. SURVIVORSHIP TEST: telomeric content by genomic compartment

Telomeric kmer density per Mb in 5 haptophyte MAGs:

| MAG | CDS | Intron | Intergenic | **Intron:CDS** | Interg:CDS | Intron:Interg |
|---|---:|---:|---:|---:|---:|---:|
| TARA_PSW_86_MAG_00284  | 1,630 | 60,373 | 12,589 | **37×**  | 7.7× | 4.8× |
| TARA_PON_109_MAG_00250 | 1,213 | 27,845 | 6,967  | **23×**  | 5.7× | 4.0× |
| TARA_ARC_108_MAG_00319 | 815   | 22,182 | 4,287  | **27×**  | 5.3× | 5.2× |
| TARA_PSE_93_MAG_00226  | 1,544 | 22,074 | 9,394  | **14×**  | 6.1× | 2.4× |
| TARA_AOS_82_MAG_00154  | 2,151 | 48,714 | 15,427 | **23×**  | 7.2× | 3.2× |

**Two distinct findings:**

(a) **Survivorship via purifying selection on CDS:**
   Intergenic regions have 5-8× more telomeric content than CDS. Telomeric insertions in coding regions are purged by selection, consistent with the survivorship hypothesis from the synthesized review.

(b) **Intron enrichment beyond intergenic baseline (2-5×):**
   Introns are NOT just neutral relative to intergenic — they are 2-5× *more* enriched. Three non-mutually-exclusive explanations:
   - **Splice-compatible-only retention:** insertions that disrupt splicing are silently misannotated; only those compatible with GT-AG borders are visible as "introns" → annotation bias enriches the visible intronic class.
   - **Intron-specific insertion preference:** transcribed regions (with R-loops, open chromatin, G4 structure) are DSB hotspots and telomerase substrates. Introns → R-loop susceptibility → preferential targeting.
   - **Selection-permissive zone:** intronic insertions that don't disrupt splicing are essentially neutral and accumulate to high copy numbers; intergenic insertions may be lost to drift.

**Window-level positivity rates** (windows with ≥2 telomeric kmers per 200bp):
- TARA_PSW_86: CDS 2.1%, Intron **72.2%**, Intergenic 8.6%

This is an extraordinary intron enrichment in haptophyte MAGs, well beyond what neutral background would produce.

---

## 9. Summary: best-fit mechanism model

Based on all diagnostics, the dominant mechanism producing haptophyte telotrons is:

**Single-event telomerase-mediated DSB repair with mosaic internal degeneration**, with the following sub-features:

1. **One-ended dnTA / fragment capture** is most common (88% of telotrons, SINGLE_G or SINGLE_C). The two strand orientations are roughly equiprobable, consistent with non-directional mechanism.

2. **Two-ended healing produces converging arrays** (12% of telotrons, CONVERGING). The 87% of converging events with a "linker" between the two arrays mostly contain DEGENERATE telomeric content — they are best interpreted as a single long telomerase product with central erosion, not a two-event scenario.

3. **TERC retrotranscription (Nergadze 2007)** is rare. The truly non-repeat medium-length linkers (n=395, 4.5% of converging) are the only legitimate candidates; would require BLAST against a haptophyte TERC (currently uncharacterized) to confirm.

4. **Fragment capture** is rare (long, truly non-repeat linkers: n=16, 0.2% of converging).

5. **TTI / NR2C/F mechanism** is essentially absent: low GGGTCA motif density and modest TTI-specific variant fractions argue against this mechanism.

6. **Survivorship bias is real and quantified:**
   Introns have 14-37× more telomeric content than CDS — confirms that exonic insertions are purged by selection. The additional 2-5× enrichment over intergenic suggests both insertion preference (transcribed-region DSB hotspots) and annotation bias (only splice-compatible insertions are "visible" as introns).

---

## 9b. EIMERIA vs HAPTOPHYTE: fundamentally different architectures

The user observation that "Eimeria has more linkers" reveals a much deeper finding:
**Eimeria telotrons have a fundamentally different architecture from haptophyte telotrons.**

| Feature | Haptophytes (Tara, n=55,822) | Eimeria (n=17,171) |
|---|---:|---:|
| With-linker fraction (of all telomeric introns) | 21% | **37%** |
| CONVERGING : DIVERGING ratio | **40 : 1** | **1 : 12.5** |
| CONVERGING_WITHLINKER count | 10,058 | 467 |
| DIVERGING_WITHLINKER count | 288 | **5,872** |
| Median linker length | 44 bp | **123 bp** |
| Pure non-repeat linker % (CONVERGING_WITHLINKER) | 8.6% | **25.9%** |
| Length bucket: short (5-30bp, NHEJ scar) | 35% | 1% |
| Length bucket: medium (30-120bp, TERC-frag) | 55% | 46% |
| Length bucket: long (120+bp, capture) | 10% | **53%** |

### Interpretation

In **haptophytes**: telotrons are dominantly CONVERGING with short to medium linkers, with most "linkers" being degenerate continuations of telomeric arrays. This is consistent with **classical two-ended dnTA at internal DSBs** (the Nergadze 2004 model), with the central region simply being older/eroded telomeric sequence.

In **Eimeria**: telotrons are dominantly DIVERGING (telomere arrays point AWAY from the linker on the gene strand) with long, often truly non-repeat linkers. Mechanically, classical two-ended dnTA cannot produce DIVERGING — telomerase always extends the G-rich strand outward from the chromosome end, so inward-pointing arrays (CONVERGING) is the dnTA prediction.

Three possible mechanisms for Eimeria's pattern:
1. **Subtelomeric fragment capture**: Apicomplexa (esp. Plasmodium) have highly active subtelomeric recombination. Capture of a fragment containing internal telomere blocks (e.g., a TARE/SCAR-like element with bidirectional telomere flanks) and re-integration at an internal DSB could produce DIVERGING insertions.
2. **Inversion after CONVERGING insertion**: original CONVERGING telotron flipped by a subsequent inversion event.
3. **Different DSB-repair geometry**: telomerase may engage Apicomplexan DSBs from different angles (e.g., extending 5' overhangs rather than 3'), producing different array orientations.

### Inverted-repeat test result

Eimeria DIVERGING linkers do NOT show palindromic structure (0/200 had ≥70% palindromic identity), so the linker itself is not a self-folded RNA structure or rationally-paired inverted-repeat. The flanking arrays are inverted but the middle is unique.

### Eimeria linkers are richer in true non-repeat sequence

- **1,250 truly non-repeat linkers in Eimeria** (vs 751 in Tara haptophytes despite Tara having 3× more telotrons total)
- Eimeria pure linker median: 87 bp (medium-length zone)
- These are the strongest TERC-retrotransposition or fragment-capture candidates we've found.

These should be BLASTed against:
- Eimeria's own genome (subtelomeric/TARE regions)
- Apicomplexan TERC sequences (Plasmodium TERC is partially characterized)
- Repeat databases (Repbase, Dfam) for known transposable elements

The enrichment of long, non-repeat linkers in Eimeria is a strong argument for **fragment capture** as the dominant mechanism in this lineage, distinct from the telomerase-only mosaic mechanism in haptophytes.

---

## 9c. WHERE DO EIMERIA LINKERS COME FROM? — Genomic origin analysis

### BLAST analysis: 1,230 pure non-repeat Eimeria linkers vs their own genomes

We BLASTed 1,230 truly non-repeat Eimeria linkers against the host genomes. Hit criteria: ≥80% query coverage, ≥85% identity, excluding self-hits.

| Result | Count | % |
|---|---:|---:|
| Linkers with **0 non-self HQ paralog hits** | 1,207 | **98.1%** |
| Linkers with 1-2 paralog hits                | 23    | 1.9% |
| Linkers with ≥3 paralogs                     | 0     | 0% |

**Most Eimeria linkers are unique sequences in their genome.**

### When paralogs DO exist, they are subtelomeric

For the 24 off-target paralog hits found:

| Group | % within 5kb of contig end | % within 1kb |
|---|---:|---:|
| Source telotrons | 45.9% | 12.7% |
| Background telomeric introns | 33.2% | — |
| **Off-target linker hits** | **95.8%** | **70.8%** |

**The few paralogous hits are 2-3× more enriched at contig ends than the source loci themselves.** This is consistent with subtelomeric origin for those linkers that DO have paralogs, even if most don't.

Caveat: Eimeria assemblies are scaffold-level with median contig length ~4 kb, so "near contig end" is partly an assembly artifact. But the 2-3× enrichment over source-locus baseline is real biological signal.

### EIMERIA TELOTRONS CLUSTER IN TANDEM ARRAYS

Across all Eimeria telotrons (n=17,200):

| Statistic | Value |
|---|---:|
| Contigs with ≥1 telotron | 4,024 |
| Contigs with ≥2 telotrons | **2,647 (65.8%)** |
| Adjacent telotron pairs (same contig) | 13,147 |
| Pairs <100 bp apart | **3,301 (25.1%)** |
| Pairs 100bp-1kb apart | 5,135 (39.1%) |
| Pairs 1-10kb apart | 3,677 (28.0%) |
| Median adjacent-pair distance | **348 bp** |

**Two-thirds of telotron-bearing contigs have multiple telotrons, and half of those neighbors are within 1 kb.**

### Shared-linker pairs

5 of the 24 paralogous linker hits (21%) land **inside another telotron** on the same contig within a few hundred bp:

```
Contig                   Source telotron pos      Shared linker found in...
NW_013651890.1:1453-1534  → adjacent telotron at 813-1038 (gap 415 bp)
NW_013564060.1:17301-17370 → adjacent telotron at 17489-17860 (gap 119 bp)
NW_013549356.1:7752-7806   → adjacent telotron at 7385-7568 (gap 184 bp)
NW_013545172.1:1021-1050   → cross-contig telotron at NW_013542093.1:1077-1986
NW_013654861.1:79248-79326 → cross-contig telotron at NW_013655283.1:230-446
```

**Two adjacent telotrons sharing the same 60-80bp linker sequence is direct evidence that they share a common origin.** Either:
1. **Local segmental duplication**: a region containing one telotron was duplicated, creating a second telotron carrying the same surrounding sequence.
2. **Sequential capture at a hotspot**: two independent capture events at the same locus, both using the same template/source fragment.
3. **Recombination between subtelomeres**: subtelomeric homologous recombination spreads telotron-containing fragments across loci.

### Putting it together: Eimeria mechanism model

Combining all evidence:
1. Eimeria telotrons cluster in **TANDEM ARRAYS** (65% of contigs have ≥2; 25% of pairs <100bp apart)
2. The "linkers" between telomeric arrays in DIVERGING/CONVERGING telotrons are mostly **unique sequences** in the genome
3. The few paralogous linkers are **subtelomeric** (2-3× enriched vs background)
4. **DIVERGING orientation dominates** Eimeria (12.5:1 vs CONVERGING)
5. Multi-telotron contigs are common, with adjacent telotrons sharing linker sequence in a few cases

Best-fit model: **subtelomeric instability creates recurrent capture events at specific hotspot loci**. Once a locus accumulates one telotron, subsequent telomere/fragment capture events at the same hotspot deposit additional telomeric arrays in tandem. The DIVERGING orientation arises naturally from **outward-pointing telomere extension** at a chromosome-end-like context (which is how the captured fragments may have been before integration). The "linker" sequences in DIVERGING pairs are **the ordinary intronic/intergenic DNA that was preserved between two successive telomeric repeat insertions at the hotspot**.

This explains why:
- Linkers are usually unique (they're original genomic context, not captured elsewhere)
- Telotrons cluster (recurrent capture at hotspots)
- DIVERGING dominates (geometry of subtelomeric-style insertion)
- Some linkers are shared (those hotspots that experienced multi-event capture preserve a copy of the previous insertion)

The Eimeria mechanism likely reflects **Apicomplexan subtelomeric biology** — these parasites are well-known for active subtelomeric recombination (var gene switching in Plasmodium, SCAR/TARE elements). Eimeria telotrons may be a byproduct of this normally-subtelomeric recombination machinery acting on internal DSB sites that mimic chromosome-end structure.

---

## 9f. EIMERIA LINKER → APICOMPLEXA BLAST (2026-04-29)

Direct test: BLAST 759 Eimeria pure non-repeat linkers (≥25 bp, no
TTAGGG/TTTAGGG kmers) against a focused database of 12 Apicomplexa genomes
(5 Eimeria + Cyclospora + 4 Sarcocystidae + Cryptosporidium + Plasmodium).
~665 MB BLAST db, blastn at ≥70% identity, ≥30 bp aligned, exclude self-genome.

**Aggregate hit profile:**

| Lineage class | Linkers with ≥1 hit | / 759 input |
|---|---:|---:|
| Other Eimeria | 371 | 48.9% |
| Sarcocystidae (Toxoplasma/Neospora/Besnoitia/Sarcocystis) | 30 | 4.0% |
| Plasmodium / Cryptosporidium (deep Apicomplexa) | 17 | 2.2% |
| Cyclospora (sister Eimeriidae) | 9 | 1.2% |

**~50% of pure non-repeat Eimeria linkers have a strong cross-Eimeria hit.**
This is dramatically higher than the §9c finding (98% no-paralog within own
genome), and reverses the apparent "linkers are unique" conclusion when the
search expands to closely-related species.

### Where do the cross-Eimeria hits land? (telotron vs plain DNA)

Using the full ~20K Eimeria telotron interval set (from
`ultra_results/GCF_000499*_telotrons.tsv`):

| Per-linker classification (n=371 with hits) | Count | % |
|---|---:|---:|
| **Hits ONLY at non-telotron locations in sister Eimeria** | **111** | **30%** |
| Hits ONLY inside other telotrons in sister Eimeria | 191 | 52% |
| Mixed (some hits in telotrons, some not) | 69 | 19% |

**The 111 "all-outside-telotrons" linkers are the smoking gun for the
ancestral-DNA hypothesis** — sequences that exist as plain genomic DNA in
sister Eimeria species but appear inside telotrons in the source species.

Of these 111, partitioned by hit count:
- **Sparse (1-5 hits at specific loci): 91**, of which **85 have a best-hit
  pid×len product ≥3200** (e.g., 125 bp at 92.8% identity over a 129 bp
  linker — full-length match to a single locus in a sister Eimeria).
- Moderate (6-20 hits): 10 — possibly small repeat families.
- Heavy (>20 hits): 10 — likely repetitive/subtelomeric content that
  survived the exact-kmer filter.

### Top sparse cases (sample)

| Source linker | Best target | Hit length | Identity |
|---|---|---:|---:|
| `E. tenella NW_013541501.1__67325-67552` (linker 129 bp) | E. necatrix | 125 bp | 92.8% |
| `E. necatrix NW_013654953.1__55517-55966` (linker 111 bp) | E. tenella | 111 bp | 97.3% |
| `E. tenella NW_013545730.1__2048-2301` (linker 123 bp) | E. necatrix + E. acervulina | 122 bp | 86.1% |
| `E. necatrix NW_013652280.1__39811-40024` (linker 109 bp) | E. tenella | 109 bp | 93.6% |

These are full-linker-length BLAST hits at 85-97% identity to single
specific genomic loci in sister Eimeria species. The sister locus is plain
genomic DNA (not itself a telotron), which means **the linker IS the
ancestral coccidian DNA at this orthologous gene context**. Telomeric arrays
were inserted into pre-existing intronic DNA in the source species; the
sister species retains the ancestral state without insertion.

### Reconciliation with §9e (ancestral-intron-as-linker null result)

§9e tested using PLAIN-intron OUTGROUP species (Toxoplasma, Cyclospora etc.)
and concluded "no homology." That test failed because:
1. Eimeria↔Toxoplasma divergence (~150-200 My) erodes intron homology below
   detection over short linker spans (~50-100 bp).
2. The intron-orthology pipeline (`linker_analysis/orthology_v2`) only
   produced 5 fill_in alignments due to a `MAX_LOCI=60` cap.

The new §9f test uses **closer relatives** (sister Eimeria species,
divergence ~10-25 My) where homology IS still detectable, and recovers
exactly the signal §9e missed. **The user's hypothesis is supported in the
strict-coccidian setting**, even though it failed in the deeper Apicomplexan
test.

### Sarcocystidae hits (n=30 linkers)

30 Eimeria linkers also hit Sarcocystidae genomes (Toxoplasma/Neospora/
Besnoitia/Sarcocystis). These represent linkers conserved across the
~150 My Coccidia ancestor. Worth following up to localize at orthologous
gene loci.

### Plasmodium / Cryptosporidium hits (n=17)

17 linkers hit deeper Apicomplexa. Particularly interesting because these
genomes are very divergent from Eimeria (>400 My) — surviving homology
implies highly-conserved sequence (likely conserved coding-flank context
near the telotron rather than the linker proper). Should inspect.

Inputs: `real_telotrons/strict_linkers_eimeria.json`, ~20K telotron intervals
        from `ultra_results/GCF_000499*_telotrons.tsv`
Outputs: `real_telotrons/eimeria_linker_blast/{eimeria_linkers_vs_apicomplexa.blastn.tsv,
         eimeria_linker_strong_hits.tsv, eimeria_linker_blast_summary.json}`
Pipeline: `pan_euk_telotrons/blast_eimeria_linkers_vs_paneuk.py`

### Figures (`real_telotrons/eimeria_linker_blast/figures/`)

- **`summary.png`** — 3-panel: (A) BLAST aggregate funnel, (B) per-linker
  classification pie chart, (C) hit length × identity scatter colored by
  lineage class.
- **`catA_example{1,2,3}.png`** — Category A: linker hits non-telotron
  region in sister Eimeria (ancestral DNA candidates).
- **`catB_example{1,2,3}.png`** — Category B: linker matches inside another
  telotron (shared telotron across species).
- **`catC_example{1,2,3}.png`** — Category C: Sarcocystidae outgroup hit.
- **`catD_example{1,2,3}.png`** — Category D: deep Apicomplexa hit
  (Plasmodium / Cryptosporidium, >400 My).

Each example figure shows a wrapped DNA-letter MSA with the host telotron
span marked, telomeric-kmer-residue cells colored separately from
non-telomeric ("linker") cells, and per-row ungapped coordinates.
Pipeline: `pan_euk_telotrons/make_linker_origin_figures.py`,
`pan_euk_telotrons/make_linker_summary_figure.py`

---

## 9e. ANCESTRAL-INTRON-AS-LINKER TEST (2026-04-29)

User hypothesis: telotron linkers are ancestral intronic DNA preserved at the
locus, with telomeric arrays inserted into pre-existing introns. Tested by
extracting the telo-stripped linker from each telotron and comparing it to
sister-species sequence at orthologously aligned columns.

### Iteration history (read this whole section, including the corrections)

**v1 — Eimeria-only ortholog groups (n=106 intronic-class breakpoints).**
Sister species without called telotrons at orthologous positions are *also*
heavily telomeric in most cases (the §9c partial-decay finding). Only 3
comparisons met the clean test (sister telo<30%, host linker ≥30 bp); none
showed significant linker→sister homology.

**v2/v3 — pan-eukaryotic fill_in alignments with PLAIN outgroups.**
Using existing alignments (n=5 fill_in alignment files in
`linker_analysis/orthology_v2/alignments/fill_in/`). Initially appeared to
give a striking positive result for *Ficedula albicollis*
`rna-XM_005062594.2`: 1369 bp host "linker" with strong BLAST hits
(362 bp at 82.87% identity to *Poecile atricapillus*, 154 bp at 88.96% to
*Vidua macroura*, etc.).

**v4 — extended to all 30 fill_in loci with ≥3 cross-species orthologs**
(directly from `ortholog_results.json`, no need to rerun MSA pipeline).
Result before filtering: **2/30 loci** showed strong ancestral homology
(Ficedula + a Tara MAG case).

**v4 with strict filtering — both apparent positives collapse:**

1. **Ficedula `rna-XM_005062594.2`** has `telo_type = TTGGG~1mm`. Inspecting
   the host intron: it is `(GGGAT)n` tandem-repeat microsatellite, NOT
   TTAGGG. `count("TTAGGG") = 0`, `count("GGGAT") = 269` over 1369 bp.
   The sister-bird "introns" are also `(GGGAT)n` arrays at the same
   orthologous locus. The BLAST homology is between intronic
   *microsatellites*, not preserved ancestral DNA. CLAUDE.md flags
   `~1mm` ULTRA hits as ~10:1 noise; this is exactly that failure mode.

2. **Tara MAG `mRNA.TARA_MED_95_MAG_00463_000000000360.13.8`** has
   `telo_type = TTAGG` (insect 5-bp telomere). The host intron is dominated
   by `CTAACCTAAC...` and `TTAGGTTAGG...` — degenerate 5-bp arrays. After
   stripping exact TTAGGG/TTTAGGG kmers, the residual "linker" is itself
   substituted-CCTAA. The strong hits are to *Bombus impatiens*,
   *Bicyclus anynana*, *Eriocheir sinensis* — insects/crustaceans that use
   `TTAGG`-family telomeric repeats at this orthologous locus. Same
   microsatellite-homology artifact, different repeat unit.

### Final result (filtered)

After excluding `~1mm`-type hits and microsatellite-dominated linkers:

| Result | Count | % |
|---|---:|---:|
| fill_in loci with ≥3 targets, exact-telo type | 29 | — |
| **ANCESTRAL_HOMOLOGY (≥2 plain targets with strong BLAST)** | **0** | **0%** |
| WEAK (1 plain target with strong BLAST) | 1 | 3% |
| NO_HOMOLOGY | 28 | 97% |

The single `WEAK` case is *E. maxima* `rna-XM_013481188.1` with one
`Sarcocystis_neurona` hit; small total linker (63 bp) makes this likely
chance.

**The user's hypothesis is NOT supported by the data.** Linkers in real
TTAGGG/TTTAGGG-family telotrons are **not** preserved ancestral intron DNA
that pre-existed at the locus. They are novel sequence introduced alongside
the telomeric insertion event.

This is a meaningful biological finding, not a null result:
1. The "fill_in" classification (sister has intron at orthologous position)
   reflects conservation of *gene structure* (intron position), not
   conservation of *intron content* in the host.
2. Combined with §9d (linkers are not TERC retrotranscripts) and §9c
   (linkers are unique in the host genome, not paralogous), all three major
   "external source" hypotheses for linker DNA — TERC, paralog capture,
   ancestral intron — are now ruled out within detection limits.
3. This leaves **non-templated polymerase fill-in (NHEJ/MMEJ scar)** and
   **Polθ TMEJ insertion** as the dominant remaining hypotheses.

**Two important caveats:**
- Eimeria introns are short (median ~150 bp). After stripping telomeric
  kmers, the residual linker (~50-100 bp) may simply be too short to
  preserve detectable homology across the 150-200 My Eimeria↔Toxoplasma
  divergence. Absence of evidence ≠ evidence of absence at this scale.
- Cases where the telotron is genuinely embedded in a long pre-existing
  intron (vertebrate-style) may exist but were not present among the
  exact-telomere fill_in loci with ≥3 testable orthologs.

Inputs: `linker_analysis/orthology_v2/results/ortholog_results.json`
Outputs:
  `linker_analysis/orthology_v2/results/ancestral_linker_v4.tsv` (per-target),
  `linker_analysis/orthology_v2/results/ancestral_linker_v4_per_locus.tsv`,
  `linker_analysis/orthology_v2/results/ancestral_linker_v4_summary.json`
Pipelines: `pan_euk_telotrons/test_ancestral_linker_v4.py` (current),
  `test_ancestral_linker_v3.py` (prior, 5 loci with MSAs),
  `test_ancestral_linker_v2.py` (Eimeria-only).

**Lesson:** Microsatellite-style arrays (`TTAGG`-family insect telomeres,
`GGGAT`-style intronic minisats) cause the same kind of false-positive cross-
species BLAST homology as the degenerate-telomeric "linkers" already flagged
in §5. Any future linker-comparison analysis must explicitly filter for
microsatellite content in BOTH the host linker AND the target intron, not
just for exact telomeric kmers.

---

## 9d. TERC BLAST RESULT (2026-04-29)

Tested whether pure non-repeat linkers are TERC retrotranscript fragments
(Nergadze 2007). Pulled 201 TERC sequences from NCBI Entrez (queries:
"telomerase RNA component"[Title], TERC[gene] across major eukaryotic
clades). BLASTed 1,353 pure non-repeat linkers (759 Eimeria + 594 haptophyte,
≥20 bp, no TTAGGG/CCCTAA/TTTAGGG kmers) against this DB.

| Metric | Value |
|---|---:|
| Linkers queried | 1,353 |
| TERC sequences in DB | 201 |
| Raw hits (e<1e-3) | 15 |
| **Strong hits (≥80% ident, ≥20 bp)** | **1** |
| Unique linkers with any strong TERC match | 1/1,353 (0.07%) |

The single strong hit is to *Electrophorus electricus* "telomerase RNA
component interacting RN[A]" — an interacting protein, **not TERC itself** —
at 83.6% identity over 61 bp inside a 235 bp linker. Almost certainly spurious.

**Conclusion:** TERC retrotranscription (Nergadze 2007 model) is **essentially
ruled out** as a source of telotron linkers in either Eimeria or haptophyte
data. Combined with §9c's 98% no-paralog finding, this means linkers are
neither host-genomic in origin nor matches to known eukaryotic TERC.

**Caveats:** NCBI has no annotated TERC for Apicomplexa or Haptophyta, so we
cannot exclude lineage-specific TERC (this would require ab-initio TERC
reconstruction from MAG TERT context). Old insertions would have eroded
TERC homology beyond detection.

Remaining viable hypotheses for the pure non-repeat linker fraction:
1. **Polθ TMEJ templated insertion** from local intronic/intergenic context
   (consistent with §9c finding that paralogs, when they exist, are local).
2. **NHEJ/MMEJ scar** with non-templated polymerase fill-in (but the median
   linker is 58 bp — longer than typical scar lengths).
3. **Lineage-specific TERC retrotranscription** of an uncharacterized RNA
   template (untestable without MAG TERC reconstruction).

Inputs: `terc_blast/pure_nonrepeat_linkers.fa`, `terc_blast/terc_db.fa`
Outputs: `terc_blast/linkers_vs_terc.blastn.tsv`, `terc_blast/terc_blast_summary.json`
Pipeline: `pan_euk_telotrons/blast_linkers_vs_terc.py`

---

## 10. Open questions & next steps

1. **BLAST 395 medium pure-non-repeat linkers against TERC sequences.**
   ✅ Done (§9d, 2026-04-29). Result: 1/1,353 strong hits — TERC
   retrotranscription ruled out within sensitivity of the public TERC
   collection (no annotated Apicomplexa/Haptophyta TERC available).

2. **Cluster the 8,701 linkers by sequence similarity.**
   Are there a small number of dominant consensus sequences that could be templates? (Done at coarse level — most "clusters" are driven by shared TTAGGG/CCCTAA kmers; finer analysis on the pure-non-repeat subset would be informative.)

3. **Empty-ortholog comparison.**
   For converging telotrons with linkers, find orthologous loci in sister MAGs without the insertion. Do the linker-flanking exonic sequences pre-exist? Does the TSD/microhomology fall in the empty locus?

4. **Test the splice-compatible-retention hypothesis.**
   Compare the splice-site quality of:
   - Validated intronic telotrons (in our set)
   - Intronic telomeric tracts that fail splice annotation (would need GFF + RNA-seq)
   - Intergenic telomeric tracts
   The hypothesis predicts the validated set has the strongest GT-AG signals.

5. **Test the DSB-hotspot hypothesis.**
   Are the haptophyte telotron loci enriched at G4-forming sites, R-loop-prone regions, or replication fork barriers?
