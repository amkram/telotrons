# Manuscript Validation Report: Telotron Claims vs Data

## Summary

Systematically validated all quantitative claims in the manuscript against the underlying data files. Of ~55 distinct claims checked, **most are verified**, but **several critical discrepancies** require attention.

---

## VERIFIED CLAIMS (No Issues)

| Claim | Manuscript | Data | Status |
|-------|-----------|------|--------|
| Total MAGs | 257 | 257 | ✓ |
| Supergroups | 7 | 7 | ✓ |
| MAGs with ≥50 telo introns | 30 | 30 | ✓ |
| 4 supergroups in top 30 | Cryptista, Haptista, Opisthokonta, Stramenopiles | Matches | ✓ |
| Primary MAG total introns | 58,890 | 58,890 | ✓ |
| High-purity in primary MAG | 844 | 844 (146 very_strict + 698 strict) | ✓ |
| Very high-purity (>95%) in primary | 146 | 146 | ✓ |
| At 85–95% in primary | 698 | 698 | ✓ |
| Telo-containing (10–85%) in primary | 32,974 | 32,974 | ✓ |
| Haptophyta burden | 89.5% | 89.5% (46,121/51,532) | ✓ |
| Converging telotrons in primary | 462 (54.7%) | 462 (54.7%) | ✓ |
| Splicing: 2,542 genes both types | 2,542 | Verified | ✓ |
| Splicing: 46.7% vs 41.8% | 46.7 vs 41.8 | Verified | ✓ |
| Fisher's P | 3.0 × 10⁻¹³ | Verified | ✓ |
| Orthologs: 52 GT-AG | 52 | Verified | ✓ |
| Orthologs: 116 pairwise | 116 | Verified | ✓ |
| Confirmation rate | 6.2% (52/844) | 6.16% | ✓ |
| TSDs at 71% | 71% | Verified | ✓ |
| CTAA boundary enrichment | 87–325× | 4.4–12.2% vs ~0.04% background | ✓ |
| TERT 10 of 20 MAGs | 5 Tier A + 5 Tier B | Verified | ✓ |
| Burden-youth Spearman | r = 0.73, P = 3.8 × 10⁻¹⁰ | Verified | ✓ |
| Within-gene G4 Wilcoxon | P = 5.0 × 10⁻¹² | Verified | ✓ |
| G4 enrichment | 13% higher | Verified | ✓ |
| Within-Haptista Pif1 | ρ = −0.82, P = 0.023, n = 7 | Exact match | ✓ |
| Purity medians | 96.5% very high, 89.7% high | Consistent with data | ✓ |
| ~50,000 degraded (abstract) | ~50,000 | ~49,000 (51,532 − 2,371) | ✓ |

---

## CRITICAL DISCREPANCIES

### 1. Strand Configuration Breakdown (PRIMARY MAG)

**Severity: HIGH** — Two of three per-category numbers are wrong.

| Configuration | Manuscript | Data | Difference |
|--------------|-----------|------|------------|
| Converging | 462 (54.7%) | 462 (54.7%) | ✓ Match |
| Template-only | **372 (44.1%)** | **344 (40.8%)** | **+28 overcounted** |
| Coding-strand-only | **10 (1.2%)** | **30 (3.6%)** | **−20 undercounted** |
| Chimeric/other | (not mentioned) | 8 (0.9%) | **Missing category** |
| **Total** | 844 | 844 | ✓ Match |

The non-converging total is correct (382), but the split between template-only and coding-strand-only is wrong. The manuscript overcounts template-only by 28 and undercounts coding-strand-only by 20, with 8 chimeric entries unaccounted for.

**Impact:** The "1.2%" coding-strand-only claim is central to the splice-signal-noise argument. The actual rate is 3.6% — still rare, but 3× higher than stated.

### 2. Template-Strand Preference Ratio

**Severity: MEDIUM**

| Metric | Manuscript | Data |
|--------|-----------|------|
| Template-strand preference | "8.8:1" | 11.5:1 (high-purity) or 7.5:1 (all telo) |
| Template-strand bias | "81%" in haptophytes | 92.0% (high-purity) or 88.3% (all) |

The "8.8:1" ratio doesn't match any available computation. The high-purity ratio is 11.5:1 (344:30); the all-intron ratio is 7.5:1 (26,218:3,484). Neither is 8.8:1.

The "81% template-strand bias" claim also doesn't match: high-purity single-strand class is 92% template, all-intron single-strand is 88%.

### 3. GT-AG Splice Site Rate

**Severity: MEDIUM** — Numerator and denominator differ from data.

| Metric | Manuscript | Data |
|--------|-----------|------|
| GT-AG rate | 674/721 = 93.5% | 800/844 = 94.8% |
| GT donors | "89% across 8 MAGs" | 94.8% in primary |
| AG acceptors | "95% across 8 MAGs" | 95.7% in primary |

The "721 with annotated splice boundaries" denominator doesn't match the data — all 844 high-purity telotrons have splice site annotations. The actual GT-AG rate is 94.8%, not 93.5%.

### 4. 2,371 High-Purity Across 8 MAGs

**Severity: MEDIUM** — Cannot fully verify.

The 844 in the primary MAG is confirmed. The 2,371 total across 8 MAGs and "248 very high-purity" cannot be verified from available data because the strand_deep_analysis_v2.pkl contains only primary MAG data (40,164 entries). The telo_purity.npy shows 2,239 with purity > 0.85 across all 257 MAGs (not just 8), but structural filtering reduces this differently.

### 5. "248 very high-purity (>95%): 146 primary + 102 remaining"

**Severity: LOW** — 146 in primary is confirmed. 102 in other 7 MAGs is unverifiable from current data.

---

## CLAIMS WITH DATA GAPS (Cannot Verify)

### 6. Repair Gene Correlations Across 48 MAGs

**Severity: HIGH** — The data to verify these doesn't exist in the available files.

| Claim | Manuscript | Available Data |
|-------|-----------|----------------|
| Ku70 ρ = −0.67, P = 2.1 × 10⁻⁷ | "across 48 MAGs" | Only 8 MAGs in repair_gene_hmm.tsv |
| Pif1 ρ = −0.41, P = 4.3 × 10⁻³ | "across 48 MAGs" | Only 8 MAGs |
| Within-Haptista Ku70 | ρ = −0.68, P = 0.002, n = 18 | No Ku70 data for 18 Haptista |
| 6.3% vs 5.8% repair disruption | — | No data found |

The repair gene HMM analysis was only run on the top 8 MAGs. Correlations across 48 MAGs and within 18 Haptista MAGs require repair gene data for many more MAGs than currently exist. Either this analysis was done externally or the numbers are from a different computation.

### 7. "53 Contain Telomeric-Repeat Introns"

**Severity: LOW-MEDIUM** — Technically incorrect as stated.

All 257 MAGs have total_telos > 0 (range: 1–19,558). The "53" corresponds to total_telos ≥ 20, but the manuscript says "53 contain telomeric-repeat introns" without specifying a threshold. Should say "53 carry ≥20 telomeric-repeat introns" or similar.

### 8. "48 Have Sufficient Gene Model Counts"

**Severity: LOW** — All 257 MAGs have gene counts. The criterion for "sufficient" is undefined.

---

## TAXONOMY DISCREPANCIES

| Metric | Manuscript | per_mag_telotron_statistics (257 MAGs) | telotron_metrics_detailed (605 MAGs) |
|--------|-----------|--------------------------------------|--------------------------------------|
| Phyla | 18 | 17 | — |
| Classes | 32 | 25 | 32 |
| Genera | "up to 49" | 25 | 49 |

The "32 classes" and "49 genera" match the 605-MAG dataset, not the 257 MAGs described in the sentence. The 257 MAGs have only 25 classes and 25 genera. Either the sentence should reference the larger dataset or the counts should be corrected to 25/25.

The 18 vs 17 phyla discrepancy is unexplained — no data source gives exactly 18.

---

## ANALYSIS CROSS-CHECK ISSUES

The ANALYSIS_RESULTS_EXACT_NUMBERS.txt (generated during this session) reports **837** high-purity telotrons with a completely different strand breakdown (549 single_rev, 58 single_fwd, 101 converging, 129 chimeric = 837). This is because that analysis used a simple purity > 0.85 filter on a different subset (the strand_deep_analysis had the telotrons matched to a different classification), NOT the tier-based structural criteria (very_strict + strict = 844). The per-configuration GT-AG rates reported in that analysis should not be used.

---

## RECOMMENDED FIXES (Priority Order)

1. **Strand counts**: Correct template-only to 344, coding-strand-only to 30, add note about 8 chimeric/ambiguous
2. **Template-strand ratio**: Correct "8.8:1" to the actual ratio (11.5:1 in high-purity or 7.5:1 in all)
3. **Template-strand bias**: Correct "81%" to actual value from the data
4. **GT-AG rate**: Verify the 721 denominator — if all 844 have annotations, report 800/844 = 94.8%
5. **Taxonomy**: Correct to 17 phyla, 25 classes, 25 genera for the 257 MAGs (or clarify the 605-MAG source)
6. **"53 contain"**: Add threshold qualifier (≥20 telo introns)
7. **Repair gene correlations**: Run HMM analysis on additional MAGs to verify n=48 and n=18 claims, or note that these are from a separate analysis
