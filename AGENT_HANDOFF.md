# Agent Handoff: Telotron Project

**Date:** 2026-04-06
**Repo:** https://github.com/amkram/telotrons
**Working directory:** `/Users/alex/telotrons/`

---

## What This Project Is

We discovered **telotrons** — introns composed of tandem telomeric repeats (TTAGGG, TTTAGGG, etc.) that are created when telomerase repairs double-strand breaks in coding sequence. The inserted repeat array gets spliced out as an intron. First found in marine eukaryotic MAGs (Tara Oceans), now confirmed across multiple eukaryotic lineages via a pan-eukaryotic survey of NCBI reference genomes.

A Nature Genetics manuscript is in preparation.

---

## Current State (what's running / what just finished)

### ULTRA Survey v2 (IN PROGRESS)
- **Script:** `pan_euk_telotrons/run_ultra_survey.py`
- **Status:** ~842/1697 genomes complete, ~855 remaining (mostly large vertebrate/plant genomes that fail to download)
- **Process:** `nohup python3 run_ultra_survey.py >> /tmp/ultra_survey.log 2>&1 &`
- **Config:** 3 workers, 16 ULTRA threads each
- **Checkpoint:** `pan_euk_telotrons/ultra_checkpoint.json`
- **Output:** `pan_euk_telotrons/ultra_results/` — per-species TSV of telomeric-repeat introns

The v2 survey uses ULTRA (a fast tandem repeat finder) instead of custom hexamer scanning. It searches **13 known telomeric repeat types** (TTAGGG, TTTAGGG, TTAGG, TTTGGG, TTGGGG, etc.) and their 1-mismatch variants. The old v1 survey only searched TTAGGG.

**To check progress:**
```bash
python3 -c "import json; ckpt=json.load(open('pan_euk_telotrons/ultra_checkpoint.json')); print(f'{len(ckpt)}/1697')"
tail -5 /tmp/ultra_survey.log
```

**To resume after interruption:**
```bash
cd pan_euk_telotrons
# Clear download failures (they retry):
python3 -c "
import json
ckpt = json.load(open('ultra_checkpoint.json'))
ckpt = {k:v for k,v in ckpt.items() if v.get('status') == 'ok'}
json.dump(ckpt, open('ultra_checkpoint.json','w'), indent=1)
"
nohup python3 run_ultra_survey.py >> /tmp/ultra_survey.log 2>&1 &
```

### What HASN'T been done yet on the ULTRA results
1. **Post-filtering:** The raw ULTRA hits include massive noise (~1mm false positives). Need to filter by:
   - TR coverage ≥ 50% or ≥85% of intron length
   - Exact consensus match to a known telo repeat (not just 1-mismatch)
   - Canonical splice sites (GT-AG or CT-AC)
2. **Misannotation filtering** (flanking exon checks) — script exists at `pan_euk_telotrons/filter_misannotation.py`
3. **Rebuild `real_telotrons/` folder** with the ULTRA-based results
4. **Regenerate the figure** (`telotron_map_all_species.png`)

---

## Key Results So Far

### Validated telotrons (from v1 survey + misannotation filter)
- **186 validated telotrons across 59 species** in `pan_euk_telotrons/real_telotrons/`
- **Eimeria (Apicomplexa):** 103 telotrons across 5 species — the strongest signal outside Haptista
  - Uses TTTAGGG (7-mer, plant/apicomplexan type), not TTAGGG
  - 9/9 E. tenella splice junctions validated by mRNA in NCBI (100% match to XM_ accessions)
  - 41/92 genes have functional protein annotations (RNA polymerase, importin, kinases, etc.)
  - Flanking exons are non-telomeric, genes are multi-exon, all protein-coding
- **Monocercomonoides:** 289 raw → 1 after filter (288 were misannotated split telo arrays)
- **Most insect/vertebrate hits:** misannotations (telomeric exons) or 1-mismatch noise

### Human genome search
- **T2T-CHM13:** 0 telotrons at ≥85%. Best candidate 57.6% (degraded ITS in NR_137167.1 on chr5)
- **1KG ONT Vienna SVs:** 14 high-purity telomeric insertions across 967 humans
  - 4 fall within gene introns (ITPR1, CTNND2, IRAG2, SCGB1C2)
  - All are **insertions INTO existing introns** — NOT telotrons (no splice signals from repeat)
  - Template-strand bias (86%), no TSDs, consistent with telomerase DSB repair
- Results in `human_telotron_survey/`

### Tara Oceans MAGs (the original discovery)
- 844 high-purity telotrons in primary MAG (TARA_PSW_86_MAG_00284)
- 32,974 telo-containing introns across 257 MAGs
- Data in `tara_oceans_euk_mags/` (92 GB — gitignored)
- Pipeline docs: `tara_oceans_euk_mags/LLM_ONBOARDING.md`

---

## Directory Structure

```
telotrons/
├── .gitignore
├── AGENT_HANDOFF.md              ← you are here
├── plan.md                       # Pan-eukaryotic survey plan
│
├── # Manuscript & reviews
├── TRI_manuscript_draft*.docx    # Nature Genetics manuscript drafts
├── TRI_Extended_Data.*           # Extended data figures
├── REVIEW_NatGenetics_v3.md      # Latest internal review (minor revisions)
├── MANUSCRIPT_AUDIT_REPORT.md    # Claim-by-claim audit with verified numbers
├── LITERATURE_VALIDATION_REPORT.md
├── STRUCTURAL_CRITIQUE.md
├── critical_review_nat_genet.md  # Pre-submission review
├── VALIDATION_REPORT.md
│
├── # Figures (ED = Extended Data)
├── ED_Fig*.png                   # Extended Data figures
├── Figure1_composite.*           # Main figures
├── Figure2_composite.*
│
├── # Pan-eukaryotic survey
├── pan_euk_telotrons/
│   ├── run_survey.py             # v1 survey (TTAGGG-only hexamer scan)
│   ├── run_ultra_survey.py       # v2 survey (ULTRA, 13 repeat types) ← CURRENT
│   ├── filter_misannotation.py   # Validates against GFF gene structure
│   ├── validate_telotrons.py     # Basic validation (splice, purity, tandem)
│   ├── rescreen_multirepeat.py   # Re-classify by correct repeat type
│   ├── checkpoint.json           # v1 survey results (2312 genomes, gitignored)
│   ├── ultra_checkpoint.json     # v2 survey results (in progress, gitignored)
│   ├── telotrons/                # v1 raw candidate TSVs (144 species)
│   ├── ultra_results/            # v2 raw candidate TSVs (gitignored, large)
│   ├── validated_results/        # Misannotation filter output
│   ├── real_telotrons/           # Final validated telotrons (186 across 59 species)
│   │   ├── SUMMARY.tsv
│   │   ├── all_real_telotrons.json
│   │   ├── telotron_map_all_species.png  # Main figure
│   │   └── *_real_telotrons.tsv  # Per-species TSVs with flanks + sequences
│   └── bin/
│       └── ultra                 # ULTRA binary (compiled from source)
│
├── # Human genome search
├── human_telotron_survey/
│   ├── run_human_survey.py       # Search annotated human assemblies
│   ├── 1kg_ont_vienna/
│   │   ├── scan_sv_telotrons.py  # Search 1KG ONT SVs for telo insertions
│   │   └── results/              # SV scan results
│   └── results/                  # Per-assembly results
│
├── # Tara Oceans data (92 GB, gitignored)
├── tara_oceans_euk_mags/         # MAG genomes, intron catalogs, validation
│
├── # Analysis scripts
├── analyze_terminals_minimal.py  # Terminal sequence analysis
├── analyze_telotron_terminals*.py
└── *.py, *.txt, *.md             # Various analysis outputs
```

---

## Known Issues in the Manuscript

See `MANUSCRIPT_AUDIT_REPORT.md` for full details. Critical fixes needed:

| Issue | Current | Correct |
|---|---|---|
| TSD percentage | 71% | ~17% (or remove) |
| Phase chi-square | 324.4, n=218 | 626.7, n=844 |
| Phase 5 enrichment | 58%, 3.5-fold | 47.6%, 2.9-fold |
| Template-only count | 372 (44.1%) | 344 (40.8%) |
| Coding-only count | 10 (1.2%) | 30 (3.6%) |
| GT-AG rate | 674/721 = 93.5% | 800/844 = 94.8% |
| "TTAA in opisthokonts" | stated | unverifiable (no data) |
| Format | Brief Communication | Should be Article (exceeds word limits) |

---

## Key Tools & Dependencies

| Tool | Location | Purpose |
|---|---|---|
| ULTRA | `pan_euk_telotrons/bin/ultra` | Fast tandem repeat detection |
| NCBI datasets | `../pan_euk_telotrons_datasets` (binary) | Download genomes + GFF3 |
| Python 3 | system | All scripts |
| matplotlib | pip | Figure generation |
| numpy | pip | Data processing |

ULTRA was compiled from https://github.com/TravisWheelerLab/ULTRA (v1.2.1).

---

## Telomeric Repeat Types Searched

| Repeat | Length | Organisms |
|---|---:|---|
| TTAGGG | 6 | Vertebrates, fungi, most eukaryotes |
| TTTAGGG | 7 | Plants, Apicomplexa (Eimeria, Toxoplasma) |
| TTAGG | 5 | Insects |
| TTTGGG | 6 | Paramecium, some protists |
| TTGGGG | 6 | Tetrahymena |
| TTTGGGG | 7 | Vorticella |
| TTTTGGGG | 8 | Oxytricha, Euplotes |
| TTTTAGGG | 8 | Genlisea, some Charophytes |
| TTTTTTAGGG | 10 | Chara |
| TTAGGC | 6 | Nematodes (C. elegans) |
| TAGGG | 5 | Giardia |
| TTCAGG | 6 | Chlamydomonas |
| TTGGG | 5 | Colpidium |

For each, all rotations + reverse complement are searched (e.g., TTAGGG → 12 variants).

---

## What To Do Next

### Immediate (when ULTRA survey finishes)
1. Post-filter ULTRA results: keep only introns with ≥50% TR coverage from an exact telo consensus (drop all `~1mm` noise)
2. Run `filter_misannotation.py` on the filtered set
3. Rebuild `real_telotrons/` with updated results
4. Regenerate figure
5. Push to GitHub

### For the manuscript
1. Add pan-eukaryotic results (Eimeria as second lineage with telotrons)
2. The Eimeria finding changes "telotrons are haptophyte-specific" to "telotrons occur across eukaryotic supergroups"
3. TTTAGGG repeat in Eimeria → different TERC template boundary signature (GTTT|AG vs GGTT|AG)
4. Consider switching from Brief Communication to Article format
5. Apply all corrections from MANUSCRIPT_AUDIT_REPORT.md

### Stretch goals
- Retry failed genome downloads (vertebrates, large plants) with longer timeouts
- Run ULTRA survey on remaining ~276 genomes that the v1 survey also failed on
- Build phylogenetic figure showing telotron distribution across tree of life
- Search for telotrons in pangenome/T2T assemblies of model organisms
- Investigate whether Eimeria telotrons are ortholog-confirmed de novo insertions

---

## Gotchas & Lessons Learned

1. **Misannotation is the #1 false positive source.** Gene predictors split long telomeric arrays into fake exons + introns. Always check flanking exon telomeric content.

2. **Species-specific repeats matter.** Eimeria uses TTTAGGG not TTAGGG. The v1 survey found them by accident (TTAGGG is a substring). Always search with the correct repeat.

3. **1-mismatch matching produces massive noise.** In the ULTRA survey, `~1mm` hits outnumber exact hits ~10:1. Post-filter aggressively.

4. **Human "telotrons" are passengers.** The 4 intronic telomeric insertions found in 1KG ONT data are insertions INTO existing introns, not introns CREATED BY the insertion. CCCTAA repeats lack GT and AG, so they can't create splice signals.

5. **Parallel NCBI downloads fail.** >10 concurrent `datasets download` calls overwhelm NCBI. Use ≤5 workers for download-heavy runs.

6. **Large genome downloads time out.** Vertebrate genomes (>1 GB compressed) need 3600s+ timeout. The `download_genome()` function in all scripts has a `timeout` parameter.

7. **ULTRA is fast but memory-hungry with many threads.** Memory scales linearly with thread count. Use `-t 16` max, not `-t 64`.
