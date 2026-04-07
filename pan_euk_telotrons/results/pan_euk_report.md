# Pan-Eukaryotic Telotron Survey — All NCBI RefSeq Eukaryotes

**Date**: 2026-04-02
**Genomes in catalog**: 2,312 annotated eukaryotic RefSeq assemblies
**Successfully processed**: 1,697 (73.4%)
**Failed/skipped**: 615 (mostly download timeouts on genomes >3 GB)
**Total introns scanned**: 334,142,092
**Telotrons found (>=85% purity)**: 1,798
**Species with telotrons**: 144

---

## Top Species by Telotron Count

| # | Species | Group | Telotrons | Introns | Top 3' motif | Top 5' motif |
|---|---|---|---:|---:|---|---|
| 1 | *Rhipicephalus sanguineus* | Tick (Arachnida) | **660** | 237,840 | GGTT\|AG (73%) | GT\|TAGG |
| 2 | *Monocercomonoides exilis* | Metamonada (no mito) | **289** | 35,319 | GGTT\|AG (99%) | GT\|TAGG |
| 3 | *Eriocheir sinensis* | Crab (Crustacea) | **68** | 628,769 | GGTT\|AG (59%) | GT\|TAGG |
| 4 | *Eimeria necatrix* | Apicomplexa | **60** | 27,957 | GTTT\|AG (40%) | GT\|TTAG |
| 5 | *Spodoptera frugiperda* | Moth (Lepidoptera) | **55** | 325,622 | GGTT\|AG (82%) | GT\|TTGG |
| 6 | *Tympanuchus pallidicinctus* | Grouse (Aves) | **48** | 608,596 | GGTT\|AG (69%) | GT\|TAGG |
| 7 | *Pogona vitticeps* | Lizard (Reptilia) | **44** | 735,196 | AGGT\|AG (95%) | GT\|AGGT |
| 8 | *Takifugu rubripes* | Pufferfish | **44** | 656,232 | ATTC\|AG (23%) | GT\|ACTG |
| 9 | *Takifugu flavidus* | Pufferfish | **41** | 681,213 | GGTT\|AG (32%) | GT\|TAGG |
| 10 | *Cataglyphis hispanica* | Ant (Hymenoptera) | **36** | 192,831 | GGTT\|AG (100%) | GT\|TGGG |
| 11 | *Lagopus muta* | Ptarmigan (Aves) | **34** | 633,847 | GGTT\|AG (85%) | GT\|TAGG |
| 12 | *Eimeria mitis* | Apicomplexa | **28** | 33,471 | GTTT\|AG (36%) | GT\|TTAG |
| 13 | *Anser brachyrhynchus* | Goose (Aves) | **25** | 867,386 | GGTT\|AG (76%) | GT\|TAGG |
| 14 | *Leptidea sinapis* | Butterfly (Lepidoptera) | **24** | 205,167 | GGTT\|AG (79%) | GT\|TAGG |
| 15 | *Lucilia sericata* | Blowfly (Diptera) | **23** | 128,403 | TTTC\|AG (36%) | GT\|GAGT |
| 16 | *Botryllus schlosseri* | Tunicate (Tunicata) | **17** | 321,995 | AGGT\|AG (71%) | GT\|AGGT |

## Taxonomic Distribution

### High-confidence telotron lineages (>5 telotrons in >=1 species)

| Supergroup | Clade | Species (hits/total) | Total telotrons |
|---|---|---|---:|
| **Metamonada** | Oxymonadida | 1/1 | 289 |
| **Alveolata** | Apicomplexa (Eimeria) | 5/~15 | 122 |
| **Metazoa** | Arachnida (ticks) | 3/~10 | 663 |
| **Metazoa** | Crustacea | 2/~20 | 70 |
| **Metazoa** | Insecta (Lepidoptera) | 8/~50 | 113 |
| **Metazoa** | Insecta (Hymenoptera) | 5/~30 | 43 |
| **Metazoa** | Insecta (Diptera) | 2/~30 | 24 |
| **Metazoa** | Tunicata | 2/~5 | 24 |
| **Metazoa** | Actinopterygii (fish) | ~30/~200 | ~180 |
| **Metazoa** | Aves (birds) | 5/~100 | ~115 |
| **Metazoa** | Reptilia (lizards) | 3/~20 | ~46 |
| **Metazoa** | Cnidaria (corals) | 5/~10 | 15 |
| **Metazoa** | Syngnathiformes (pipefish) | 6/~8 | 19 |
| **Haptista** | Haptophyta | 1/1 | 3 |
| **Cryptista** | Cryptophyta | 1/1 | 5 |
| **Excavata** | Heterolobosea | 1/1 | 1 |
| **Amoebozoa** | Discosea | 1/1 | 2 |
| **Fungi** | Ascomycota | 3/~500 | 5 |
| **Archaeplastida** | Angiosperm | ~10/~200 | ~20 |

### Notable absences (0 telotrons)
- **Homo sapiens** (2.1M introns), **Mus musculus** (1.4M), **Equus caballus** (856k)
- **Drosophila melanogaster** and all other *Drosophila* spp.
- **Saccharomyces cerevisiae**, **Arabidopsis thaliana**
- **Plasmodium falciparum**, **Tetrahymena thermophila**, **Paramecium tetraurelia**
- All trypanosomatids (*T. brucei*, *Leishmania*)

## Boundary Motif Patterns

The dominant 3' acceptor motif (4bp before AG) across species:

| Motif | Interpretation | Species showing it |
|---|---|---|
| **GGTT\|AG** | Partial TTAGGG (template strand reads ...CCCTAA → GGTT before AG) | *Rhipicephalus* (73%), *Monocercomonoides* (99%), *Cataglyphis* (100%), *Lagopus* (85%), most TTAGGG organisms |
| **AGGT\|AG** | Partial TTAGGG (different reading frame) | *Pogona* (95%), *Botryllus* (71%), *Neurospora*, *Podospora* |
| **GTTT\|AG** | Partial TTTGGG (Eimeria/Paramecium-type repeat) | All *Eimeria* spp. (36-40%), *Naegleria*, *Malus*, *Diospyros* |
| **CCCT\|AG** | Partial CCCTAA (= TTAGGG complement) | *Emiliania huxleyi* only |

### Key observation
The **GGTT|AG** motif dominates in TTAGGG-repeat organisms — this is the template-strand reading of the telomeric repeat terminating at the splice acceptor. The Tara Oceans finding of CCCTAACC|AG (= GGTT on the opposite strand) is confirmed across diverse eukaryotes.

The *Eimeria* cluster consistently shows **GTTT|AG**, matching their TTTGGG telomeric repeat. This is the first cross-kingdom confirmation that the boundary motif is repeat-type-specific.

## Highlight Discoveries

### 1. Rhipicephalus sanguineus (brown dog tick) — 660 telotrons
The single largest telotron load in any reference genome. Ticks are known for massive genome expansion via repetitive elements; this suggests telomeric repeat invasion of introns may be a major contributor. Two other tick species (*Dermacentor silvarum*, *D. variabilis*) also show hits.

### 2. Monocercomonoides exilis — 289 telotrons
The only known eukaryote without mitochondria. 289 near-perfect TTTGGG-repeat introns (99.7% show GGTT|AG boundary). This organism lacks many standard DNA repair pathways — consistent with the Tara Oceans finding that repair-deficient lineages accumulate telotrons.

### 3. Eimeria spp. (5 species, 9-60 each) — 122 total
Apicomplexan parasites with TTTGGG repeats. Consistent GTTT|AG boundary matching their repeat type. Suggests telotron insertion is an active process in this lineage.

### 4. Vertebrate telotrons exist but are rare
Pufferfish (*Takifugu*, 44-41), birds (*Tympanuchus* 48, *Lagopus* 34, *Anser* 25), lizards (*Pogona* 44), and scattered fish species. These are genuine high-purity (>85%) telomeric repeat arrays in annotated introns — structurally distinct from the degraded ITSs previously reported in mammals.

### 5. Insect enrichment
Lepidoptera (8 species, 113 telotrons), Hymenoptera (5 species, 43), Diptera (2 species, 24), plus beetles and bugs. Insects use the 5-mer TTAGG repeat, but many of these hits match TTAGGG — possibly reflecting ancestral repeat contamination or recent repeat expansion.

### 6. Syngnathiformes (pipefish/seahorses) — 6 species, 19 telotrons
Disproportionate enrichment in this fish order, which is known for extreme genome restructuring (loss of immune genes, novel reproductive biology).

## Predictions Tested

| Prediction from plan.md | Result |
|---|---|
| CCCTAACC\|AG in TTAGGG organisms | **Confirmed** — GGTT\|AG (equivalent) is the dominant boundary in most TTAGGG species |
| Distinct signatures in variant-repeat organisms | **Confirmed** — *Eimeria* shows GTTT\|AG matching TTTGGG; *Naegleria* shows GTTT\|AG matching TTTGGG |
| Orientation specificity | Partially testable — most high-count species show template-strand bias |
| Human/horse ITSs show telotron boundary | **Not confirmed** — zero human/horse/mouse introns pass 85% purity |
| Plasmodium telotrons | **Not found** — only 1 intron at >10%, none at >85% |
| Tetrahymena telotrons | **Not found** — 0 hits despite 94k introns |

## Files

- [pan_euk_hits.tsv](pan_euk_hits.tsv) — All 144 species with telotrons
- [pan_euk_summary.tsv](pan_euk_summary.tsv) — All 2,312 genomes with status
- [telotrons/](../telotrons/) — Per-species telotron sequence files
- [checkpoint.json](../checkpoint.json) — Resumable pipeline state
