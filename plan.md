# Pan-Eukaryotic Telotron Survey: NCBI Analysis Plan

## Objective

Extend the TERC boundary / CCCTAACC|AG analysis from Tara Oceans MAGs (606 MAGs, 34,343 telotrons) to **all annotated eukaryotic genomes** available through NCBI, covering model organisms and diverse lineages where telomeric-repeat introns have never been surveyed.

---

## 1. Install NCBI Datasets CLI

```bash
# On Linux (x86_64)
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
chmod +x datasets

# Or via conda
conda install -c conda-forge ncbi-datasets-cli
```

## 2. Build the Target Species List

### Strategy
Query NCBI for all eukaryotic RefSeq genomes with gene annotations (GFF3). We want breadth across the eukaryotic tree, not just model organisms.

```bash
# Get metadata for ALL annotated eukaryotic genomes
datasets summary genome taxon 2759 \
  --annotated \
  --assembly-source refseq \
  --assembly-level chromosome,complete \
  --as-json-lines | \
  dataformat tsv genome \
    --fields organism-name,assminfo-accession,assminfo-level,annotinfo-name,organism-tax-id,assmstats-total-sequence-len \
  > euk_genomes_annotated.tsv

# Also get scaffold-level assemblies (many protists are scaffold-level only)
datasets summary genome taxon 2759 \
  --annotated \
  --assembly-source refseq \
  --assembly-level scaffold \
  --as-json-lines | \
  dataformat tsv genome \
    --fields organism-name,assminfo-accession,assminfo-level,annotinfo-name,organism-tax-id,assmstats-total-sequence-len \
  >> euk_genomes_annotated.tsv
```

### Priority species (ensure these are included)

| Supergroup | Key species | Accession | Rationale |
|---|---|---|---|
| **Haptista** | *Emiliania huxleyi* CCMP1516 | GCF_000372725.1 | Closest reference to Phaeocystis; known repetitive genome |
| **Haptista** | *Phaeocystis globosa* | GCA_947623565.1 | Chromosome-scale; primary telotron lineage |
| **Haptista** | *Chrysochromulina tobin* | GCA_001275005.1 | Second haptophyte reference |
| **Stramenopiles** | *Thalassiosira pseudonana* | GCF_000149405.2 | Model diatom |
| **Stramenopiles** | *Phaeodactylum tricornutum* | GCF_000150955.2 | Model diatom |
| **Stramenopiles** | *Nannochloropsis gaditana* | GCF_000240725.1 | Eustigmatophyte |
| **Alveolata** | *Plasmodium falciparum* | GCF_000002765.6 | Known de novo telomere addition at DSBs |
| **Alveolata** | *Tetrahymena thermophila* | GCF_000189635.1 | Telomerase was discovered here |
| **Alveolata** | *Paramecium tetraurelia* | GCF_000165425.1 | Intron-rich ciliate |
| **Alveolata** | *Symbiodinium microadriaticum* | GCF_001939145.1 | Dinoflagellate (known introner activity) |
| **Archaeplastida** | *Arabidopsis thaliana* | GCF_000001735.4 | Plant model, TTTAGGG repeat |
| **Archaeplastida** | *Chlamydomonas reinhardtii* | GCF_000002595.2 | Green alga |
| **Archaeplastida** | *Micromonas pusilla* | GCF_000151265.2 | Prasinophyte (in Tara dataset) |
| **Archaeplastida** | *Ostreococcus tauri* | GCF_000214015.3 | Smallest free-living eukaryote |
| **Opisthokonta** | *Homo sapiens* | GCF_000001405.40 | 31% of ITSs in introns |
| **Opisthokonta** | *Mus musculus* | GCF_000001635.27 | Known ITS insertions |
| **Opisthokonta** | *Equus caballus* | GCF_002863925.1 | Known ITS insertions (Nergadze 2004) |
| **Opisthokonta** | *Saccharomyces cerevisiae* | GCF_000146045.2 | Pif1/Ku70 genetics model |
| **Opisthokonta** | *Neurospora crassa* | GCF_000182925.2 | Fungal model |
| **Amoebozoa** | *Dictyostelium discoideum* | GCF_000004695.1 | AG-rich telomeres |
| **Excavata** | *Trypanosoma brucei* | GCF_000002445.2 | TTAGGG repeat |
| **Excavata** | *Leishmania major* | GCF_000002725.2 | Kinetoplastid |
| **Rhizaria** | *Bigelowiella natans* | GCA_000002455.1 | Chlorarachniophyte |

### Broader sweep
```bash
# Get ALL annotated eukaryotic genomes (expect ~2,000+ species as of 2025)
datasets summary genome taxon 2759 --annotated --assembly-source refseq \
  --as-json-lines > all_euk_refseq.jsonl

# Parse to get accessions
cat all_euk_refseq.jsonl | \
  python3 -c "
import json, sys
for line in sys.stdin:
    r = json.loads(line)
    for a in r.get('reports', [r]):
        org = a.get('organism', {})
        ai = a.get('assembly_info', {})
        print(f\"{ai.get('assembly_accession','')}\t{org.get('organism_name','')}\t{ai.get('assembly_level','')}\")
" > all_euk_accessions.tsv
```

## 3. Download Genomes + GFF3 Annotations

```bash
# Download in batches of ~50 genomes at a time
mkdir -p genomes

# For each accession:
while read acc name level; do
  datasets download genome accession "$acc" \
    --include genome,gff3 \
    --filename "genomes/${acc}.zip"
  unzip -o "genomes/${acc}.zip" -d "genomes/${acc}"
done < priority_accessions.tsv

# Or bulk download with dehydration for large sets:
datasets download genome taxon 2759 \
  --annotated --assembly-source refseq \
  --include genome,gff3 \
  --dehydrated \
  --filename all_euk_dehydrated.zip

unzip all_euk_dehydrated.zip -d all_euk/
datasets rehydrate --directory all_euk/
```

## 4. Extract Introns from GFF3

For each genome, extract all annotated introns (gaps between exons within the same mRNA):

```python
#!/usr/bin/env python3
"""extract_introns.py — Extract intron sequences from genome + GFF3."""
import re, sys, os
from collections import defaultdict

def load_fasta(path):
    seqs = {}; curr = None; parts = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if curr: seqs[curr] = ''.join(parts)
                curr = line[1:].strip().split()[0]; parts = []
            else: parts.append(line.strip().upper())
    if curr: seqs[curr] = ''.join(parts)
    return seqs

def extract_introns(gff_path, fasta_path, out_tsv):
    genome = load_fasta(fasta_path)

    # Parse GFF: collect exons per mRNA
    mrna_exons = defaultdict(list)
    mrna_strand = {}
    mrna_contig = {}

    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            if parts[2] not in ('exon', 'CDS'): continue

            attrs = dict(re.findall(r'(\w+)=([^;]+)', parts[8]))
            parent = attrs.get('Parent', '')
            contig = parts[0]
            start = int(parts[3]) - 1  # 0-based
            end = int(parts[4])
            strand = parts[6]

            mrna_exons[parent].append((start, end))
            mrna_strand[parent] = strand
            mrna_contig[parent] = contig

    with open(out_tsv, 'w') as out:
        out.write('genome_id\tcontig\tstart\tend\tlength\tstrand\tmrna_id\tdonor\tacceptor\tintron_seq\n')

        for mrna, exons in mrna_exons.items():
            exons.sort()
            contig = mrna_contig[mrna]
            strand = mrna_strand[mrna]
            if contig not in genome: continue
            seq = genome[contig]

            for i in range(len(exons) - 1):
                istart = exons[i][1]
                iend = exons[i+1][0]
                if iend <= istart: continue

                intron_seq = seq[istart:iend]
                if len(intron_seq) < 20: continue

                donor = intron_seq[:2]
                acceptor = intron_seq[-2:]

                genome_id = os.path.basename(fasta_path).replace('.fna', '').replace('.fa', '')
                out.write(f'{genome_id}\t{contig}\t{istart}\t{iend}\t{iend-istart}\t{strand}\t{mrna}\t{donor}\t{acceptor}\t{intron_seq}\n')

if __name__ == '__main__':
    extract_introns(sys.argv[1], sys.argv[2], sys.argv[3])
```

```bash
# Run on all genomes
for acc_dir in genomes/GCF_*/; do
    acc=$(basename "$acc_dir")
    gff=$(find "$acc_dir" -name "*.gff" | head -1)
    fna=$(find "$acc_dir" -name "*.fna" | head -1)
    [ -z "$gff" ] || [ -z "$fna" ] && continue
    python3 extract_introns.py "$gff" "$fna" "introns/${acc}_introns.tsv"
done
```

## 5. Screen Introns for Telomeric Content

```python
#!/usr/bin/env python3
"""screen_telotrons.py — Classify introns by telomeric hexamer coverage."""
import re, sys, os
from collections import Counter

# All 12 rotations of TTAGGG and CCCTAA
TELO_HEX = set()
for base in ['TTAGGG', 'CCCTAA']:
    for i in range(6):
        TELO_HEX.add(base[i:] + base[:i])

# Also check common variant repeats
VARIANT_REPEATS = {
    'TTTAGGG': 'plant-type (7mer)',      # Archaeplastida
    'TTAGG':   'insect-type (5mer)',      # Insects
    'TTTTGGGG': 'Tetrahymena-type (8mer)', # Some ciliates
    'TTTTAGGG': 'Paramecium-type (8mer)',  # Some ciliates
}

def telo_coverage(seq, hexset=TELO_HEX):
    """Fraction of bases covered by any hexamer in set."""
    n = len(seq)
    if n == 0: return 0.0
    covered = [False] * n
    for h in hexset:
        for m in re.finditer(h, seq):
            for i in range(m.start(), m.end()):
                covered[i] = True
    return sum(covered) / n

def strand_ratio(seq):
    """TTAGGG-rotation vs CCCTAA-rotation counts."""
    fwd = sum(len(re.findall(r, seq)) for r in ['TTAGGG','TAGGGT','AGGGTT','GGGTTA','GGTTAG','GTTAGG'])
    rev = sum(len(re.findall(r, seq)) for r in ['CCCTAA','CCTAAC','CTAACC','TAACCC','AACCCT','ACCCTA'])
    total = fwd + rev
    if total == 0: return 0.5, 'unknown'
    frac = fwd / total
    if frac > 0.7: return frac, 'coding-strand'
    elif frac < 0.3: return frac, 'template-strand'
    else: return frac, 'converging'

def screen_file(intron_tsv, out_tsv):
    with open(intron_tsv) as f, open(out_tsv, 'w') as out:
        header = f.readline().strip()
        out.write(header + '\ttelo_coverage\tfwd_frac\torientation\tdonor_8\tacceptor_8\n')

        for line in f:
            parts = line.strip().split('\t')
            seq = parts[-1].upper()
            cov = telo_coverage(seq)

            if cov < 0.10: continue  # Skip non-telomeric

            frac, orient = strand_ratio(seq)

            # Normalize to GT-AG if needed
            d2 = seq[:2]; a2 = seq[-2:]
            if d2 == 'CT' and a2 == 'AC':
                seq = seq[::-1].translate(str.maketrans('ACGT','TGCA'))
                d2 = seq[:2]; a2 = seq[-2:]

            donor_8 = seq[:8] if len(seq) >= 8 else seq[:len(seq)]
            acceptor_8 = seq[-8:] if len(seq) >= 8 else seq[-len(seq):]

            out.write(line.strip() + f'\t{cov:.4f}\t{frac:.4f}\t{orient}\t{donor_8}\t{acceptor_8}\n')

if __name__ == '__main__':
    screen_file(sys.argv[1], sys.argv[2])
```

## 6. TERC Boundary Analysis

```python
#!/usr/bin/env python3
"""terc_boundary_analysis.py — Analyze splice junction motifs across all genomes."""
import sys, os, re
from collections import Counter, defaultdict
import numpy as np

def rc(seq):
    return seq[::-1].translate(str.maketrans('ACGT','TGCA'))

def analyze_boundaries(telotron_files, out_prefix):
    """Analyze CCCTAACC|AG and GT|TAGG signatures across all genomes."""

    # Collect motifs per genome and per orientation
    results = defaultdict(lambda: {
        'n_total': 0, 'n_gt_ag': 0,
        'donor_4': Counter(), 'acceptor_4': Counter(),
        'by_orient': defaultdict(lambda: {'n': 0, 'acceptor_4': Counter()})
    })

    for tfile in telotron_files:
        genome_id = os.path.basename(tfile).replace('_telotrons.tsv', '')

        with open(tfile) as f:
            header = f.readline()
            cols = header.strip().split('\t')
            ci = {c: i for i, c in enumerate(cols)}

            for line in f:
                parts = line.strip().split('\t')
                cov = float(parts[ci['telo_coverage']])
                if cov < 0.50: continue

                seq = parts[ci['intron_seq']].upper()
                orient = parts[ci['orientation']]

                results[genome_id]['n_total'] += 1

                # Normalize to GT-AG
                d2 = seq[:2]; a2 = seq[-2:]
                if d2 == 'CT' and a2 == 'AC':
                    seq = rc(seq)
                    d2 = seq[:2]; a2 = seq[-2:]

                if d2 != 'GT' or a2 != 'AG': continue
                results[genome_id]['n_gt_ag'] += 1

                # 4bp after GT, 4bp before AG
                if len(seq) >= 8:
                    d4 = seq[2:6]
                    a4 = seq[-6:-2]
                    results[genome_id]['donor_4'][d4] += 1
                    results[genome_id]['acceptor_4'][a4] += 1
                    results[genome_id]['by_orient'][orient]['n'] += 1
                    results[genome_id]['by_orient'][orient]['acceptor_4'][a4] += 1

    # Write summary
    with open(f'{out_prefix}_summary.tsv', 'w') as out:
        out.write('genome_id\tn_total\tn_gt_ag\t'
                  'top_donor_motif\tdonor_pct\t'
                  'top_acceptor_motif\tacceptor_pct\t'
                  'AACC_AG_pct\t'
                  'AACC_template_pct\tAACC_converging_pct\tAACC_coding_pct\n')

        for gid, data in sorted(results.items()):
            if data['n_gt_ag'] < 10: continue

            top_d = data['donor_4'].most_common(1)
            top_a = data['acceptor_4'].most_common(1)
            total = sum(data['acceptor_4'].values())
            aacc_pct = 100 * data['acceptor_4'].get('AACC', 0) / total if total > 0 else 0

            # Per-orientation AACC%
            orient_pcts = {}
            for o in ['template-strand', 'converging', 'coding-strand']:
                od = data['by_orient'][o]
                ot = sum(od['acceptor_4'].values())
                orient_pcts[o] = 100 * od['acceptor_4'].get('AACC', 0) / ot if ot > 0 else 0

            out.write(f"{gid}\t{data['n_total']}\t{data['n_gt_ag']}\t"
                      f"{top_d[0][0] if top_d else 'NA'}\t"
                      f"{100*top_d[0][1]/total:.1f}%\t" if top_d and total > 0 else "NA\t"
                      f"{top_a[0][0] if top_a else 'NA'}\t"
                      f"{100*top_a[0][1]/total:.1f}%\t" if top_a and total > 0 else "NA\t"
                      f"{aacc_pct:.1f}%\t"
                      f"{orient_pcts.get('template-strand',0):.1f}%\t"
                      f"{orient_pcts.get('converging',0):.1f}%\t"
                      f"{orient_pcts.get('coding-strand',0):.1f}%\n")

if __name__ == '__main__':
    import glob
    files = glob.glob(sys.argv[1])
    analyze_boundaries(files, sys.argv[2])
```

## 7. Run the Full Pipeline

```bash
#!/bin/bash
# run_pan_eukaryotic_telotron_survey.sh

set -euo pipefail

OUTDIR="pan_euk_telotrons"
mkdir -p "$OUTDIR"/{genomes,introns,telotrons,results}

# Step 1: Download genomes
echo "=== Downloading genomes ==="
while IFS=$'\t' read -r acc name; do
    echo "  Downloading $name ($acc)..."
    datasets download genome accession "$acc" \
        --include genome,gff3 \
        --filename "$OUTDIR/genomes/${acc}.zip" 2>/dev/null || continue
    unzip -qo "$OUTDIR/genomes/${acc}.zip" -d "$OUTDIR/genomes/${acc}" 2>/dev/null
done < priority_accessions.tsv

# Step 2: Extract introns
echo "=== Extracting introns ==="
for acc_dir in "$OUTDIR"/genomes/GC*/; do
    acc=$(basename "$acc_dir")
    gff=$(find "$acc_dir" -name "*.gff" 2>/dev/null | head -1)
    fna=$(find "$acc_dir" -name "*.fna" 2>/dev/null | head -1)
    [ -z "$gff" ] || [ -z "$fna" ] && continue
    echo "  $acc..."
    python3 extract_introns.py "$gff" "$fna" "$OUTDIR/introns/${acc}_introns.tsv"
done

# Step 3: Screen for telomeric content
echo "=== Screening for telotrons ==="
for intron_file in "$OUTDIR"/introns/*_introns.tsv; do
    acc=$(basename "$intron_file" _introns.tsv)
    echo "  $acc..."
    python3 screen_telotrons.py "$intron_file" "$OUTDIR/telotrons/${acc}_telotrons.tsv"
done

# Step 4: TERC boundary analysis
echo "=== TERC boundary analysis ==="
python3 terc_boundary_analysis.py "$OUTDIR/telotrons/*_telotrons.tsv" "$OUTDIR/results/pan_euk"

echo "=== Done ==="
echo "Results in $OUTDIR/results/"
```

## 8. Key Predictions to Test

### The CCCTAACC|AG signature should be:
1. **Present in any TTAGGG-repeat organism** — Haptophyta, vertebrates (human, horse, mouse), trypanosomes, any species with TTAGGG telomeres and annotated introns containing telomeric repeats
2. **Absent in organisms with different repeat units** — *Arabidopsis* (TTTAGGG), *Tetrahymena* (TTGGGG), *Paramecium* (TTTGGG) should show analogous but distinct boundary signatures matching THEIR TERC templates
3. **Orientation-specific** — template-only introns should show the highest AACC|AG frequency; coding-strand introns should show 0%
4. **Enriched at 2nt overshoot** — the partial repeat should consistently be 2nt, not 1 or 3, across all species

### Expected novel findings across eukaryotes:
- **Human/horse**: Nergadze (2004, 2007, 2021) documented ITS at DSB sites in primates/mammals. 31% of human ITSs are intronic. These should show the CCCTAACC|AG signature if they were telomerase-derived.
- **Plasmodium**: Known de novo telomere addition at subtelomeric DSBs (Calhoun 2017). Check for telomeric introns.
- **Tetrahymena**: Telomerase discovered here. TTGGGG repeat → expect analogous boundary signature with 2nt overshoot matching Tetrahymena TERC.
- **Trypanosoma**: TTAGGG repeats → direct comparison to haptophyte telotrons.

## 9. Variant Repeat Analysis

For organisms with non-TTAGGG telomeres, build variant hexamer sets:

| Organism | Repeat | Hexamer set | Expected 3' motif |
|---|---|---|---|
| Plants | TTTAGGG | 14 rotations (7-mer) | TTTAG + 2nt overshoot + AG |
| *Tetrahymena* | TTGGGG | 12 rotations (6-mer) | Depends on TERC template |
| Insects | TTAGG | 10 rotations (5-mer) | Depends on TERC template |
| *Dictyostelium* | AG(1-8) | Variable | Variable |

## 10. Expected Output

### Table: Pan-eukaryotic telotron survey
| Species | Supergroup | Repeat | Total introns | Telo introns (≥50%) | AACC\|AG % | Template-only % | Converging % |
|---|---|---|---|---|---|---|---|

### Figure: Cross-species TERC boundary conservation
- Panel a: Phylogenetic tree with telotron counts per species
- Panel b: 3' acceptor motif frequency across species (heatmap)
- Panel c: Orientation specificity (template vs coding) across species
- Panel d: 2nt overshoot consistency across repeat types

---

## Estimated Resources

- **Storage**: ~500 GB for all annotated eukaryotic genomes + GFF3
- **Compute**: ~4-8 hours on 16-core machine for full intron extraction + screening
- **Priority run** (23 species above): ~20 GB, ~30 minutes

## Quick Start (Priority Species Only)

```bash
# 1. Create accession list
cat > priority_accessions.tsv << 'EOF'
GCF_000372725.1	Emiliania_huxleyi
GCF_000149405.2	Thalassiosira_pseudonana
GCF_000150955.2	Phaeodactylum_tricornutum
GCF_000002765.6	Plasmodium_falciparum
GCF_000189635.1	Tetrahymena_thermophila
GCF_000001735.4	Arabidopsis_thaliana
GCF_000002595.2	Chlamydomonas_reinhardtii
GCF_000001405.40	Homo_sapiens
GCF_000001635.27	Mus_musculus
GCF_002863925.1	Equus_caballus
GCF_000146045.2	Saccharomyces_cerevisiae
GCF_000002445.2	Trypanosoma_brucei
GCF_000004695.1	Dictyostelium_discoideum
GCF_000165425.1	Paramecium_tetraurelia
GCF_001939145.1	Symbiodinium_microadriaticum
GCF_000240725.1	Nannochloropsis_gaditana
GCF_000151265.2	Micromonas_pusilla
GCF_000214015.3	Ostreococcus_tauri
GCF_000182925.2	Neurospora_crassa
GCF_000002725.2	Leishmania_major
EOF

# 2. Run
bash run_pan_eukaryotic_telotron_survey.sh
```
