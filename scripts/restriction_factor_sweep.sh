#!/bin/bash
# Differential tblastn sweep of telomere-restriction / DSB-repair factors across
# the 14 apicomplexan genomes (4 outgroups + 5 telotron-poor sisters + 5 Eimeria).
# Goal: find a factor PRESENT in telotron-poor sisters but ABSENT in Eimeria.
# Rad50/Mre11 = positive controls (conserved; must hit everywhere).
set -uo pipefail
cd /scratch1/alex/telotrons
ENV_BIN=/home/alex/.snakemake-envs/telotrons/7aa0f2645427a08bdf054e56b1233123_/bin
export PATH="$ENV_BIN:$PATH"
EUTILS="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT=analysis/restriction_factors/sweep
mkdir -p "$OUT/seeds" "$OUT/db"

# genome_id  group  label
GENOMES=$(cat <<'EOF'
GCF_000002765.6 outgroup Plasmodium_falciparum
GCF_000003225.4 outgroup Theileria_annulata
GCF_000165395.2 outgroup Babesia_bovis
GCF_000165345.1 outgroup Cryptosporidium_parvum
GCF_000006565.2 poor Toxoplasma
GCF_000208865.1 poor Neospora
GCF_002563875.1 poor Besnoitia
GCF_002600585.1 poor Cystoisospora
GCF_002999335.1 poor Cyclospora
GCF_000499385.1 rich Eimeria_necatrix
GCF_000499425.1 rich Eimeria_acervulina
GCF_000499545.2 rich Eimeria_tenella
GCF_000499605.1 rich Eimeria_maxima
GCF_000499745.2 rich Eimeria_mitis
EOF
)

# Build (or reuse) a nucleotide blastdb per genome
while read gid grp label; do
  db="$OUT/db/$gid"
  if [ ! -f "$db.nin" ] && [ ! -f "$db.nsq" ]; then
    fa=$(ls raw/refseq/$gid/*_genomic.fna 2>/dev/null | head -1)
    [ -z "$fa" ] && { zcat raw/refseq/$gid/*_genomic.fna.gz > "$OUT/db/$gid.fna"; fa="$OUT/db/$gid.fna"; }
    makeblastdb -in "$fa" -dbtype nucl -out "$db" >/dev/null 2>&1
  fi
done <<< "$GENOMES"
echo "[dbs ready]"

# factor  ncbi_search_term
fetch_seeds () {
  local factor="$1" term="$2"
  local f="$OUT/seeds/${factor}.faa"
  [ -s "$f" ] && return
  local q=$(python3 -c "import urllib.parse,sys;print(urllib.parse.quote(sys.argv[1]))" "$term")
  local ids=$(curl -s "${EUTILS}/esearch.fcgi?db=protein&term=${q}&retmax=20" | grep -oP '(?<=<Id>)[0-9]+(?=</Id>)' | tr '\n' ',' | sed 's/,$//')
  [ -n "$ids" ] && curl -s "${EUTILS}/efetch.fcgi?db=protein&id=${ids}&rettype=fasta&retmode=text" -o "$f"
  echo "  seeds $factor: $(grep -c '>' "$f" 2>/dev/null || echo 0)"
}

ORG='(Saccharomyces cerevisiae[Organism] OR Homo sapiens[Organism] OR Trypanosoma brucei[Organism] OR Arabidopsis thaliana[Organism] OR Schizosaccharomyces pombe[Organism])'
echo "[fetching seeds]"
fetch_seeds Rad50_CTRL   "RAD50[Protein Name] AND $ORG NOT partial[All Fields]"
fetch_seeds Mre11_CTRL   "MRE11[Protein Name] AND $ORG NOT partial[All Fields]"
fetch_seeds Nbs1         "(nibrin[Protein Name] OR NBS1[Protein Name] OR XRS2[Protein Name]) AND $ORG"
fetch_seeds ATR          "(Serine/threonine-protein kinase ATR[Protein Name] OR MEC1[Protein Name] OR Rad3[Protein Name]) AND $ORG"
fetch_seeds Ku70         "(XRCC6[Protein Name] OR Ku70[Protein Name] OR HDF1[Protein Name]) AND $ORG"
fetch_seeds Ku80         "(XRCC5[Protein Name] OR Ku80[Protein Name] OR HDF2[Protein Name] OR YKU80[Protein Name]) AND $ORG"
fetch_seeds Stn1         "(STN1[Protein Name] OR CST complex subunit STN1[Protein Name]) AND $ORG"
fetch_seeds Ten1         "(TEN1[Protein Name] OR CST complex subunit TEN1[Protein Name]) AND $ORG"
fetch_seeds Cdc13_CTC1   "(CDC13[Protein Name] OR CTC1[Protein Name]) AND $ORG"
fetch_seeds Pot1         "(protection of telomeres protein 1[Protein Name] OR POT1[Protein Name]) AND $ORG"
fetch_seeds Rap1         "(RAP1[Protein Name] OR repressor activator protein 1[Protein Name]) AND $ORG"
fetch_seeds Pif1         "(PIF1[Protein Name] OR ATP-dependent DNA helicase PIF1[Protein Name]) AND $ORG"

echo "[tblastn sweep]  (best E per factor x genome)"
printf "%-12s" "factor" > "$OUT/matrix.tsv"
while read gid grp label; do printf "\t%s" "$label" >> "$OUT/matrix.tsv"; done <<< "$GENOMES"
printf "\n" >> "$OUT/matrix.tsv"

for sf in "$OUT"/seeds/*.faa; do
  [ -s "$sf" ] || continue
  factor=$(basename "$sf" .faa)
  printf "%-12s" "$factor" >> "$OUT/matrix.tsv"
  while read gid grp label; do
    best=$(tblastn -query "$sf" -db "$OUT/db/$gid" -evalue 1 -seg no -max_target_seqs 3 -num_threads 8 \
           -outfmt "6 evalue" 2>/dev/null | sort -g | head -1)
    printf "\t%s" "${best:-NA}" >> "$OUT/matrix.tsv"
  done <<< "$GENOMES"
  printf "\n" >> "$OUT/matrix.tsv"
  echo "  done $factor"
done
echo "[wrote $OUT/matrix.tsv]"
