# Telotron survey pipeline: scan eukaryotic introns (NCBI RefSeq + Tara SMAGs)
# for telomeric-repeat arrays, filter, analyze, plot, and package outputs.


import os

configfile: "config.yaml"


TARA_BASE = config["tara_base"]
REFSEQ_URL = config["refseq_url"]
THREADS = int(config["threads"])
ENV = config["env"]
REFSEQ_GROUPS = "|".join(config["refseq_groups"])
ACCESSIONS = config.get("accessions") or []
CANONICAL = config.get("canonical_telomere_motifs", {}) or {}
CANONICAL_BY_GENOME = CANONICAL.get("by_genome", {}) or {}
CANONICAL_BY_GROUP = CANONICAL.get("by_group", {}) or {}

# CLI accepts either a comma-separated string or a YAML list.
_motifs = config["telomere_motifs"]
TELOMERE_MOTIFS = ",".join(_motifs) if isinstance(_motifs, list) else _motifs

SCAN_MIN_FRAC = config["scan"]["min_repeat_frac"]
SCAN_MAX_FLANK_FRAC = config["scan"]["max_flank_repeat_frac"]
SCAN_MIN_INTRON_LEN = config["scan"].get("min_intron_len", 30)
FILTER_MIN_FRAC = config["filter"]["min_repeat_frac"]
FILTER_BIDIR_FRAC = config["filter"]["bidir_min_repeat_frac"]
FILTER_BIDIR_HITS = config["filter"]["bidir_min_hits"]
# Consolidated "what counts as a telotron" thresholds — see top of config.yaml
# for the full threshold map. These keys default to the script-internal values
# when the config block is missing (backward-compatible).
ARCH_MIN_ARRAY = int(config.get("architecture", {}).get("min_array_bp", 18))
INTERSTITIAL_MIN_ARRAY_LEN = int(config.get("interstitial", {}).get("min_array_len", 30))
INTERSTITIAL_MIN_ITS_UNITS = int(config.get("interstitial", {}).get("min_its_units", 4))
FILTER_FLAGS = " ".join(
    flag
    for flag, on in [
        ("--collapse-unique-loci", config["filter"]["collapse_unique_loci"]),
        ("--require-terminal-motif-match", config["filter"]["require_terminal_motif_match"]),
        ("--require-canonical-splice", config["filter"]["require_canonical_splice"]),
    ]
    if on
)

FIGURES = ["work/results/figures/telotron_counts.png"]
MEME_ENV = "envs/meme.yaml"

# Sentinels (path = ".../<name>/.done" or the canonical output file).
EXTRACT_FASTA_SENTINEL          = "work/results/telotron_fasta/.done"
CTRL_FLANKED_SENTINEL           = "work/results/non_telotron_fasta/.done"
TERT_DEEP_HOMOLOGY_TSV          = "work/results/tert_deep_homology/confirmed_tert.tsv"
TELOTRON_ORTHO_SENTINEL         = "work/results/telotron_orthologs/.done"
TELOTRON_ORTHO_LOCI_PDF         = "work/results/figures/telotron_ortholog_loci.pdf"
TELOTRON_ORTHO_V2               = "work/results/telotron_orthologs_v2"
TELOTRON_ORTHO_TEXT             = TELOTRON_ORTHO_V2 + "/locus_text/.textdump.done"

# Per-locus 5'/3' end orientation. Hybrid (GC, CG) = G-rich one end + C-rich
# the other (strand switch within an intron). Interstitial arrays are
# single-direction (GG or CC only). Non-telotron control introns aren't stratified.
TELOTRON_CAT = {
    "GG": ("architecture", "GT-F-AG"),
    "CC": ("architecture", "GT-R-AG"),
    "GC": ("architecture", "GT-F-R-AG,GT-F-linker-R-AG"),
    "CG": ("architecture", "GT-R-linker-F-AG"),
}
INTERSTITIAL_CAT = {c: ("cat5_3", c) for c in ("GG", "GC", "CG", "CC")}


# Default target = core survey through `package` + `pipeline_report`. Named
# subgraph targets: `telomerase_search`, `orthologs`, `architecture_analyses`,
# `analysis_arm`, `telogator2`. Kick off individual rules by name for reruns.
rule all:
    input:
        "work/results/telotron_pipeline_outputs.zip",
        "work/results/pipeline_report.html",


rule telomerase_search:
    input:
        TERT_DEEP_HOMOLOGY_TSV,


# Telotron-host-gene ortholog alignment (protein-space MSA + telotron-vs-intron DNA)
# + the compiled per-locus figure PDF.
rule orthologs:
    input:
        TELOTRON_ORTHO_SENTINEL,
        TELOTRON_ORTHO_LOCI_PDF,
        TELOTRON_ORTHO_TEXT,


# Build the unified genome manifest from RefSeq assembly summary + Tara SMAGs index.
# Keeps only annotated euk lineages (group, col 25) with a non-"na" FTP path (col 20).
# If config["accessions"] is non-empty, subset to just those genome_ids.
rule manifests:
    output:
        refseq="work/manifests/refseq_euk.tsv",
        tara="work/manifests/tara_mags.tsv",
        all="work/manifests/all_genomes.tsv",
    params:
        accessions=ACCESSIONS,
    shell:
        r"""
        mkdir -p manifests
        curl -L --fail -s {REFSEQ_URL} \
          | awk -F'\t' 'BEGIN{{OFS="\t"}} NR==1 || /^#/ {{next}} \
              $25 ~ /^({REFSEQ_GROUPS})$/ && $20!="na" \
              {{print $1,$8,$25,$20}}' \
          > {output.refseq}

        curl -L --fail -s {TARA_BASE}/SMAGs_v1_individual.gff.tar.gz \
          | tar -tzf - \
          | sed -n 's#.*/\(TARA_.*MAG_[0-9][0-9]*\)\.gmove\.gff#\1#p' \
          | sort -u \
          | awk 'BEGIN{{OFS="\t"}} {{print $1,"Tara Oceans MAG","tara","genoscope"}}' \
          > {output.tara}

        printf "genome_id\torganism\tgroup\tsource\n" > {output.all}
        cat {output.refseq} {output.tara} >> {output.all}

        # Optional accession whitelist: keep header + matching genome_ids.

        # (subset is applied below; canonical_motifs rule consumes the final {output.all})
        if [ -n "{params.accessions}" ]; then
            printf '%s\n' {params.accessions} > work/manifests/.accessions.txt
            awk 'NR==FNR{{keep[$1];next}} FNR==1 || $1 in keep' \
                work/manifests/.accessions.txt {output.all} > {output.all}.subset
            mv {output.all}.subset {output.all}
            awk 'NR==FNR{{keep[$1];next}} $1 in keep' \
                work/manifests/.accessions.txt {output.refseq} > {output.refseq}.subset
            mv {output.refseq}.subset {output.refseq}
            awk 'NR==FNR{{keep[$1];next}} $1 in keep' \
                work/manifests/.accessions.txt {output.tara} > {output.tara}.subset
            mv {output.tara}.subset {output.tara}
        fi
        """


# Tara SMAGs ship as two monolithic ~50 GB tarballs (contigs + GFF). Stream-extract only
# the MAGs listed in work/manifests/tara_mags.tsv so the test set doesn't blow up disk.
rule tara_archives:
    input:
        "work/manifests/tara_mags.tsv",
    output:
        fna="data/raw/tara/.fna.done",
        gff="data/raw/tara/.gff.done",
    shell:
        r"""
        mkdir -p data/raw/tara
        mapfile -t MAGS < <(awk -F'\t' '{{print $1}}' {input})
        if [ "${{#MAGS[@]}}" -eq 0 ]; then
            touch {output.fna} {output.gff}; exit 0
        fi
        FNA_PATS=(); GFF_PATS=()
        for m in "${{MAGS[@]}}"; do
            FNA_PATS+=("Contigs/${{m}}.*")
            GFF_PATS+=("Genes/GFF/${{m}}.gmove.*")
        done
        # Stream both tarballs concurrently — independent, no shared output paths.
        ( curl -L --fail -s {TARA_BASE}/SMAGs_contigs_individual.fna.tar.gz \
            | tar --wildcards -xzf - -C data/raw/tara "${{FNA_PATS[@]}}" ) &
        ( curl -L --fail -s {TARA_BASE}/SMAGs_v1_individual.gff.tar.gz \
            | tar --wildcards -xzf - -C data/raw/tara "${{GFF_PATS[@]}}" ) &
        wait
        touch {output.fna} {output.gff}
        """


# Derive per-genome FNA and GFF URLs from the RefSeq FTP path.
rule refseq_urls:
    input:
        "work/manifests/refseq_euk.tsv",
    output:
        "work/manifests/refseq_urls.tsv",
    shell:
        r"""
        awk -F'\t' 'BEGIN{{OFS="\t"}} {{
          sub(/\/$/, "", $4);
          n=split($4,a,"/"); base=a[n];
          print $1, $4"/"base"_genomic.fna.gz", $4"/"base"_genomic.gff.gz"
        }}' {input} > {output}
        """


# Fan out per-genome curls. NCBI rate-limits aggressive parallel pulls, so cap the
# inner pool at 8 regardless of THREADS. curl retries handle transient 503s.
# `exit 255` aborts xargs on a persistent failure after retries.
rule download_refseq:
    input:
        "work/manifests/refseq_urls.tsv",
    output:
        touch("data/raw/refseq/.done"),
    threads: 8
    shell:
        r"""
        mkdir -p data/raw/refseq
        awk -F'\t' '{{print $1"\t"$2"\n"$1"\t"$3}}' {input} \
          | xargs -P 8 -n2 sh -c \
              'mkdir -p data/raw/refseq/$0 && curl -L --fail --retry 5 --retry-delay 3 --retry-all-errors -s "$1" -o data/raw/refseq/$0/$(basename "$1") || exit 255'
        touch {output}
        """


# Build a `genome_id → motif` TSV from the curated literature mapping in config.yaml.
# Scan stage uses this in preference to contig-end scanning.
rule canonical_motifs:
    input:
        manifest="work/manifests/all_genomes.tsv",
    output:
        "work/manifests/canonical_motifs.tsv",
    params:
        by_genome=CANONICAL_BY_GENOME,
        by_group=CANONICAL_BY_GROUP,
    run:
        import csv
        with open(input.manifest) as f, open(output[0], "w", newline="") as out:
            w = csv.writer(out, delimiter="\t")
            w.writerow(["genome_id", "motif"])
            for r in csv.DictReader(f, delimiter="\t"):
                gid = r["genome_id"]
                grp = r["group"]
                motif = params.by_genome.get(gid)
                if motif is None:
                    motif = params.by_group.get(grp, "")
                w.writerow([gid, motif or ""])


# Core scan: per-intron motif coverage, splice signals, orientation,
# distance-to-end. Broad sieve (low --min-repeat-frac) so linker-architecture
# loci survive into the TSV; strict cutoffs live in filter_final.
# --max-flank-repeat-frac 0.25 + --min-intron-len 30 reject misannotated subtelomere
# introns (flanks themselves telomeric) and 1-bp degenerate GFF intron annotations.
rule scan_all:
    input:
        manifest="work/manifests/all_genomes.tsv",
        canonical="work/manifests/canonical_motifs.tsv",
        tara=["data/raw/tara/.fna.done", "data/raw/tara/.gff.done"],
        refseq="data/raw/refseq/.done",
    output:
        loci="work/results/all_telotron_loci.tsv",
        introns="work/results/all_introns_scanned.tsv",
        summary="work/results/all_species_raw_summary.tsv",
    threads: THREADS
    conda:
        ENV
    shell:
        r"""
        mkdir -p results
        python scripts/scan_telotrons.py \
            --manifest {input.manifest} \
            --canonical-motifs {input.canonical} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --motifs {TELOMERE_MOTIFS} \
            --min-repeat-frac {SCAN_MIN_FRAC} --max-flank-repeat-frac {SCAN_MAX_FLANK_FRAC} \
            --min-intron-len {SCAN_MIN_INTRON_LEN} \
            --threads {threads} \
            --loci {output.loci} --introns {output.introns} --summary {output.summary}
        """


# Tighten candidates → final set. --require-terminal-motif-match enforces that
# the intronic motif matches the genome's actual telomere motif (kills cross-motif noise).
rule filter_final:
    input:
        loci="work/results/all_telotron_loci.tsv",
        introns="work/results/all_introns_scanned.tsv",
        summary="work/results/all_species_raw_summary.tsv",
    output:
        final="work/results/final_telotron_set.tsv",
        species="work/results/final_species_summary.tsv",
        neg="work/results/final_negative_controls.tsv",
    conda:
        ENV
    shell:
        r"""
        python scripts/filter_final_set.py \
            --loci {input.loci} --introns {input.introns} --summary {input.summary} \
            --min-repeat-frac {FILTER_MIN_FRAC} \
            --bidir-min-repeat-frac {FILTER_BIDIR_FRAC} \
            --bidir-min-hits {FILTER_BIDIR_HITS} \
            {FILTER_FLAGS} \
            --final {output.final} --species {output.species} --negatives {output.neg}
        """


# Collapse doubly-assembled exons / close gene copies by aligning per-locus
# 250 bp upstream + 250 bp downstream flanks within each species (all-vs-all
# blastn). Two loci are duplicates if either flank shares >=100 bp at >=90%
# identity; the longest-intron member of each component is kept.
rule dedup_telotrons:
    input:
        final="work/results/final_telotron_set.tsv",
        tara=["data/raw/tara/.fna.done"],
        refseq="data/raw/refseq/.done",
    output:
        final="work/results/final_telotron_set_dedup.tsv",
        log="work/results/dedup_log.tsv",
    threads: THREADS
    conda:
        ENV
    shell:
        r"""
        python scripts/dedup_telotrons.py \
            --final {input.final} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --out-final {output.final} --out-log {output.log} \
            --threads {threads}
        """


# Classify each (deduplicated) telotron locus into one of:
# GT-F-AG, GT-R-AG, GT-F-R-AG, GT-F-linker-R-AG, GT-R-linker-F-AG, Other.
# Emits the long-form boundary-kmer-by-architecture TSV used by the per-species
# figures (donor / acceptor / linker_left / linker_right 6-mers).
rule classify_architecture:
    input:
        final="work/results/final_telotron_set_dedup.tsv",
        tara=["data/raw/tara/.fna.done"],
        refseq="data/raw/refseq/.done",
    output:
        loci="work/results/final_telotron_set_architecture.tsv",
        kmers="work/results/boundary_kmers_by_architecture.tsv",
    conda:
        ENV
    shell:
        r"""
        python scripts/classify_telotron_architecture.py \
            --final {input.final} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --min-array {ARCH_MIN_ARRAY} \
            --out-loci {output.loci} --out-kmers {output.kmers}
        """


# Boundary k-mer enrichment, length-matched distance-to-end test, architecture summary.
rule analyze:
    input:
        final="work/results/final_telotron_set_dedup.tsv",
        introns="work/results/all_introns_scanned.tsv",
        species="work/results/final_species_summary.tsv",
        neg="work/results/final_negative_controls.tsv",
    output:
        kmers="work/results/boundary_kmer_enrichment.tsv",
        dist="work/results/distance_to_end.tsv",
        arch="work/results/architecture_summary.tsv",
    threads: THREADS
    conda:
        ENV
    shell:
        r"""
        python scripts/analyze_telotrons.py \
            --final {input.final} --introns {input.introns} --species {input.species} \
            --boundary-kmers {output.kmers} --distance {output.dist} --architecture {output.arch} \
            --threads {threads}
        """


# Four headline figures (counts, orientation, boundary k-mers, distance scatter).
rule figures:
    input:
        species="work/results/final_species_summary.tsv",
        final="work/results/final_telotron_set.tsv",
        kmers="work/results/boundary_kmer_enrichment.tsv",
        dist="work/results/distance_to_end.tsv",
        arch="work/results/architecture_summary.tsv",
    output:
        FIGURES,
    conda:
        ENV
    shell:
        r"""
        mkdir -p work/results/figures
        python scripts/plot_telotrons.py \
            --species {input.species} --final {input.final} \
            --kmers {input.kmers} --distance {input.dist} --architecture {input.arch} \
            --outdir work/results/figures
        """


# Per-species, per-architecture FASTA + flanked-text extracts.
# Flanked lines: [LEFT100] [INTRON] [RIGHT100], or for linker archs
# [LEFT100] [ARRAY1] [LINKER] [ARRAY2] [RIGHT100]. One subdir per species,
# one file per architecture.
rule extract_telotron_fasta:
    input:
        arch="work/results/final_telotron_set_architecture.tsv",
        tara=["data/raw/tara/.fna.done"],
        refseq="data/raw/refseq/.done",
    output:
        touch(EXTRACT_FASTA_SENTINEL),
    conda:
        ENV
    shell:
        r"""
        mkdir -p work/results/telotron_fasta work/results/telotron_flanked
        python scripts/extract_telotron_fasta.py \
            --final {input.arch} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --fasta-dir work/results/telotron_fasta \
            --flanked-dir work/results/telotron_flanked
        touch {output}
        """


# Per-genome ORF mask: 6-frame ATG→stop scan, ORFs >= min_orf_nt (default 450,
# at which shuffled E. necatrix coverage is ~2% — the noise floor). Used by
# find_interstitial_arrays on top of the annotated gene/intron mask.
rule make_unannotated_masks:
    input:
        manifest="work/results/all_species_raw_summary.tsv",
        tara=["data/raw/tara/.fna.done"],
        refseq="data/raw/refseq/.done",
    output:
        touch("work/results/masks/.done"),
    params:
        outdir="work/results/masks",
        min_orf=450,
    threads: THREADS
    conda:
        ENV
    shell:
        r"""
        mkdir -p {params.outdir}
        python scripts/make_unannotated_mask.py \
            --manifest {input.manifest} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --outdir {params.outdir} \
            --min-orf-nt {params.min_orf}
        """


# Interstitial telomeric arrays: non-terminal (>=5 kb from contig end),
# non-genic, non-intronic. Exclude (annotated genes ∪ introns) ∪ 6-frame
# ORF mask (work/results/masks/{gid}.bed, ORFs >= 450 nt). Emits 5'/3'
# boundary 6- and 12-mers.
rule find_interstitial_arrays:
    input:
        manifest="work/results/all_species_raw_summary.tsv",
        tara=["data/raw/tara/.fna.done"],
        refseq="data/raw/refseq/.done",
        masks="work/results/masks/.done",
    output:
        "work/results/interstitial_arrays.tsv",
    params:
        mask_dir="work/results/masks",
    threads: THREADS
    conda:
        ENV
    shell:
        r"""
        python scripts/find_interstitial_arrays.py \
            --manifest {input.manifest} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --motifs {TELOMERE_MOTIFS} \
            --mask-dir {params.mask_dir} \
            --min-array-len {INTERSTITIAL_MIN_ARRAY_LEN} \
            --out {output} --threads {threads}
        """


# STREME: --dna, w 4-12, 10 motifs, p<0.05, deterministic seed.
STREME_ARGS = "--dna --minw 4 --maxw 12 --nmotifs 10 --thresh 0.05 --seed 1"


# ── TERT by deep homology (miniprot + Pfam domain architecture) ─────────────
# BLAST/tblastn fail on apicomplexan TERT. Field-standard route: intron-aware
# miniprot with apicomplexan TERT seeds, confirmed by RBD (PF12009) UPSTREAM of
# RT (PF00078) -- rules out Eimeria retroelement RTs. Iterative re-seeding from
# confirmed proteins. Toxo gondii = positive control (sister; TERT unannotated).
TERT_GENOME_IDS = config.get("tert_genome_ids") or [
    "GCF_000499385.1",   # Eimeria necatrix
    "GCF_000499425.1",   # Eimeria acervulina
    "GCF_000499545.2",   # Eimeria tenella
    "GCF_000499605.1",   # Eimeria maxima
    "GCF_000499745.2",   # Eimeria mitis
    "GCF_000006565.2",   # Toxoplasma gondii (sister-taxon positive control)
]


rule fetch_tert_seeds_hmms:
    output:
        seeds="work/results/tert_deep_homology/refs/tert_seeds.faa",
        trbd="work/results/tert_deep_homology/refs/PF12009.hmm",
        rt="work/results/tert_deep_homology/refs/PF00078.hmm",
    shell:
        r"""
        python scripts/fetch_tert_seeds_hmms.py \
            --outdir work/results/tert_deep_homology/refs
        """


rule find_tert:
    input:
        seeds="work/results/tert_deep_homology/refs/tert_seeds.faa",
        trbd="work/results/tert_deep_homology/refs/PF12009.hmm",
        rt="work/results/tert_deep_homology/refs/PF00078.hmm",
        refseq="data/raw/refseq/.done",
        tara="data/raw/tara/.fna.done",
    output:
        tsv=TERT_DEEP_HOMOLOGY_TSV,
        faa="work/results/tert_deep_homology/all_confirmed.faa",
    params:
        gids=",".join(TERT_GENOME_IDS),
    threads: THREADS
    conda:
        ENV
    shell:
        r"""
        python scripts/find_tert_deep_homology.py \
            --genome-ids {params.gids} \
            --seeds {input.seeds} --trbd-hmm {input.trbd} --rt-hmm {input.rt} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --outdir work/results/tert_deep_homology \
            --iterate 3 --threads {threads}
        """


# ── Telotron ortholog alignment ──────────────────────────────────────────────
# Per telotron-host gene: intron-aware miniprot maps the ortholog in each panel
# genome, aligns it in protein space, then DNA-aligns the telotron vs the
# orthologous intron. Resolves fill-vs-create. Frameshifts/stops at the homologous
# junction are flagged (not corrected).
_TELO_ORTHO = config.get("telotron_ortholog", {}) or {}
_TELO_ORTHO_FOCAL_IDS = _TELO_ORTHO.get("focal_ids", []) or []
_TELO_ORTHO_PANEL_IDS = _TELO_ORTHO.get("ortholog_ids", []) or []
_TELO_ORTHO_FOCAL = ",".join(_TELO_ORTHO_FOCAL_IDS)
_TELO_ORTHO_PANEL = ",".join(_TELO_ORTHO_PANEL_IDS)
# panel-split for the canonical locus_text: within-Eimeria = focal genomes,
# outgroup = the telotron-negative sisters (ortholog_ids minus the cross-Eimeria).
_TELO_WITHIN_PANEL = ",".join(_TELO_ORTHO_FOCAL_IDS)
_TELO_OUTGROUP_PANEL = ",".join([x for x in _TELO_ORTHO_PANEL_IDS if x not in set(_TELO_ORTHO_FOCAL_IDS)])


rule telotron_ortholog_align:
    input:
        final="work/results/final_telotron_set_architecture.tsv",
        tara=["data/raw/tara/.fna.done"],
        refseq="data/raw/refseq/.done",
    output:
        touch(TELOTRON_ORTHO_SENTINEL),
    params:
        focal=_TELO_ORTHO_FOCAL,
        panel=_TELO_ORTHO_PANEL,
        min_array=_TELO_ORTHO.get("min_array_bp", 50),
        min_id=_TELO_ORTHO.get("min_identity", 0.20),
        tol=_TELO_ORTHO.get("intron_tol", 8),
    threads: THREADS
    conda:
        ENV
    shell:
        r"""
        mkdir -p work/results/telotron_orthologs
        python scripts/telotron_ortholog_align.py \
            --final {input.final} \
            --focal-ids {params.focal} --ortholog-ids {params.panel} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --outdir work/results/telotron_orthologs \
            --min-array-bp {params.min_array} --min-ident {params.min_id} \
            --intron-tol {params.tol} --threads {threads}
        touch {output}
        """


# Compiled per-locus figure: one page per telotron with the flanking-exon protein
# alignment (zoomed to the intron) over the telotron-vs-orthologous-intron DNA alignment.
rule plot_telotron_ortholog_loci:
    input:
        TELOTRON_ORTHO_SENTINEL,
    output:
        TELOTRON_ORTHO_LOCI_PDF,
    conda:
        ENV
    shell:
        r"""
        mkdir -p work/results/figures
        python scripts/plot_telotron_ortholog_loci.py \
            --ortho-dir work/results/telotron_orthologs \
            --out {output} --window 28 --max-rows 12
        """


# Per-locus TEXT dump: unaligned/aligned × DNA/aa (flanks + intron), grouped into
# category folders (intron_present_nontelo / telotron_present / intron_absent / uncertain).
rule telotron_ortholog_textdump:
    input:
        TELOTRON_ORTHO_SENTINEL,
    output:
        touch(TELOTRON_ORTHO_TEXT),
    params:
        within=_TELO_WITHIN_PANEL,
        outgroup=_TELO_OUTGROUP_PANEL,
        v2=TELOTRON_ORTHO_V2,
    threads: THREADS
    conda:
        ENV
    shell:
        r"""
        # one canonical panel-split locus_text (within_eimeria + outgroup), reading the
        # shared align artifacts under work/results/telotron_orthologs.
        python scripts/telotron_ortholog_textdump.py \
            --ortho-dir work/results/telotron_orthologs --threads {threads} \
            --panel-ids {params.within} \
            --out-dir {params.v2}/locus_text/within_eimeria
        python scripts/telotron_ortholog_textdump.py \
            --ortho-dir work/results/telotron_orthologs --threads {threads} \
            --panel-ids {params.outgroup} \
            --out-dir {params.v2}/locus_text/outgroup
        """


# ── Splice-signal sequence logos (±10 bp of GT donor and AG acceptor) ──────
#   {set}='telotron'      + {cat}=''                 -> all telotrons
#   {set}='telotron'      + {cat}='_GG/_GC/_CG/_CC'  -> architecture-filtered
#   {set}='non_telotron'  + {cat}=''                 -> non-telotron control
rule plot_splice_signal_logos:
    input:
        sentinel=lambda w: CTRL_FLANKED_SENTINEL if w.set == "non_telotron" else EXTRACT_FASTA_SENTINEL,
        ann=lambda w: "work/results/final_telotron_set_architecture.tsv" if w.cat else [],
    output:
        touch("work/results/figures/{set}_splice_logos{cat}/.done"),
    wildcard_constraints:
        set=r"non_telotron|telotron",
        cat=r"|_GG|_GC|_CG|_CC",
    params:
        flanked_dir=lambda w: f"work/results/{w.set}_flanked",
        filter_args=lambda w: (
            f"--annotation-tsv work/results/final_telotron_set_architecture.tsv "
            f"--filter-col {TELOTRON_CAT[w.cat[1:]][0]} --filter-value \"{TELOTRON_CAT[w.cat[1:]][1]}\""
            if w.set == "telotron" and w.cat else ""
        ),
    conda: ENV
    shell:
        r"""
        mkdir -p work/results/figures/{wildcards.set}_splice_logos{wildcards.cat}
        python scripts/plot_splice_signal_logos.py \
            --flanked-dir {params.flanked_dir} \
            --outdir work/results/figures/{wildcards.set}_splice_logos{wildcards.cat} \
            {params.filter_args}
        """


# Definitely-non-telotron intron control set: introns with telomeric_frac < 10%
# from the same positive species, sampled per genome. Produced TSV mirrors the
# architecture-table schema so extract_telotron_fasta.py can ingest it.
rule build_non_telotron_controls:
    input:
        introns="work/results/all_introns_scanned.tsv",
        final="work/results/final_telotron_set.tsv",
    output:
        tsv="work/results/non_telotron_controls.tsv",
    conda: ENV
    shell:
        r"""
        python scripts/build_non_telotron_controls.py \
            --introns {input.introns} --final {input.final} \
            --out {output.tsv} --max-frac 0.01 --n-per-species 5000 \
            --all-genomes
        """


# Per-species and per-architecture FASTAs + flanked .txt for the control set,
# under work/results/non_telotron_fasta and work/results/non_telotron_flanked (architecture
# is always "control"). Reuses the telotron extractor.
rule extract_non_telotron_fasta:
    input:
        controls="work/results/non_telotron_controls.tsv",
        tara=["data/raw/tara/.fna.done"],
        refseq="data/raw/refseq/.done",
    output:
        touch(CTRL_FLANKED_SENTINEL),
    conda: ENV
    shell:
        r"""
        mkdir -p work/results/non_telotron_fasta work/results/non_telotron_flanked
        python scripts/extract_telotron_fasta.py \
            --final {input.controls} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --fasta-dir work/results/non_telotron_fasta \
            --flanked-dir work/results/non_telotron_flanked
        """


# Bundle all final TSVs + figures + manifest into a single zip for sharing.
PACKAGE_INPUTS = [
    "work/results/final_telotron_set.tsv",
    "work/results/final_telotron_set_dedup.tsv",
    "work/results/final_telotron_set_architecture.tsv",
    "work/results/final_species_summary.tsv",
    "work/results/final_negative_controls.tsv",
    "work/results/boundary_kmer_enrichment.tsv",
    "work/results/boundary_kmers_by_architecture.tsv",
    "work/results/distance_to_end.tsv",
    "work/results/architecture_summary.tsv",
    "work/results/dedup_log.tsv",
    "work/results/interstitial_arrays.tsv",
    "work/results/linker_blast_hits_own_genome.tsv",
    "work/results/linker_blast_hits_all_genomes.tsv",
    *FIGURES,
    "work/manifests/all_genomes.tsv",
]


rule package:
    input:
        PACKAGE_INPUTS,
    output:
        "work/results/telotron_pipeline_outputs.zip",
    shell:
        "zip -qj {output} {input}"


# ── Self-contained HTML report ──────────────────────────────────────────────
# Aggregates the pipeline outputs into a single HTML file with sortable tables,
# inline base64 figures, colored MSAs, and collapsible sections. Inputs are the
# same as `package` plus a handful of figure/MSA roots that the script samples.
rule pipeline_report:
    input:
        PACKAGE_INPUTS,
        EXTRACT_FASTA_SENTINEL,
    output:
        "work/results/pipeline_report.html",
    conda: ENV
    shell:
        r"""
        python scripts/build_pipeline_report.py \
            --results work/results \
            --out {output}
        """


# ── Architecture / linker / interstitial-ITS analyses ──────────────────────
# Aggregate target: `architecture_analyses`. Downstream of telotron_ortholog +
# interstitial-array arms; chains off the v3 locus_text sentinel.

V2_LOCUS_TEXT = TELOTRON_ORTHO_TEXT   # textdump sentinel; the scripts read the dir themselves


rule interstitial_ortholog_textdump:
    input:
        "work/results/interstitial_arrays.tsv",
        "data/raw/refseq/.done",
    output:
        directory("work/results/interstitial_orthologs/locus_text"),
    threads: THREADS
    conda:
        ENV
    shell:
        r"""
        python scripts/interstitial_ortholog_textdump.py \
            --threads {threads} \
            --min-its-units {INTERSTITIAL_MIN_ITS_UNITS}
        """


rule linker_segmentation:
    input:
        V2_LOCUS_TEXT,
    output:
        "work/results/mechanism_deepdive/linker_segmentation.tsv",
        "work/results/mechanism_deepdive/linker_segmentation.jsonl",
        "work/results/mechanism_deepdive/architecture_per_locus.tsv",
    conda:
        ENV
    shell:
        r"""
        python scripts/linker_segmentation.py
        """


rule cluster_linkers:
    input:
        "work/results/mechanism_deepdive/linker_segmentation.tsv",
    output:
        "work/results/mechanism_deepdive/linker_clusters.tsv",
        "work/results/mechanism_deepdive/linker_recurrence_summary.txt",
    conda:
        ENV
    shell:
        r"""
        python scripts/cluster_linkers.py
        """


rule mask_telotron_arrays:
    input:
        V2_LOCUS_TEXT,
    output:
        "work/results/telotron_masked/telotron_masked.fasta",
        "work/results/telotron_masked/telotron_masked.msa.fasta",
        "work/results/telotron_masked/telotron_masked_msa.txt",
        "work/results/telotron_masked/telotron_masked_msa.html",
    conda:
        ENV
    shell:
        r"""
        python scripts/mask_telotron_arrays.py
        """


rule architecture_analyses:
    input:
        rules.interstitial_ortholog_textdump.output,
        rules.cluster_linkers.output,
        rules.mask_telotron_arrays.output,


# ── Assembly-based telomere boundary detection ──────────────────────────────
# Field-standard motif sliding-window scan on assembled long-read contigs.
# Reports the most-internal motif position contiguous with each contig end.
# Config block `telomere_boundaries` sets per-species assembly + motif.
_TELOBND = config.get("telomere_boundaries", {}) or {}
_TELOBND_ASM = _TELOBND.get("assemblies", {}) or {}


# ── Long-read data + Telogator2 allele-specific telomere-length analysis ────
# data/raw/longread/manifest.tsv is the source of truth. Telogator2 consumes
# the nanopore subset, anchored to the custom Eimeria subtelomere reference at
# work/results/telogator2_ref/. Replaces the retired TARPON pipeline.
# Targets: `longread_data` (download all), `telogator2` (run on every ONT run).

import csv as _csv

LR_MANIFEST = "data/raw/longread/manifest.tsv"

def _read_lr_manifest():
    import os
    if not os.path.exists(LR_MANIFEST):
        return []
    with open(LR_MANIFEST) as fh:
        return list(_csv.DictReader(fh, delimiter="\t"))

LR_RUNS = _read_lr_manifest()
TG2_RUNS = [r for r in LR_RUNS if r["platform"] == "nanopore"]

# accession -> list of local fastq paths (handles multi-file runs e.g. SRR24971026)
_TG2_FILES = {}
for _r in TG2_RUNS:
    _TG2_FILES.setdefault(_r["accession"], []).append("data/raw/longread/" + _r["fastq_local"])

# accession -> species (one row per accession is fine, multi-file runs share species)
_TG2_SPECIES = {r["accession"]: r["species"] for r in TG2_RUNS}

# URL lookup keyed by local destination path
_LR_URL_BY_DEST = {("data/raw/longread/" + r["fastq_local"]): ("https://" + r["fastq_url"]) for r in LR_RUNS}


# ===================================================================================
# Nucleosome flanking-site inputs. Builds insertion-site + flank FASTAs used by
# nucleosome_features (composition/periodicity panel) and nucleosome_withingene
# (within-gene sibling-intron control). The NuPoP occupancy arm was retired
# 2026-07-14 — memory (telotron_nucleosome_nupop_2026-06-08) recorded the
# linker interpretation as artifact; only the composition/periodicity signals
# in nucleosome_features held up.
rule nucleosome_inputs:
    input:
        arch="work/results/final_telotron_set_architecture.tsv",
        controls="work/results/non_telotron_controls.tsv",
        tara=["data/raw/tara/.fna.done"],
        refseq="data/raw/refseq/.done",
    output:
        manifest="work/results/nucleosome/manifest.tsv",
        control_manifest="work/results/nucleosome/control_manifest.tsv",
    conda:
        ENV
    shell:
        r"""
        python scripts/nucleosome_inputs.py \
            --table {input.arch} --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --out work/results/nucleosome
        python scripts/nucleosome_control_inputs.py \
            --controls {input.controls} --telo-manifest {output.manifest} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara --out work/results/nucleosome
        """


# ── Downstream analysis arm (was hand-run; now wired) ───────────────────────
# These rules cover the previously DAG-disconnected analysis scripts whose inputs
# ARE regenerable from the core pipeline (the architecture table, the controls
# table, the RefSeq GFFs, the nucleosome manifests). They write to the script's
# own output paths (mostly analysis/ + work/results/...). Scripts that depend on
# non-regenerable external data (ONT/Hi-C/RNA-seq reads, the work/old/ archive,
# hand-built telogator2 refs) are NOT wired — they are documented under the
# `analysis_external_inputs` note below, because a rule with an input no other
# rule can produce would only ever error.

ARCH_TSV = "work/results/final_telotron_set_architecture.tsv"

# RNA-seq gene-coverage arm (config-driven; SRA accessions in config["rnaseq"]).
_RNASEQ = config.get("rnaseq", {}) or {}
_RNASEQ_SP = _RNASEQ.get("species", {}) or {}
RNASEQ_THREADS = int(_RNASEQ.get("threads", THREADS))
RNASEQ_COV = {sp: f"work/results/rnaseq/{sp}_gene_cov.tsv" for sp in _RNASEQ_SP}


# Per-species RNA-seq -> gene depth-sum TSV (samtools bedcov). Replaces the hand-run
# work/necatrix_rnaseq/run_pipeline.sh and the ephemeral /tmp/eten_gene_cov.tsv.
rule rnaseq_gene_coverage:
    input:
        genome=lambda w: _RNASEQ_SP[w.species]["genome"],
        gff=lambda w: _RNASEQ_SP[w.species]["gff"],
        refseq="data/raw/refseq/.done",
    output:
        "work/results/rnaseq/{species}_gene_cov.tsv",
    params:
        srr=lambda w: ",".join(_RNASEQ_SP[w.species]["srr"]),
        maxsize=lambda w: _RNASEQ_SP[w.species].get("max_size", "30g"),
    threads: RNASEQ_THREADS
    conda:
        ENV
    shell:
        r"""
        python scripts/rnaseq_gene_coverage.py \
            --species {wildcards.species} \
            --genome {input.genome} --gff {input.gff} \
            --srr {params.srr} --max-size {params.maxsize} \
            --threads {threads} --out {output}
        """


# Insertion-site composition/periodicity feature panel (telomere-masked flanks,
# BH-FDR). Reads the nucleosome manifests + per-locus FASTAs from nucleosome_inputs.
rule nucleosome_features:
    input:
        manifest="work/results/nucleosome/manifest.tsv",
        control_manifest="work/results/nucleosome/control_manifest.tsv",
    output:
        "work/results/nucleosome/nucleosome_feature_summary.png",
    conda:
        ENV
    shell:
        r"""
        python scripts/nucleosome_features.py
        """


# Within-gene control: telotron introns vs same-gene sibling introns vs random
# non-host introns (separates gene-class confound from local targeting). Prints
# to stdout; we tee it to a log file so the rule has a concrete output.
rule nucleosome_withingene:
    input:
        arch=ARCH_TSV,
        controls="work/results/non_telotron_controls.tsv",
        refseq="data/raw/refseq/.done",
        tara=["data/raw/tara/.fna.done"],
    output:
        "work/results/nucleosome/withingene_control.txt",
    conda:
        ENV
    shell:
        r"""
        python scripts/nucleosome_withingene.py | tee {output}
        """




# Length distribution by (lineage x architecture) — BH-FDR corrected, single-MAG
# caveat in the output. Per-arm burst-length subfigure is emitted from the same
# TSVs by the same script.
rule length_distribution_by_arch:
    input:
        ARCH_TSV,
    output:
        summary="work/results/length_distribution/summary_by_lineage_arch.tsv",
        linker="work/results/length_distribution/linker_vs_nolinker.tsv",
        per_arch="work/results/length_distribution/psw_vs_eimeria_per_arch.tsv",
        hist="work/results/length_distribution/length_hist_by_lineage_arch.pdf",
        per_arm="work/results/length_distribution/per_arm_burst_length.pdf",
    conda:
        ENV
    shell:
        r"""
        TELOTRON_ROOT="$(pwd)" python scripts/length_distribution_by_arch_lineage.py
        """


# Expression-vs-telotron-presence figure. Gene-level coverage comes from the
# reproducible rnaseq_gene_coverage rule; the optional locus-level panel reads
# data/raw/rnaseq_splice_2026/per_locus_counts.tsv (manual build) and is skipped
# gracefully when absent. One consolidated script (was 2: presence + final).
rule telotron_expr_figures:
    input:
        eten="work/results/rnaseq/eten_gene_cov.tsv",
        necatrix="work/results/rnaseq/necatrix_gene_cov.tsv",
        arch=ARCH_TSV,
        refseq="data/raw/refseq/.done",
    output:
        "work/results/figures/telotron_expression.png",
    conda:
        ENV
    shell:
        r"""
        mkdir -p work/results/figures
        python scripts/telotron_expression.py
        """


# Host-gene-class dissection (host vs disjoint non-host intron-count + gene length;
# per-intron rate vs gene intron count + 5'->3' position). One consolidated
# analysis for the intron-rich long-gene bias (was 2 rules: telotron_gene_bias
# + telotron_per_intron, merged 2026-07-14).
rule telotron_gene_class:
    input:
        arch=ARCH_TSV,
        refseq="data/raw/refseq/.done",
    output:
        "work/results/figures/telotron_gene_class.png",
    conda:
        ENV
    shell:
        r"""
        mkdir -p work/results/figures
        python scripts/telotron_gene_class.py
        """


# Aggregate target for the wired analysis arm.
rule analysis_arm:
    input:
        rules.nucleosome_features.output,
        rules.nucleosome_withingene.output,
        rules.telotron_gene_class.output,
        rules.length_distribution_by_arch.output,
        rules.telotron_expr_figures.output,
#
# detect_telomere_boundaries.py is fully CLI-parameterised and driven by the
# telomere_boundaries rule.
