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
NUPOP_ENV = "envs/nupop.yaml"

# Sentinels (path = ".../<name>/.done" or the canonical output file).
TERMINAL_DENSITY_SENTINEL       = "work/results/figures/terminal_motif/.done"
ARCH_KMER_SENTINEL              = "work/results/figures/boundary_kmers_by_arch/.done"
INTERSTITIAL_KMER_SENTINEL      = "work/results/figures/interstitial_boundary_kmers/.done"
INTERSTITIAL_LOGO_SENTINEL      = "work/results/figures/interstitial_boundary_logos/.done"
PIPELINE_STAGES_SENTINEL        = "work/results/figures/pipeline_stages/.done"
EXTRACT_FASTA_SENTINEL          = "work/results/telotron_fasta/.done"
MSA_SENTINEL                    = "work/results/msa_regions/.done"
CTRL_FLANKED_SENTINEL           = "work/results/non_telotron_fasta/.done"
TELOTRON_BOUNDARY_KMER_SENTINEL = "work/results/figures/telotron_boundary_kmers/.done"
CTRL_BOUNDARY_KMER_SENTINEL     = "work/results/figures/non_telotron_boundary_kmers/.done"
TELOTRON_SPLICE_LOGO_SENTINEL   = "work/results/figures/telotron_splice_logos/.done"
CTRL_SPLICE_LOGO_SENTINEL       = "work/results/figures/non_telotron_splice_logos/.done"
COMPOSITE_KMER_SENTINEL         = "work/results/figures/composite_boundary_kmers/.done"
COMPOSITE_LOGO_SENTINEL         = "work/results/figures/composite_boundary_logos/.done"
COMPOSITE_5P3P_KMER_SENTINEL    = "work/results/figures/composite_boundary_kmers_5p3p/.done"
COMPOSITE_5P3P_LOGO_SENTINEL    = "work/results/figures/composite_boundary_logos_5p3p/.done"
ARRAY_LEN_DIST_PNG              = "work/results/figures/array_length_distribution.png"
ARRAY_LEN_DIST_MIN40_PNG        = "work/results/figures/array_length_distribution_min40.png"
STREME_TELO_SENTINEL            = "work/results/streme/telotrons/streme.xml"
STREME_LINKER_SENTINEL          = "work/results/streme/linkers/streme.xml"
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


# Default target = core survey through `package` + `pipeline_report`. Exploratory
# branches are on-demand targets: `sequences` (FASTAs + MSAs), `figures_extra`
# (boundary-kmer/logo/composite/density), `motif_discovery` (STREME on telotrons
# + linkers), `telomerase_search` (TERT deep-homology), `orthologs` (host-gene
# ortholog align + per-locus PDF + textdump), `architecture_analyses` (mechanism
# deep-dive), `telogator2` (long-read telomere lengths), or `everything`.
rule all:
    input:
        "work/results/telotron_pipeline_outputs.zip",
        "work/results/pipeline_report.html",


rule sequences:
    input:
        EXTRACT_FASTA_SENTINEL,
        MSA_SENTINEL,


rule figures_extra:
    input:
        TERMINAL_DENSITY_SENTINEL,
        ARCH_KMER_SENTINEL,
        INTERSTITIAL_KMER_SENTINEL,
        INTERSTITIAL_LOGO_SENTINEL,
        ARRAY_LEN_DIST_PNG,
        ARRAY_LEN_DIST_MIN40_PNG,
        PIPELINE_STAGES_SENTINEL,
        TELOTRON_SPLICE_LOGO_SENTINEL,
        CTRL_SPLICE_LOGO_SENTINEL,
        CTRL_BOUNDARY_KMER_SENTINEL,
        TELOTRON_BOUNDARY_KMER_SENTINEL,
        COMPOSITE_KMER_SENTINEL,
        COMPOSITE_LOGO_SENTINEL,
        COMPOSITE_5P3P_KMER_SENTINEL,
        COMPOSITE_5P3P_LOGO_SENTINEL,


rule motif_discovery:
    input:
        STREME_TELO_SENTINEL,
        STREME_LINKER_SENTINEL,


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


# Core survey + every exploratory branch (the pre-refactor default target).
rule everything:
    input:
        rules.all.input,
        rules.sequences.input,
        rules.figures_extra.input,
        rules.motif_discovery.input,
        rules.telomerase_search.input,
        rules.orthologs.input,


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


# Per-species "where do telomere repeats sit?" plots — one figure per species,
# one panel per top-N contig, x = genomic coordinate, y = rolling hit density of the
# species's identified terminal motif. Telotron candidate loci overlaid as ticks.
rule terminal_motif_figures:
    input:
        summary="work/results/all_species_raw_summary.tsv",
        loci="work/results/all_telotron_loci.tsv",
        tara=["data/raw/tara/.fna.done"],
        refseq="data/raw/refseq/.done",
    output:
        touch(TERMINAL_DENSITY_SENTINEL),
    threads: THREADS
    conda:
        ENV
    shell:
        r"""
        mkdir -p work/results/figures/terminal_motif
        python scripts/plot_terminal_motif_density.py \
            --summary {input.summary} --loci {input.loci} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --motifs {TELOMERE_MOTIFS} \
            --outdir work/results/figures/terminal_motif \
            --threads {threads}
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


# BLAST linker sequences from linker-architecture telotrons (1) against their own
# species genome and (2) against a concatenated all-species DB. Outputs hit TSVs
# annotated with the query genome_id so cross-species recurrences are findable.
rule blast_linkers:
    input:
        arch="work/results/final_telotron_set_architecture.tsv",
        tara=["data/raw/tara/.fna.done"],
        refseq="data/raw/refseq/.done",
    output:
        own="work/results/linker_blast_hits_own_genome.tsv",
        all="work/results/linker_blast_hits_all_genomes.tsv",
    threads: THREADS
    conda:
        ENV
    shell:
        r"""
        mkdir -p work/results/blast_linkers
        python scripts/blast_linkers.py \
            --arch-tsv {input.arch} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --workdir work/results/blast_linkers \
            --out-own {output.own} --out-all {output.all} \
            --threads {threads}
        """


# Per-species × per-architecture MAFFT MSAs of the locus regions
# (upstream_50 | intron|arm1[+linker+arm2]|arm1[+arm2] | downstream_50).
# Both raw and homopolymer-compressed alignments; combined.aln.txt collapses
# everything into one fixed-width view with spaces at region boundaries.
rule msa_telotron_regions:
    input:
        arch="work/results/final_telotron_set_architecture.tsv",
        tara=["data/raw/tara/.fna.done"],
        refseq="data/raw/refseq/.done",
    output:
        touch(MSA_SENTINEL),
    threads: 4
    conda:
        ENV
    shell:
        r"""
        mkdir -p work/results/msa_regions
        python scripts/msa_telotron_regions.py \
            --arch-tsv {input.arch} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --min-array {ARCH_MIN_ARRAY} \
            --outdir work/results/msa_regions --threads {threads}
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


# Per-species figure: top boundary 6-mers (GT-XXXX donor, XXXX-AG acceptor,
# and 4+2 linker-edge 6-mers) split by architecture.
rule plot_boundary_kmers_by_arch:
    input:
        kmers="work/results/boundary_kmers_by_architecture.tsv",
        species="work/results/final_species_summary.tsv",
    output:
        touch(ARCH_KMER_SENTINEL),
    conda:
        ENV
    shell:
        r"""
        mkdir -p work/results/figures/boundary_kmers_by_arch
        python scripts/plot_boundary_kmers_by_arch.py \
            --kmers {input.kmers} --species {input.species} \
            --outdir work/results/figures/boundary_kmers_by_arch
        touch {output}
        """


# Per-species 5'/3' boundary 6-mers and 12-mers of interstitial (non-genic,
# non-intronic, non-terminal) telomeric repeat arrays.
#   {cat}=''                       -> all arrays
#   {cat}='_GG'/'_GC'/'_CG'/'_CC'  -> filter by 5'/3' orientation category
rule plot_interstitial_boundary_kmers:
    input:
        arrays="work/results/interstitial_arrays.tsv",
    output:
        touch("work/results/figures/interstitial_boundary_kmers{cat}/.done"),
    wildcard_constraints:
        cat=r"|_GG|_GC|_CG|_CC",
    params:
        filter_args=lambda w: (
            f"--filter-col {INTERSTITIAL_CAT[w.cat[1:]][0]} --filter-value \"{INTERSTITIAL_CAT[w.cat[1:]][1]}\""
            if w.cat in ("_GG", "_GC", "_CG", "_CC") else ""
        ),
    conda: ENV
    shell:
        r"""
        mkdir -p work/results/figures/interstitial_boundary_kmers{wildcards.cat}
        python scripts/plot_interstitial_boundary_kmers.py \
            --arrays {input.arrays} \
            --outdir work/results/figures/interstitial_boundary_kmers{wildcards.cat} \
            {params.filter_args}
        """


# Per-species sequence logos (information bits, logomaker) of the 5'-flank
# + first repeat unit and last repeat unit + 3'-flank for interstitial arrays.
#   {cat}=''                       -> all arrays
#   {cat}='_GG'/'_GC'/'_CG'/'_CC'  -> filter by 5'/3' orientation category
rule plot_interstitial_boundary_logos:
    input:
        arrays="work/results/interstitial_arrays.tsv",
    output:
        touch("work/results/figures/interstitial_boundary_logos{cat}/.done"),
    wildcard_constraints:
        cat=r"|_GG|_GC|_CG|_CC",
    params:
        filter_args=lambda w: (
            f"--filter-col {INTERSTITIAL_CAT[w.cat[1:]][0]} --filter-value \"{INTERSTITIAL_CAT[w.cat[1:]][1]}\""
            if w.cat else ""
        ),
    conda: ENV
    shell:
        r"""
        mkdir -p work/results/figures/interstitial_boundary_logos{wildcards.cat}
        python scripts/plot_interstitial_boundary_logos.py \
            --arrays {input.arrays} \
            --outdir work/results/figures/interstitial_boundary_logos{wildcards.cat} \
            --flank-len 10 \
            {params.filter_args}
        """


rule plot_array_length_distribution:
    """Per-species histograms: ITS array_len | telotron telomeric_bases |
    non-telotron intron_len. {suffix}='' = unfiltered, '_min40' = drop <40 bp."""
    input:
        interstitial="work/results/interstitial_arrays.tsv",
        telotrons="work/results/final_telotron_set.tsv",
        non_telotrons="work/results/non_telotron_controls.tsv",
    output:
        "work/results/figures/array_length_distribution{suffix}.png",
    wildcard_constraints:
        suffix=r"|_min40",
    params:
        min_len=lambda w: "--min-len 40" if w.suffix == "_min40" else "",
    conda:
        ENV
    shell:
        r"""
        python scripts/plot_array_length_distribution.py \
            --interstitial {input.interstitial} \
            --telotrons {input.telotrons} \
            --non-telotrons {input.non_telotrons} \
            {params.min_len} \
            --out {output}
        """


# Composite figure (PDF + tall PNG): per-species rows, intron architecture
# boundary panel on the left, interstitial array boundary panel on the right.
rule plot_combined_boundary_kmers:
    input:
        species="work/results/final_species_summary.tsv",
        intron_dir=ARCH_KMER_SENTINEL,
        interst_dir=INTERSTITIAL_KMER_SENTINEL,
    output:
        pdf="work/results/figures/boundary_kmers_combined.pdf",
        png="work/results/figures/boundary_kmers_combined.png",
    conda:
        ENV
    shell:
        r"""
        python scripts/plot_combined_boundary_kmers.py \
            --species {input.species} \
            --intron-dir work/results/figures/boundary_kmers_by_arch \
            --interstitial-dir work/results/figures/interstitial_boundary_kmers \
            --out-pdf {output.pdf} --out-png {output.png}
        """


# Pipeline-stage summary figures (manifests, scan, filter funnel, distance,
# architecture).
rule plot_pipeline_stages:
    input:
        manifest="work/manifests/all_genomes.tsv",
        raw="work/results/all_species_raw_summary.tsv",
        sp="work/results/final_species_summary.tsv",
        dist="work/results/distance_to_end.tsv",
        arch="work/results/final_telotron_set_architecture.tsv",
    output:
        touch(PIPELINE_STAGES_SENTINEL),
    conda:
        ENV
    shell:
        r"""
        mkdir -p work/results/figures/pipeline_stages
        python scripts/plot_pipeline_stages.py \
            --manifests {input.manifest} \
            --raw-summary {input.raw} \
            --species-final {input.sp} \
            --distance {input.dist} \
            --architecture-loci {input.arch} \
            --outdir work/results/figures/pipeline_stages
        touch {output}
        """


# ── STREME (MEME suite) de novo motif discovery ──────────────────────────────
# Three positive sets: telotron introns, non-telomeric introns, and linkers.
# STREME generates its own shuffled control.

rule build_streme_inputs:
    input:
        telotrons_done=EXTRACT_FASTA_SENTINEL,
        introns="work/results/all_introns_scanned.tsv",
        linkers="work/results/blast_linkers/linker_queries/_all_linkers.fa",
    output:
        telo="work/results/streme_inputs/telotrons.fa",
        linkers="work/results/streme_inputs/linkers.fa",
    conda: ENV
    shell:
        r"""
        python scripts/build_streme_inputs.py \
            --telotron-fasta-dir work/results/telotron_fasta \
            --introns-tsv {input.introns} \
            --linkers-fa {input.linkers} \
            --refseq-dir data/raw/refseq --tara-dir data/raw/tara \
            --outdir work/results/streme_inputs \
            --n-non-telo 0
        """


# STREME: --dna, w 4-12, 10 motifs, p<0.05, deterministic seed.
STREME_ARGS = "--dna --minw 4 --maxw 12 --nmotifs 10 --thresh 0.05 --seed 1"

# {set} = telotrons | linkers
rule streme:
    input: "work/results/streme_inputs/{set}.fa"
    output: "work/results/streme/{set}/streme.xml"
    wildcard_constraints: set=r"telotrons|linkers"
    conda: MEME_ENV
    shell:
        r"""
        streme {STREME_ARGS} --p {input} --oc work/results/streme/{wildcards.set}
        """


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


# Per-species 2x2 boundary k-mer rank plots from the first40/last40 columns.
#   {set}='telotron'      + {cat}=''                 -> all telotrons
#   {set}='telotron'      + {cat}='_GG/_GC/_CG/_CC'  -> architecture-filtered
#   {set}='non_telotron'  + {cat}=''                 -> non-telotron control
rule plot_intron_boundary_kmers:
    input:
        introns=lambda w: (
            "work/results/non_telotron_controls.tsv" if w.set == "non_telotron"
            else "work/results/final_telotron_set_architecture.tsv"
        ),
    output:
        touch("work/results/figures/{set}_boundary_kmers{cat}/.done"),
    wildcard_constraints:
        set=r"non_telotron|telotron",
        cat=r"|_GG|_GC|_CG|_CC",
    params:
        label=lambda w: (
            "non-telotron control (telomeric_frac<0.10)" if w.set == "non_telotron"
            else (f"telotrons 5'/3' = {w.cat[1:]}" if w.cat else "telotrons")
        ),
        filter_args=lambda w: (
            f"--filter-col {TELOTRON_CAT[w.cat[1:]][0]} --filter-value \"{TELOTRON_CAT[w.cat[1:]][1]}\""
            if w.set == "telotron" and w.cat else ""
        ),
    conda: ENV
    shell:
        r"""
        mkdir -p work/results/figures/{wildcards.set}_boundary_kmers{wildcards.cat}
        python scripts/plot_intron_boundary_kmers.py \
            --introns {input.introns} \
            --outdir work/results/figures/{wildcards.set}_boundary_kmers{wildcards.cat} \
            --label "{params.label}" \
            {params.filter_args}
        """


# Composite per-species figures: side-by-side panels.
#   {mode}='kmers' or 'logos'
#   {split}=''       -> 3 panels (interstitial | non_telotron control | telotron)
#   {split}='_5p3p'  -> 9 panels stratified by 5'/3' end-orientation
# Panel paths are asymmetric (interstitial uses _logos for mode=logos but
# telotron/non_telotron use _splice_logos); _composite_panel_args() handles it.
def _composite_panel_args(mode, split):
    interst_dir = f"work/results/figures/interstitial_boundary_{'kmers' if mode == 'kmers' else 'logos'}"
    non_telo_dir = ("work/results/figures/non_telotron_boundary_kmers"
                    if mode == "kmers" else "work/results/figures/non_telotron_splice_logos")
    telo_dir = ("work/results/figures/telotron_boundary_kmers"
                if mode == "kmers" else "work/results/figures/telotron_splice_logos")
    suffix_kmers = " arrays" if mode == "kmers" else " array boundary"
    if split == "":
        panels = [
            (interst_dir, f"interstitial{suffix_kmers}"),
            (non_telo_dir, "non-telotron control introns" if mode == "kmers" else "non-telotron intron splice signal"),
            (telo_dir, "telotrons" if mode == "kmers" else "telotron splice signal"),
        ]
    else:  # _5p3p
        cats = [("GG", "5'=G 3'=G"), ("GC", "5'=G 3'=C (hybrid)"),
                ("CG", "5'=C 3'=G (hybrid)"), ("CC", "5'=C 3'=C")]
        panels = [(f"{interst_dir}_{c}", f"interstitial {l}") for c, l in cats]
        panels.append((non_telo_dir, "non-telotron control" if mode == "kmers" else "non-telotron control splice signal"))
        panels += [(f"{telo_dir}_{c}", f"telotron {l}") for c, l in cats]
    return " ".join(f'--panel "{p}:{lab}"' for p, lab in panels)

def _composite_inputs(mode, split):
    args = _composite_panel_args(mode, split).split('"')[1::2]
    return [a.split(":", 1)[0] + "/.done" for a in args]


rule plot_composite_per_species:
    input:
        sentinels=lambda w: _composite_inputs(w.mode, w.split),
    output:
        touch("work/results/figures/composite_boundary_{mode}{split}/.done"),
    wildcard_constraints:
        mode=r"kmers|logos",
        split=r"|_5p3p",
    params:
        panel_args=lambda w: _composite_panel_args(w.mode, w.split),
    conda: ENV
    shell:
        r"""
        mkdir -p work/results/figures/composite_boundary_{wildcards.mode}{wildcards.split}
        python scripts/plot_composite_per_species.py \
            --outdir work/results/figures/composite_boundary_{wildcards.mode}{wildcards.split} \
            {params.panel_args}
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
    "work/results/figures/boundary_kmers_combined.pdf",
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
        MSA_SENTINEL,
        PIPELINE_STAGES_SENTINEL,
        TELOTRON_BOUNDARY_KMER_SENTINEL,
        CTRL_BOUNDARY_KMER_SENTINEL,
        ARCH_KMER_SENTINEL,
        COMPOSITE_KMER_SENTINEL,
        COMPOSITE_LOGO_SENTINEL,
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


rule ortholog_review_html:
    input:
        V2_LOCUS_TEXT,
    output:
        "work/results/telotron_orthologs_v2/review.html",
    conda:
        ENV
    shell:
        r"""
        python scripts/build_ortholog_review.py
        """


rule architecture_analyses:
    input:
        rules.interstitial_ortholog_textdump.output,
        rules.cluster_linkers.output,
        rules.mask_telotron_arrays.output,
        rules.ortholog_review_html.output,


# ── Assembly-based telomere boundary detection ──────────────────────────────
# Field-standard motif sliding-window scan on assembled long-read contigs.
# Reports the most-internal motif position contiguous with each contig end.
# Config block `telomere_boundaries` sets per-species assembly + motif.
_TELOBND = config.get("telomere_boundaries", {}) or {}
_TELOBND_ASM = _TELOBND.get("assemblies", {}) or {}


rule telomere_boundaries:
    """Per-arm telomere boundary table from assembled contigs."""
    input:
        [v["path"] for v in _TELOBND_ASM.values()],
    output:
        tsv="work/results/telomere_boundaries/per_arm.tsv",
        summary="work/results/telomere_boundaries/summary.md",
        bed="work/results/telomere_boundaries/telomeres.bed",
    params:
        assemblies=",".join(f"{k}:{v['path']}" for k, v in _TELOBND_ASM.items()),
        motifs=",".join(f"{k}:{v['motif']}" for k, v in _TELOBND_ASM.items()),
        scan_kb=_TELOBND.get("scan_kb", 20),
        window_bp=_TELOBND.get("window_bp", 200),
        density_frac=_TELOBND.get("density_frac", 0.70),
        max_gap_bp=_TELOBND.get("max_gap_bp", 100),
        min_units=_TELOBND.get("min_units", 5),
        max_mm=_TELOBND.get("max_mm_per_unit", 1),
    conda: ENV
    shell:
        r"""
        mkdir -p work/results/telomere_boundaries
        python scripts/detect_telomere_boundaries.py \
            --assemblies "{params.assemblies}" \
            --motifs "{params.motifs}" \
            --scan-kb {params.scan_kb} \
            --window-bp {params.window_bp} \
            --density-frac {params.density_frac} \
            --max-gap-bp {params.max_gap_bp} \
            --min-units {params.min_units} \
            --max-mm {params.max_mm} \
            --out-tsv {output.tsv} \
            --out-summary {output.summary} \
            --out-bed {output.bed}
        """


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


rule longread_data:
    """Download every long-read fastq listed in data/raw/longread/manifest.tsv."""
    input:
        ["data/raw/longread/" + r["fastq_local"] for r in LR_RUNS],


rule download_longread_fastq:
    """Generic per-file downloader. Resumable via wget -c."""
    output:
        "data/raw/longread/{platform}/{filename}.fastq.gz",
    params:
        url=lambda wc, output: _LR_URL_BY_DEST[str(output[0])],
    threads: 1
    resources:
        network=1,
    shell:
        r"""
        mkdir -p $(dirname {output})
        wget -q --show-progress=off -t 5 -c -O {output} '{params.url}'
        """


rule telogator2:
    """Aggregate target: run Telogator2 on every nanopore run in the manifest."""
    input:
        [f"work/results/telogator2/{_TG2_SPECIES[acc]}_{acc}/tlens_by_allele.tsv"
         for acc in _TG2_FILES],


rule telogator2_one:
    """Per-accession Telogator2 invocation.

    Telogator2 (Stephens & Kocher 2024, BMC Bioinf 25:194) is an allele-
    specific telomere-length / TVR caller for long reads. It internally
    aligns raw fastq to the supplied subtelomere reference (`-t`) with
    minimap2 and clusters reads by arm.

    Inputs:
      - nanopore fastq(s) (multi-file runs are concatenated upstream).
      - Custom Eimeria subtelomere ref (T2 deliverable, renamed by
        scripts/build_telogator2_ref.py to <species>_chr<N><p|q> form so
        that source/tg_util.py:LEXICO_2_IND parsing succeeds).
      - Maize kmer palette (canonical CCCTAAA, identical to Eimeria
        TTTAGGG revcomp). For E. maxima the canonical is TTAGGG and the
        default human kmers TSV is used instead — see params.kmers.

    Notes (cf. work/results/telogator2/telogator2_cli_spec.md):
      - PYTHONNOUSERSITE=1 is mandatory: a stale user-site Biopython 1.76
        masks the env's 1.87 and crashes Telogator2's pickling.
      - --minimap2 must be an absolute path (no default in Telogator2).
      - WGS nanopore depth is low; -n 4 from the README is appropriate
        for 30x-equivalent runs. Drop to -n 3 for lower-depth datasets.
      - There is no `-bed` flag for ITS coordinates; interstitial
        flagging is hardwired to bundled human ITS positions. We post-
        process tlens_by_allele.tsv against work/results/telogator2_ref/
        its_catalog.bed downstream (not in this rule).
    """
    input:
        fastqs=lambda wc: _TG2_FILES[wc.accession],
        # PER-SPECIES reference, not pooled. Pooling produces cross-species
        # mis-mapping (e.g. E.acervulina ONT reads anchor to E.maxima/tenella
        # subtelomeres ~95% of the time due to subtelomere homology + ONT
        # error). Per-species ref drops fail_reads_unmapped from 100% to 0%.
        ref=lambda wc: f"work/results/telogator2_ref/{wc.species[0]}{wc.species.split('_')[1]}.telogator2.fasta",
    output:
        tsv="work/results/telogator2/{species}_{accession}/tlens_by_allele.tsv",
    params:
        outdir=lambda wc: f"work/results/telogator2/{wc.species}_{wc.accession}",
        merged=lambda wc: f"work/results/telogator2/{wc.species}_{wc.accession}/_merged.fastq.gz",
        kmers=lambda wc: (
            "/scratch1/alex/telogator2/source/resources/kmers.tsv"
            if wc.species == "Eimeria_maxima"     # TTAGGG canonical
            else "/scratch1/alex/telogator2/source/resources/non-human/kmers_maize.tsv"
        ),
        env=config["telogator2"]["env"],
        repo=config["telogator2"]["repo"],
        read_type=config["telogator2"].get("read_type", "ont"),
        min_reads_per_cluster=config["telogator2"].get("min_reads_per_cluster", 2),
        # Loosened filters for shotgun WGS coverage (defaults assume capture).
        # filt-tel default 400 bp drops most shotgun telomere-spanning reads.
        filt_tel=config["telogator2"].get("filt_tel", 100),
        filt_nontel=config["telogator2"].get("filt_nontel", 500),
        filt_sub=config["telogator2"].get("filt_sub", 200),
        min_read_len=config["telogator2"].get("min_read_len", 1000),
        min_kmer_hits=config["telogator2"].get("min_kmer_hits", 4),
    threads: 16
    shell:
        r"""
        mkdir -p {params.outdir}
        if [ $(echo {input.fastqs} | wc -w) -eq 1 ]; then
            ln -sf $(realpath {input.fastqs}) {params.merged}
        else
            cat {input.fastqs} > {params.merged}
        fi

        # --debug-nosubtel skips the subtel-refinement clustering stage,
        # which crashes on non-vertebrate refs (dendrogram label/Z size
        # mismatch). Final anchoring still runs; alleles are anchored from
        # the iter-0 / iter-1 cluster representatives directly.
        PYTHONNOUSERSITE=1 micromamba run -n {params.env} \
            python {params.repo}/telogator2.py \
                -i {params.merged} \
                -o {params.outdir} \
                -t {input.ref} \
                -k {params.kmers} \
                -r {params.read_type} \
                -n {params.min_reads_per_cluster} \
                -p {threads} \
                -l {params.min_read_len} \
                -c {params.min_kmer_hits} \
                --filt-tel {params.filt_tel} \
                --filt-nontel {params.filt_nontel} \
                --filt-sub {params.filt_sub} \
                --debug-nosubtel \
                --minimap2 /home/alex/micromamba/envs/{params.env}/bin/minimap2
        """


# ===================================================================================
# Nucleosome occupancy (NuPoP), replicating Gozashti et al. 2022 (doi:10.1073/pnas.2209766119).
# Adaptive-flank NuPoP (predNuPoP species=0 model=4) over telotron insertion sites and matched
# confident non-telotron introns, run twice: with the element and with the element computationally
# removed. The composition-independent (element-removed) test asks whether telotron insertion sites
# sit in nucleosome-linker DNA (introner-like). The within-array occupancy is confounded by the
# telomeric-repeat composition and is reported but not interpreted. The Tara MAG is included via the
# adaptive flank (its short contigs fail a fixed 5 kb flank). Figure is faceted by species x architecture.
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


# Run predNuPoP on every single-sequence FASTA. NuPoP needs one sequence per file and writes the
# prediction into the input file's directory (run_nupop.R setwd's there to avoid collisions); we
# fan {threads} R workers over the file list.
rule nucleosome_nupop:
    input:
        manifest="work/results/nucleosome/manifest.tsv",
        control_manifest="work/results/nucleosome/control_manifest.tsv",
    output:
        touch("work/results/nucleosome/.nupop.done"),
    threads: THREADS
    conda:
        NUPOP_ENV
    shell:
        r"""
        ls work/results/nucleosome/seqs/with/*.fa work/results/nucleosome/seqs/without/*.fa \
           work/results/nucleosome/control_seqs/with/*.fa work/results/nucleosome/control_seqs/without/*.fa \
           > work/results/nucleosome/all_fastas.txt
        rm -rf work/results/nucleosome/chunks && mkdir -p work/results/nucleosome/chunks
        split -n l/{threads} -d work/results/nucleosome/all_fastas.txt work/results/nucleosome/chunks/chunk_
        for c in work/results/nucleosome/chunks/chunk_*; do
            Rscript scripts/run_nupop.R "$c" > "${{c}}.log" 2>&1 &
        done
        wait
        """


# Per-locus occupancy scalars, paper-style binomial tests, telotron-vs-control comparison, and the
# species x architecture figure (element-removed insertion-site occupancy profiles).
rule nucleosome_figure:
    input:
        "work/results/nucleosome/.nupop.done",
    output:
        png="work/results/nucleosome/nucleosome_occupancy.png",
        per_locus="work/results/nucleosome/per_locus.tsv",
    conda:
        ENV
    shell:
        r"""
        python scripts/nucleosome_aggregate.py
        """


# Convenience target.
rule nucleosome:
    input:
        "work/results/nucleosome/nucleosome_occupancy.png",
        "work/results/nucleosome/per_locus.tsv",


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


# Host-gene-class characterisation (intron-rich/long; host vs DISJOINT non-host;
# per-genome contrast). Panel C GO enrichment is optional (degrades if the
# work/old/ archive table is absent).
rule telotron_gene_bias:
    input:
        arch=ARCH_TSV,
        refseq="data/raw/refseq/.done",
    output:
        "analysis/telotron_gene_bias.png",
    conda:
        ENV
    shell:
        r"""
        mkdir -p analysis
        python scripts/telotron_gene_bias.py
        """


# Per-intron telotron-rate logistic model (descriptive panels + cluster caveat).
rule telotron_per_intron:
    input:
        arch=ARCH_TSV,
        refseq="data/raw/refseq/.done",
    output:
        "analysis/telotron_per_intron.png",
    conda:
        ENV
    shell:
        r"""
        mkdir -p analysis
        python scripts/telotron_per_intron.py
        """


# Length distribution by (lineage x architecture) — BH-FDR corrected, single-MAG
# caveat in the output. Repo-relative paths (TELOTRON_ROOT override available).
rule length_distribution_by_arch:
    input:
        ARCH_TSV,
    output:
        summary="work/results/length_distribution_2026-06-07/summary_by_lineage_arch.tsv",
        linker="work/results/length_distribution_2026-06-07/linker_vs_nolinker.tsv",
        per_arch="work/results/length_distribution_2026-06-07/psw_vs_eimeria_per_arch.tsv",
        hist="work/results/length_distribution_2026-06-07/length_hist_by_lineage_arch.pdf",
    conda:
        ENV
    shell:
        r"""
        TELOTRON_ROOT="$(pwd)" python scripts/length_distribution_by_arch_lineage.py
        """


rule length_per_arm_figure:
    input:
        ARCH_TSV,
    output:
        "work/results/length_distribution_2026-06-07/per_arm_burst_length.pdf",
    conda:
        ENV
    shell:
        r"""
        TELOTRON_ROOT="$(pwd)" python scripts/length_per_arm_figure.py
        """


# Schematic mechanism figures — no data inputs (hardcoded labels distilled from
# the analyses); wired so they regenerate deterministically.
rule mechanism_diagrams:
    output:
        "analysis/telotron_mechanism_diagram.png",
        "analysis/proven_mechanism.png",
        "analysis/telotron_removal_fill_diagram.png",
    conda:
        ENV
    shell:
        r"""
        mkdir -p analysis
        python scripts/plot_mechanism_diagram.py
        python scripts/plot_proven_mechanism.py
        python scripts/plot_removal_fill_diagram.py
        """


# Expression-vs-telotron-presence figures. Gene-level coverage now comes from the
# reproducible rnaseq_gene_coverage rule (not /tmp). The locus-level panel in
# telotron_expr_presence additionally reads data/raw/rnaseq_splice_2026/
# per_locus_counts.tsv (still a manual build — see the external note below), so
# that rule is wired only when both species' coverage TSVs exist; the figure
# scripts skip the locus panel gracefully if the splice-counts file is absent.
rule telotron_expr_figures:
    input:
        eten="work/results/rnaseq/eten_gene_cov.tsv",
        necatrix="work/results/rnaseq/necatrix_gene_cov.tsv",
        arch=ARCH_TSV,
        refseq="data/raw/refseq/.done",
    output:
        "analysis/telotron_expr_presence.png",
        "analysis/telotron_expr_presence_necatrix.png",
    conda:
        ENV
    shell:
        r"""
        mkdir -p analysis
        python scripts/telotron_expr_presence.py
        python scripts/telotron_expr_final.py
        """


# Aggregate target for the wired analysis arm.
rule analysis_arm:
    input:
        rules.nucleosome_features.output,
        rules.nucleosome_withingene.output,
        rules.telotron_gene_bias.output,
        rules.telotron_per_intron.output,
        rules.length_distribution_by_arch.output,
        rules.length_per_arm_figure.output,
        rules.mechanism_diagrams.output,
        rules.telotron_expr_figures.output,


# ── One remaining external-input helper (documented, NOT wired) ─
#   telogator2 reference  (build_telogator2_ref.py)
#       needs: a hand-curated cap_survey.tsv + subtelomeres FASTA (T2T
#       deliverable). The per-species *.telogator2.fasta consumed by
#       telogator2_one are hand-built from its output.
#
# All previously "not wired" exploratory arms (Hi-C, ONT, age-ladder /
# cross-strain, 1-species expression test, proteome-BLAST telomerase,
# find_gene_deep_homology, branchpoint STREME/FIMO) were retired 2026-07-14
# under the Snakefile-driven / one-analysis-per-task principle. See
# memory/extreme_deslop_2026-07-14.md for the retirement record.
#
# detect_telomere_boundaries.py is fully CLI-parameterised and driven by the
# telomere_boundaries rule.
