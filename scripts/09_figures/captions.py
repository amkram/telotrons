"""Figure captions as in the manuscript, written to figures/<stem>.caption.md."""

from pathlib import Path
import os

PACKAGE = Path(os.environ['TELO_PACKAGE'])

captions = {
    "Figure1_telotron_tree": (
        "**Figure 1. Distribution of telotrons across eukaryotes.** Filled circles mark telotron "
        "bearer groups and open squares mark controls, with genome counts in parentheses for groups "
        "of more than one genome. Every genome was screened with the published telomere motif of its "
        "lineage (supplementary table S1), and controls have no intron above 50% repeat coverage "
        "under that motif. Colors indicate lineages. The >50% counts include the >85% set. "
        "Percentages of absent introns use primary loci whose closest qualifying comparison either "
        "lacks an intron or carries an ordinary intron supported by annotation at the same aligned "
        "position and phase (2,115/2,255 overall, 93.8%). TSD percentages use all primary loci and "
        "are pooled per lineage."
    ),
    "Figure2_architecture_examples": (
        "**Figure 2. Telotron architecture and examples.** (A) Architectures of 4,331 classified "
        "telotrons in 25 bearer genomes. Bars show percentages within each lineage among classified "
        "loci, labeled with counts. (B) Each example outside Eimeria shows the closest validated "
        "ortholog lacking the intron. Eimeria shows a comparison with an ordinary intron. Exons are "
        "uppercase and introns lowercase. Gray letters show amino acids, numbers omitted bases and "
        "lines absent introns."
    ),
    "Figure3_ITS_density": (
        "**Figure 3. Interstitial telomeric sequence density in bearer genomes and controls.** Each "
        "point represents one of 25 bearers or 25 distinct controls. Circles indicate isolates, "
        "squares MAGs and colors lineages. Black bars show medians of 177.8 and 0.267 tracts per Mb "
        "(665-fold difference). Tracts overlapping annotated or projected gene spans or terminal "
        "telomeric regions were removed. Density uses the remaining A/C/G/T sequence. Lines connect "
        "seven selected lineage pairs (exact two-sided Wilcoxon W = 0, p = 0.0156). The y-axis is "
        "logarithmic above 0.1 and linear below."
    ),
    "FigureS1_unclassified_reasons": (
        "**Supplementary Figure S1. Architecture classification of telotrons.** Loci in the primary "
        "(A, >85% repeat coverage) and additional (B, >50% to ≤85%) tiers, grouped by the outcome of "
        "the architecture classifier. Classified loci (blue) carry one of the five observed "
        "architectures. Unclassified loci are grouped by the first rule they fail. The rule each "
        "locus failed is in supplementary table S2."
    ),
    "FigureS2_genome_tree": (
        "**Supplementary Figure S2. Telotrons and ortholog comparisons in individual genomes.** "
        "Expanded version of Fig. 1 with the same 25 bearers and 25 controls. The >50% counts include "
        "primary loci. Intron and TSD lengths are median [minimum–maximum] in bp. “Included in "
        "comparison” is the percentage of primary loci with a supported closest comparison showing an "
        "absent intron or an ordinary intron at the same position and phase. The column for absent "
        "introns gives absence among those comparisons. TSD percentages use all primary loci in each "
        "genome."
    ),
    "FigureS3_single_motif_introns": (
        "**Supplementary Figure S3. Introns composed of one repeat motif in bearers and controls.** "
        "Every annotated intron (≥30 bp) of the 25 bearer genomes and their 25 controls was scored "
        "against every primitive repeat motif of period 1 to 8 bp. An intron counts once when the "
        "motif covering most of it spans more than 85% of its length. Each point is one genome "
        "(circles isolates, squares MAGs, colors lineages). Boxes show the median and interquartile "
        "range, whiskers the extremes. The y-axis is logarithmic above 1 and linear below it. "
        "Telomeric motif: the lineage telomere motif or a motif within one substitution or indel of "
        "it. Counts for each bearer genome are in supplementary table S1."
    ),
    "FigureS1_architecture_by_tier": (
        "**Unnumbered figure. Repeat architectures by coverage tier.** (A) Primary loci (>85% "
        "coverage). (B) Additional loci (>50% to ≤85%). Bars show shares of classified loci in each "
        "lineage, labeled classified/all. Only linker architectures can be additional, since the "
        "others require >85% coverage. Rules and colors as in Fig. 2."
    ),
}
for stem, caption in captions.items():
    (PACKAGE / 'figures' / f'{stem}.caption.md').write_text(caption + '\n')
