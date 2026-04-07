# Strand-specific telomeric density across normalized intron coordinates

Extended Data Figure: Boundary telomeric density profiles for three telotron strand configurations, with intron coordinates normalized to [0, 1] and 50 bp exonic flanks shown in absolute coordinates.

Note: Among the 674 GT-AG telotrons in the primary MAG (93.5% of all 721 telotrons), 336 are converging and 328 are template-only (nearly equal split), with only 10 coding-strand-only and 14 GC-AG minor-class introns. All 52 telotrons with intron-free orthologs in related MAGs confirmed de novo insertion with ≥95% exonic identity; 71% carry target site duplications of 3–21 bp (mean 6.0 bp).

## The converging crossover

The most striking pattern emerges in the converging (head-to-head) configuration, which accounts for 54.7% of telotrons (n = 4,758). TTAGGG density on the gene strand peaks at ~3.5 near the 5′ splice site and decays smoothly toward the 3′ end, while CCCTAA density mirrors this from the opposite direction. The two signals cross over at approximately the normalized midpoint of the intron.

This crossover is a direct readout of the double-strand break repair mechanism. Telomerase extended the gene strand rightward from the 5′ break, depositing TTAGGG repeats readable in the FASTA 5′→3′ direction. Simultaneously, telomerase extended the anti-sense strand leftward from the 3′ break — those TTAGGG repeats, laid down 5′→3′ going right-to-left on the anti-sense strand, appear as CCCTAA when read on the gene strand. The midpoint crossover indicates that both extensions contributed roughly equal lengths on average.

The gradual sigmoidal shape of the transition reflects population-level variance: some introns received a longer gene-strand extension, others a longer anti-sense extension. Averaging across thousands of introns smooths this into the observed continuous crossover.

## Template-only: uniform CCCTAA

In the template-only configuration (40.8%, n = 21,408), nearly the entire intron shows CCCTAA on the gene strand with negligible TTAGGG signal. This indicates that only the anti-sense strand was extended by telomerase, with the gene strand subsequently filled in by DNA polymerase as the complement — hence CCCTAA throughout. The flat plateau across normalized positions 0.1–0.9 confirms a uniform repeat array rather than enrichment at one end.

## Coding-only: TTAGGG dominant

The coding-only configuration (3.6%, n = 2,152) is the mirror image: predominantly TTAGGG on the gene strand, confirming the gene strand was directly extended. A residual CCCTAA component visible in the 3′ half may represent a subset of introns where minor anti-sense extension also occurred, or degenerate repeat sequences matching CCCTAA rotations.

## Sharp splice-site boundaries

Across all three configurations, telomeric density drops from its plateau to near-zero within 1–2 bp at both the GT donor and AG acceptor sites. There is no gradual decay into exonic flanks — the arrays terminate precisely at annotated splice junctions. The faint signal in the ±50 bp exon regions (~0.03–0.05) represents the background rate of random hexamer matches in coding sequence, consistent with chance expectation for GC-rich hexamers.

## A spectrum of DSB repair outcomes

Together, the three profiles form a complete spectrum. Converging telotrons represent the canonical mechanism — both broken ends extended. Template-only represents cases where only the anti-sense 3′ overhang was engaged by telomerase. Coding-only is the rare inverse. The 11:1 ratio of template-only to coding-only suggests a bias toward anti-sense strand extension, potentially reflecting transcription-coupled asymmetry in 3′ overhang accessibility or telomerase recruitment.
