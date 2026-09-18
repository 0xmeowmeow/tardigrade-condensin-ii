---
title: "Tardigrades have a complete condensin II: two divergent CAP-H2 kleisins in three eutardigrade genomes"
author:
  - "Meow-Ludo Meow-Meow (meow@meow-meow.io)"
date: "Preprint, version 1, 19 September 2026 · doi:10.5281/zenodo.22832511"
abstract: |
  Condensin II shapes interphase genome architecture: species that lack it tend to fold their chromosomes in a Rabl-like, type-I configuration, while species that have it form chromosome territories (type II). In the study that established this, the tardigrade *Hypsibius dujardini* (now *H. exemplaris*) was the only species of 24 in which neither condensin's accessory subunits could be found. It was counted among the type-I species that lack condensin II, and the authors left open whether the subunits had diverged beyond recognition. Here we searched the annotated proteomes of three eutardigrades and found both complexes complete. Reciprocal BLAST recovers condensin I and the condensin II HEAT subunits CAP-D3 and CAP-G2. A profile hidden Markov model built from 33 animal CAP-H2 sequences finds the missing kleisin: each species has two CAP-H2 genes. Both carry the Pfam CNDH2_N and CNDH2_C domains, score against CAP-H2 and not CAP-H, and form a single tardigrade clade within CAP-H2 in a maximum-likelihood kleisin tree (ultrafast bootstrap 94 for the placement, 100 for the clade). Conserved gene order dates the duplication to before the split between *Hypsibius* and *Paramacrobiotus*. Copy a sits at an orthologous locus in all three species, and the b copies of *Hypsibius* and *Paramacrobiotus* share a conserved neighbourhood. In *H. exemplaris*, CAP-H2 is co-expressed with CAP-D3, CAP-G2, SMC2 and SMC4 across 48 developmental and adult samples; its correlations with them rank in the top 0.1–0.3% of 14,851 expressed genes. In single embryos, all condensin II subunits are expressed through the first ~46 hours after laying and fall to near zero thereafter, while CAP-H persists. The two *Ramazzottius varieornatus* CAP-H2 genes were already listed in a pan-eukaryotic survey (van Hooff et al. 2025); the *H. exemplaris* and *Paramacrobiotus* genes, the correction of the species record behind the architecture study, the paralog history (from phylogeny and gene order) and the expression data are new here. With the tardigrade reclassified, the condensin II–architecture association still holds (Fisher's exact p = 0.019, from 0.007). The tardigrade becomes the seventh species with a complete condensin II and a type-I call, a call that rests on a single manually annotated feature in whole-animal Hi-C. Early embryos, where condensin II is expressed and chromosomes remain individualised through interphase, are the place to test whether tardigrade genome architecture changes with developmental stage.
geometry: margin=2.3cm
fontsize: 10pt
linestretch: 1.15
colorlinks: true
header-includes:
  - \usepackage{booktabs}
  - \usepackage{float}
  - \floatplacement{figure}{H}
  - \usepackage{newunicodechar}
  - \newunicodechar{−}{\ensuremath{-}}
  - \newunicodechar{≤}{\ensuremath{\leq}}
  - \newunicodechar{×}{\ensuremath{\times}}
  - \newunicodechar{ρ}{\ensuremath{\rho}}
  - \newunicodechar{‡}{\ensuremath{\ddagger}}
  - \newunicodechar{†}{\ensuremath{\dagger}}
  - \newunicodechar{§}{\S}
---

# Introduction

Condensin complexes compact chromosomes for mitosis. Animals usually have two, which share the SMC2–SMC4 ATPase but differ in their kleisin (CAP-H in condensin I, CAP-H2 in condensin II) and their HEAT-repeat subunits (CAP-D2 and CAP-G; CAP-D3 and CAP-G2) [@hirano2012]. Condensin II also acts beyond mitosis. Hoencamp et al. [-@hoencamp2021] used Hi-C on 24 species across the tree of life and found two architecture types. Type-I genomes fold in a Rabl-like configuration, with centromeres and telomeres clustered; type-II genomes form chromosome territories. Condensin II status tracked the type: "Eight species lacked one or more condensin II subunit(s) (table S6) and exhibited Rabl-like features (table S3)" [@hoencamp2021, p. 986]. Deleting CAP-H2 in human cells was enough to shift them toward type-I features.

The tardigrade *Hypsibius dujardini* was one of those eight species. It stood apart even among them: SMC2 and SMC4 were identified, but "could not identify any of the accessory subunits of condensin I nor II" [@hoencamp2021, Supplementary Materials p. 11]. The authors named two explanations: the subunits "could have evolved beyond recognition with our methods", or the species had lost both complexes. They concluded that the question "will need to be investigated further". Their search used tblastn and tblastx against the nucleotide and whole-genome-shotgun databases, and their SMC hits come from the early LMYF01 and LRSR01 assemblies (Table S6). The strain they sequenced, the Sciento Z151 laboratory culture, was later described as *H. exemplaris* [@gasiorek2018]. We use that name below.

Losing condensin I as well as condensin II would be extraordinary for an animal. Kleisins evolve fast, and three annotated tardigrade proteomes are now available [@yoshida2017; @hara2021], so we reinvestigated the question.

# Results

## Condensin I and the condensin II HEAT subunits are present

We searched human reviewed sequences of all eight condensin subunits against the proteomes of *H. exemplaris* (nHd_3.1), *Ramazzottius varieornatus* (Rvar_4.0) and *Paramacrobiotus metropolitanus* (GCF_019649055.1), with blastp and, for the kleisins, five rounds of PSI-BLAST. We kept a hit only if it was a reciprocal best hit against the whole human reviewed proteome (Table 1).

SMC2, SMC4, CAP-D2, CAP-D3 and CAP-G2 are reciprocal best hits in all three species. CAP-G is a weak reciprocal best hit (E ~4 × 10^−6^) in *H. exemplaris* and *Paramacrobiotus*. CAP-H was found by PSI-BLAST in *Paramacrobiotus* (XP_055347164.1; reciprocal with human NCAPH). The *H. exemplaris* (OQV23316.1) and *R. varieornatus* (GAU89914.1) CAP-H proteins were found through it: both hit it as their best match, and a CAP-H profile built from 24 animals recovers exactly these three proteins and no others. So condensin I is complete, and condensin II has its HEAT subunits. Only CAP-H2 was missing, as it was from blastp and five rounds of PSI-BLAST in all three species.

## A CAP-H2 profile finds two genes in each species

Kleisins are among the least conserved condensin subunits, so we built a profile HMM from 33 animal CAP-H2 sequences (UniProt, one per species, no tardigrades). We searched it against the three proteomes alongside a profile of 24 CAP-H sequences. A true CAP-H2 should score against the CAP-H2 profile and not the CAP-H one.

Each species has two such proteins, and none scores against CAP-H (Table 2). Their E-values run from 4 × 10^−11^ to 2 × 10^−20^, and every one scores E = 1 on the CAP-H profile. All six carry the Pfam CAP-H2 N-terminal domain (CNDH2_N, PF06278; E 5 × 10^−6^ to 3 × 10^−9^) near their start and the C-terminal domain (CNDH2_C, PF16869; E 7 × 10^−11^ to 9 × 10^−12^) near their end. They lack the CAP-H domain Cnd2 (PF05786), which all three tardigrade CAP-H proteins carry (Fig. 1). The middle CAP-H2 domain (CNDH2_M) was not detected in any of them.

We ran three controls. Human NCAPH2 scores E 2 × 10^−216^ on the CAP-H2 profile and E 0.7 on the CAP-H profile. The CAP-H profile finds exactly the known CAP-H in each tardigrade. A version of the *H. exemplaris* proteome with every sequence shuffled gives no hit at E ≤ 10 in three replicates.

![**Domain architecture of the tardigrade kleisins.** Pfam domain envelopes (hmmscan, independent E ≤ 10^−3^) on human CAP-H2 and CAP-H and on the tardigrade proteins. The six CAP-H2 candidates share the CNDH2_N and CNDH2_C domains of human CAP-H2; the three tardigrade CAP-H proteins carry Cnd2 instead. The tardigrade CAP-H proteins are ~1,050 residues long, with an insertion between their two Cnd2 regions.](figs/fig1-domains.pdf){width=100%}

## The tardigrade genes are CAP-H2, from an old duplication

To place the candidates, we aligned them with the 57 profile seeds, the three tardigrade CAP-H proteins, and human RAD21 and yeast Scc1 (cohesin kleisins) as outgroup. We then inferred a maximum-likelihood tree (Q.YEAST+F+R5, 631 columns, 1,000 ultrafast bootstraps; Fig. 2).

All six candidates form one clade (ultrafast bootstrap, UFBoot, 100) inside CAP-H2 (UFBoot 94 for CAP-H2 plus the tardigrade clade). The three tardigrade CAP-H proteins group with CAP-H (100). Within the tardigrade CAP-H2 clade, one paralog ("a": *H. exemplaris* OWA52468.1, *R. varieornatus* GAU96458.1, *Paramacrobiotus* XP_055355123.1) is supported at 100. A second ("b": OWA51095.1, GAU88883.1) is supported at 99, and the two group together at 98. The second *Paramacrobiotus* gene (XP_055348578.1) falls outside both.

So the duplication happened before *Hypsibius* and *Ramazzottius* diverged. Whether *Paramacrobiotus* shares it, with a fast-evolving b copy, or duplicated separately is not resolved by this tree; gene order resolves it (next section). The tardigrade clade sits on a long branch (1.39 substitutions per site) and branches as sister to all other animal CAP-H2 rather than near arthropods or nematodes, which may reflect long-branch attraction.

![**Maximum-likelihood tree of CAP-H2, CAP-H and tardigrade kleisins**, rooted on the cohesin kleisins RAD21 and Scc1. Blue: CAP-H2 seeds; green: CAP-H seeds; orange: tardigrade proteins. Numbers are ultrafast bootstrap values of 70 or more. (a) and (b) mark the two *H. exemplaris* paralogs.](figs/fig2-tree.pdf){width=92%}

## Gene order shows the duplication is older than the split with *Paramacrobiotus*

To date the duplication independently of sequence, we compared gene order around each CAP-H2 locus (Table 3). For each pair of loci in two species, we counted the genes within 20 genes of one locus whose reciprocal best hit (DIAMOND, whole proteomes) lies within 20 genes of the other. We compared that count with two references: every orthologous gene pair between the same two species, and 20,000 random gene pairs.

Copy a sits at an orthologous locus in all three species. *H. exemplaris* OWA52468.1 and *R. varieornatus* GAU96458.1 share 8 neighbours, and each shares 2 with *Paramacrobiotus* XP_055355132.1. Random gene pairs share 0.03–0.06 on average, so p = 0.001–0.006.

Copy b of *H. exemplaris* (OWA51095.1) shares 7 neighbours with the second *Paramacrobiotus* gene, XP_055348578.1 (p = 0.002). That gene therefore sits at the b locus even though the tree could not place it. No a locus shares any neighbour with a b locus.

These counts are ordinary for true orthologs. Across all orthologous pairs, the median is 3–5 shared neighbours, and 27–33% share none. The *R. varieornatus* b copy (GAU88883.1) is one of those: it shares no neighbours with either b locus, although it is the reciprocal best hit of both b genes and groups with the *H. exemplaris* b copy in the tree (UFBoot 99). Its neighbourhood has probably been rearranged since the split.

So both copies were present in the common ancestor of *Hypsibius*, *Ramazzottius* and *Paramacrobiotus*. That ancestor predates the split between the superfamilies Hypsibioidea and Macrobiotoidea, which is deep within the eutardigrade order Parachela. The pairings from gene order and from sequence agree: every a–a and b–b pair is a reciprocal best hit, and no a–b pair is.

## CAP-H2 is expressed with the other condensin II subunits

Both *H. exemplaris* CAP-H2 genes (BV898_16921 and BV898_15594) are expressed. Across the 48 RNA-seq samples of Yoshida et al. [-@yoshida2017] (eggs on days 1–5, juveniles on days 1–7 and active and dried adults, in triplicate), BV898_16921 correlates with CAP-G2 (Spearman ρ 0.90), CAP-D3 (0.89), SMC4 (0.89) and SMC2 (0.87). Against all 14,851 other expressed genes, each of these correlations ranks in the top 0.1–0.3%; its correlation with CAP-H (0.67) ranks in the top 15% (Fig. 3c). BV898_15594 follows the same pattern less tightly: CAP-G2 0.85 (top 1.2%) and CAP-H 0.66 (top 20%).

The subunits peak together on egg days 1–2 and in adults, and are low in juveniles and in eggs from day 3 (Fig. 3a). Dividing cells express many genes together, so co-expression supports the pairing with condensin II but does not prove it.

In 62 single embryos [@levin2016], the condensin II subunits and SMC2 and SMC4 are expressed through the first ~46 hours after laying (median 76–223 counts per million) and drop to near zero after ~47 hours (Fig. 3b). CAP-H stays at a median of 47 counts per million. The times here use the embryo time axis corrected to ~98 hours of development; the source workbook's minutes are halved relative to real time. The drop coincides with the largest transcriptional change in *H. exemplaris* development, between egg days 2 and 3. It also falls at a library-batch boundary in this dataset, so its exact timing needs care.

![**Expression of the *H. exemplaris* condensin genes.** (a) Mean TPM per condition in Yoshida et al. 2017, scaled to each gene's maximum. E1–E5: egg days; J1–J7: juvenile days; act and tun: active and dried adults, in two sample sizes (10k and 30 animals). (b) Counts per million in 62 single embryos (Levin et al. 2016) against hours after laying; lines are 5-embryo running means. BV898_15594 has no transcript in this assembly. (c) Spearman correlation of CAP-H2 BV898_16921 with every expressed gene across the 48 Yoshida samples; vertical lines mark the condensin subunits.](figs/fig3-expression.pdf){width=100%}

## The condensin II–architecture association with the tardigrade reclassified

We transcribed Hoencamp et al.'s Table S6 (subunits found) and Table S3 (architecture type) for the 23 non-holocentric species. This reproduces their association. Of the species with a complete condensin II, 9 are type II and 6 type I; all 8 species lacking a subunit are type I (Fisher's exact test, two-sided p = 0.0072; the source reports p < 0.05).

Moving *H. exemplaris* to "complete" gives 9 type-II and 7 type-I species with a complete condensin II, against 0 and 7 among those lacking it (p = 0.019). The association survives. The tardigrade joins the frog, sea squirt, coral, comb jelly, peanut and wheat as species with a complete condensin II and a type-I call.

That call is weak. In Table S3, the tardigrade's four automatic feature calls (territories, a feature scored as "A", centromere clustering and telomere clustering) are all negative. Its type-I assignment rests on one manually annotated feature, telomere clustering. Its Hi-C also came from whole animals (GSM5182729).

# Discussion

Tardigrades are not an exception to the condensin repertoire of animals. Both complexes are complete in three eutardigrades, and condensin II has two kleisins, both diverged far enough to escape pairwise searches from human sequences. The reason the architecture study missed them is clear from its methods: nucleotide-level tblastn and tblastx with E ≤ 10^−10^, against an early *H. dujardini* assembly. Of all condensin subunits, only CAP-H2 needed a profile search to be found in tardigrades.

**Prior art.** Van Hooff et al. [-@vanhooff2025] surveyed SMC complexes across eukaryotes with profile-based orthology and phylogenies. Their phylogenetic profiles list *R. varieornatus* with two CAP-H2 genes, one CAP-G2 and three CAP-D3 copies. Their two CAP-H2 sequences are identical to our GAU96458.1 and GAU88883.1. Their species set does not include *H. exemplaris* or *Paramacrobiotus*, and their text does not discuss tardigrades or the architecture study's tardigrade record. UniProt's automatic annotation names three of our six proteins "Condensin-2 complex subunit H2 C-terminal domain-containing protein". The genes were therefore findable. What is new here:

- the genes in the species whose record underlies the published association;
- the duplication history;
- the co-expression and developmental expression;
- the consequences for the association itself.

**What the correction changes.** The condensin II–architecture association rests on species that lack condensin II all being type I, and the tardigrade was the strongest-looking member of that group, lacking every accessory subunit. Reclassified, it weakens the test (p 0.007 to 0.019) without overturning it. Van Hooff et al. already noted other exceptions: *Plasmodium falciparum* and *Dictyostelium discoideum* have condensin II and a Rabl-like organisation, and the moss *Physcomitrium patens* lacks it and forms territories [@vanhooff2025].

The tardigrade adds an animal case with an unusual twist: its condensin II subunits are strongly developmentally regulated. They are high in early embryos and in adults, and near zero in late embryos and juveniles (Fig. 3). An architecture measured in whole adults, whose tissues are mostly post-mitotic, may not describe every stage.

**A testable connection.** Papell et al. [-@papell2026] found that in early *H. exemplaris* embryos, chromosomes remain individualised through interphase, each in its own fully or partly separate lamin-lined compartment. Their study does not discuss condensins. Condensin II promotes chromosome individualisation and territories elsewhere [@hoencamp2021]. In *H. exemplaris* it is expressed in exactly the stages Papell et al. imaged, and it falls silent after ~47 hours. We propose that tardigrade genome architecture is stage-dependent: territorial, condensin II-driven compartments in early embryos, and a weakly Rabl-like configuration in adult tissue.

Three experiments would test this:

1. Hi-C on staged embryos before and after the ~47-hour drop.
2. RNAi of the CAP-H2 paralogs, which works in *H. exemplaris* embryos [@tenlen2013], followed by imaging of chromosome individualisation.
3. Localisation of CAP-H2 in early embryos.

**Limitations.** We used only annotated proteomes. Gene order places the duplication before the split between Hypsibioidea and Macrobiotoidea; a genome from another order (Apochela) or from a heterotardigrade would date it further back. The only heterotardigrade data we could examine, an *Echiniscus testudo* transcriptome (GDAL01), contains only a fragment of SMC4 among the condensin genes. The CNDH2_M domain was not detected. This may reflect divergence or a genuine loss of that region, which matters for binding to CAP-G2. Co-expression is correlative. We did not reanalyse the Hi-C data.

# Methods

**Proteomes.** We used *H. exemplaris* nHd_3.1 (GCA_002082055.1; 20,853 proteins) [@yoshida2017], *R. varieornatus* Rvar_4.0 (GCA_001949185.1; 23,007 proteins) [@yoshida2017; @hashimoto2016] and *P. metropolitanus* (GCF_019649055.1, NCBI RefSeq annotation) [@hara2021].

**Reciprocal searches.** Human reviewed UniProt sequences of SMC2 (O95347), SMC4 (Q9NTJ3), NCAPD2 (Q15021), NCAPG (Q9BPX3), NCAPH (Q15003), NCAPD3 (P42695), NCAPG2 (Q86XI2) and NCAPH2 (Q6IBW4) were used as queries: blastp (BLAST+) against each proteome, plus five PSI-BLAST rounds for the kleisins. A best hit was accepted when its best hit back against the human reviewed proteome was the query subunit.

**Profile search.** We built the profiles from 33 animal CAP-H2 and 24 animal CAP-H sequences from UniProt (one per species, no tardigrades). Each set was aligned with MAFFT (--auto) [@katoh2013], and profiles were built with pyHMMER [@larralde2023], a Python binding to HMMER3 [@eddy2011]. We then ran hmmsearch against each proteome (E ≤ 10^−3^). Hits on the CAP-H2 profile and not the CAP-H profile were scanned with the Pfam models PF06278, PF16858, PF16869 and PF05786 [@mistry2021]. The null model was the *H. exemplaris* proteome with residues shuffled within each sequence (three replicates).

**Phylogeny.** We aligned 68 sequences with MAFFT (--auto) and removed columns with more than 50% gaps, leaving 631 columns. A maximum-likelihood tree was inferred with IQ-TREE 3.1.3 [@minh2020], with ModelFinder selecting Q.YEAST+F+R5 [@minh2021] by BIC [@kalyaanamoorthy2017] and 1,000 ultrafast bootstraps [@hoang2018]. The tree was rooted on human RAD21 and yeast Scc1.

**Gene order.** Gene order came from the GenBank and RefSeq GFF files, with the longest protein per gene (ties broken by accession). Orthologs were reciprocal best hits between whole proteomes, from DIAMOND blastp (--more-sensitive, E ≤ 10^−5^, best target only). Each pair of CAP-H2 loci was scored by its anchors: genes within 20 genes of one locus whose ortholog lies within 20 genes of the other. We compared each score with the anchors of every orthologous gene pair between the two species, and with 20,000 random gene pairs.

**Expression.** Yoshida et al. [-@yoshida2017] provide transcript TPM for 48 samples (16 conditions × 3). We summed transcripts per gene and kept genes with TPM > 1 in at least 6 samples (14,852 genes). We then computed Spearman correlations between all genes and each CAP-H2 gene, and report each subunit's rank among all genes.

For the embryos, we used Levin et al.'s [-@levin2016] single-embryo counts for 62 *H. exemplaris* embryos, mapped to nHd_3.1 gene models by transcript alignment. We converted them to counts per million per embryo. Embryo times are the workbook's minutes × 2. This matches each embryo's profile to Yoshida's daily egg samples and gives ~98 hours of development.

**Fisher's exact test.** We transcribed Hoencamp et al. Table S6 (condensin II complete if CAP-H2, CAP-G2 and CAP-D3 are all listed) and the manual architecture type (Table S3). *C. elegans* (holocentric, no type) was excluded, as in the source. We used scipy.stats.fisher_exact.

**Code and data availability.** All inputs are public (accessions above). Scripts:

- `condensin_check.py`, `caph2_profile.py`, `kleisin_tree.py`, `caph2_synteny.py`, `condensin_figdata.py`, `condensin_figures.py` and `hoencamp_fisher.py` (in `tools/genome/`);
- seed sets `condensin_caph2_seeds.fa` and `condensin_caph_seeds.fa`;
- the alignment and tree.

Code, derived data, alignments and trees: <https://github.com/0xmeowmeow/tardigrade-condensin-ii>. This preprint: Zenodo, [doi:10.5281/zenodo.22832511](https://doi.org/10.5281/zenodo.22832511). Each claim is recorded with its source and page in the project's knowledge graph.

# Acknowledgements and use of AI

The analyses, figures and first draft of this manuscript were produced with an AI assistant (Claude, Anthropic) working under the author's direction. The author checked the results and takes responsibility for the content.

# Table 1. Condensin subunits in three tardigrade proteomes

| Subunit | *H. exemplaris* | *R. varieornatus* | *P. metropolitanus* |
|---|---|---|---|
| SMC2 | OQV15428.1 | GAV04166.1 | XP_055343015.1 |
| SMC4 | OWA52980.1 | GAV06634.1 | XP_055342529.1 |
| CAP-D2 | OQV13474.1 | GAU88108.1 | XP_055334720.1 |
| CAP-G | OWA50980.1 (E 4e-6) | not found by blastp\* | XP_055332312.1 (E 4e-6) |
| CAP-H | OQV23316.1† | GAU89914.1† | XP_055347164.1 (PSI-BLAST) |
| CAP-D3 | OQV21910.1 | GAU88810.1 | XP_055335315.1 |
| CAP-G2 | OQV13088.1 | GAU94400.1 | XP_055351168.1 |
| CAP-H2 | OWA52468.1, OWA51095.1‡ | GAU96458.1, GAU88883.1‡ | XP_055355132.1§, XP_055348578.1‡ |

\* Listed as present by van Hooff et al. (2025). † Found through the *Paramacrobiotus* CAP-H and the CAP-H profile. ‡ CAP-H2 profile only (Table 2). § One gene with isoforms XP_055355123/132/139.1. All other entries are reciprocal best blastp hits with the human subunit.

# Table 2. CAP-H2 profile hits

| Species | Protein | Length (aa) | E, CAP-H2 profile | E, CAP-H profile | CNDH2_N | CNDH2_C |
|---|---|---:|---:|---:|---|---|
| *H. exemplaris* | OWA52468.1 | 719 | 1.0e-15 | 1 | 5e-8 (54–164) | 2e-10 (579–685) |
| *H. exemplaris* | OWA51095.1 | 646 | 4.7e-14 | 1 | 1e-6 (9–122) | 8e-11 (512–617) |
| *R. varieornatus* | GAU96458.1 | 702 | 1.9e-20 | 1 | 3e-8 (31–148) | 9e-12 (547–652) |
| *R. varieornatus* | GAU88883.1 | 617 | 3.8e-12 | 1 | 5e-6 (14–125) | 7e-11 (498–608) |
| *P. metropolitanus* | XP_055348578.1 | 626 | 4.3e-13 | 1 | 3e-9 (12–119) | 1e-11 (493–623) |
| *P. metropolitanus* | XP_055355132.1 | 656 | 3.9e-11 | 1 | 7e-7 (18–134) | 6e-12 (507–622) |

# Table 3. Gene order around the CAP-H2 loci

Anchors: genes within 20 genes of one locus whose reciprocal best hit lies within 20 genes of the other. The p-value is the share of 20,000 random gene pairs with at least as many anchors. Pairs not shown have 0 anchors.

| Locus 1 | Locus 2 | Anchors | p (random pairs) | Reciprocal best hit |
|---|---|---:|---:|---|
| *H. exemplaris* a, OWA52468.1 | *R. varieornatus* a, GAU96458.1 | 8 | 0.001 | yes |
| *H. exemplaris* a, OWA52468.1 | *P. metropolitanus* a, XP_055355132.1 | 2 | 0.005 | yes |
| *R. varieornatus* a, GAU96458.1 | *P. metropolitanus* a, XP_055355132.1 | 2 | 0.006 | yes |
| *H. exemplaris* b, OWA51095.1 | *P. metropolitanus* b, XP_055348578.1 | 7 | 0.002 | yes |
| *H. exemplaris* b, OWA51095.1 | *R. varieornatus* b, GAU88883.1 | 0 | 1 | yes |
| *R. varieornatus* b, GAU88883.1 | *P. metropolitanus* b, XP_055348578.1 | 0 | 1 | yes |

Orthologous gene pairs in general keep a median of 3–5 anchors, and 27–33% keep none.

# References
