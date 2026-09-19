# Tardigrades have a complete condensin II

Code and data for the preprint *Tardigrades have a complete condensin II: two divergent CAP-H2 kleisins in three eutardigrade genomes* (Meow-Ludo Meow-Meow, 2026; Zenodo, all versions [doi:10.5281/zenodo.22832510](https://doi.org/10.5281/zenodo.22832510); version 2 doi:10.5281/zenodo.22838790; version 3 [doi:10.5281/zenodo.22845658](https://doi.org/10.5281/zenodo.22845658)).

Hoencamp et al. (2021, *Science* 372:984) found no condensin accessory subunits in the tardigrade *Hypsibius dujardini* (now *H. exemplaris*), and counted it among the species that lack condensin II and fold their genomes in a Rabl-like (type-I) way. We searched three eutardigrade proteomes and found both condensin complexes complete, including two divergent CAP-H2 kleisins in each species. The results rest on these analyses:

- a profile HMM for CAP-H2, with positive controls and a shuffled-proteome null;
- a maximum-likelihood kleisin phylogeny;
- conserved gene order around each paralog;
- co-expression and developmental expression;
- a recomputation of Hoencamp et al.'s condensin II–architecture test;
- structure predictions of the kleisin regions bound to the condensin ring (version 2);
- the N-terminal regions predicted again with a longer SMC2 construct, five samples each (version 3; Fig. 4c).

The two *Ramazzottius varieornatus* CAP-H2 genes were first listed by van Hooff et al. (2025, *Cell Reports* 44:115855).

## Contents

| Path | What |
|---|---|
| `docs/papers/condensin-ii/` | Manuscript (Markdown and PDF), references, figures |
| `tools/genome/condensin_check.py` | Reciprocal blastp/PSI-BLAST search for all eight condensin subunits |
| `tools/genome/caph2_profile.py` | CAP-H2 and CAP-H profile HMMs, Pfam domain scan, shuffled null, co-expression |
| `tools/genome/kleisin_tree.py` | Kleisin alignment and IQ-TREE maximum-likelihood tree |
| `tools/genome/caph2_synteny.py` | Gene-order anchors around each CAP-H2 locus, against all ortholog pairs and random pairs |
| `tools/genome/hoencamp_fisher.py` | Hoencamp et al. Table S6 transcribed; Fisher's exact test as published and with the tardigrade reclassified |
| `tools/genome/caph2_structure_inputs.py`, `caph2_structure_analysis.py` | Boltz-2 inputs (kleisin regions with SMC head constructs, controls) and the ipTM/pLDDT/interface summary |
| `tools/genome/condensin_figdata.py`, `condensin_figures.py` | Figure data and figures |
| `data/derived/condensin_caph2_seeds.fa`, `condensin_caph_seeds.fa` | Profile seed sequences (UniProt; 33 CAP-H2, 24 CAP-H; no tardigrades) |
| `data/derived/condensin_check.tsv` | Reciprocal-search results per subunit and species |
| `data/derived/condensin_caph2_profile.tsv` | Profile hits, Pfam domains, null, co-expression summary |
| `data/derived/condensin_tree/` | Kleisin sequences, trimmed alignment, tree and IQ-TREE report |
| `data/derived/condensin_synteny.tsv` | Anchors for every pair of CAP-H2 loci |
| `data/derived/caph2_structure/` | Boltz-2 inputs, construct boundaries, `summary.tsv`, and the best model of each main prediction (`models/*.pdb`); version 3's longer-SMC2 N-terminal run in `inputs_n2/`, `constructs_n2.json`, `summary_n2.tsv` (`caph2_structure_analysis.py --n2`) and `models/n2_*.pdb`, `models/ctl_n2shuf_*.pdb` |
| `data/derived/condensin_figdata.json` | Domain positions and expression profiles plotted in the figures |
| `data/derived/genome.sqlite.gz` | Expression inputs: Yoshida et al. 2017 TPM (48 samples) and Levin et al. 2016 single-embryo counts (62 embryos), with the transcript-to-gene join onto nHd_3.1 gene models |

## Reproducing

Python 3.12 with numpy, scipy, matplotlib, biopython and pyhmmer; BLAST+, MAFFT, IQ-TREE 3 and DIAMOND on the path of the `BIN` directory set at the top of each script (edit it for your system).

1. Download the public inputs into `data/raw/`:
   - *H. exemplaris* nHd_3.1: GCA_002082055.1 (protein FASTA and GFF);
   - *R. varieornatus* Rvar_4.0: GCA_001949185.1;
   - *P. metropolitanus*: GCF_019649055.1;

   all from NCBI Datasets or the NCBI FTP site, in the file names the scripts use.
2. `gunzip -k data/derived/genome.sqlite.gz`.
3. Run, in this order:
   1. `condensin_check.py`;
   2. `caph2_profile.py`;
   3. `kleisin_tree.py`;
   4. `caph2_synteny.py`;
   5. `hoencamp_fisher.py`;
   6. `condensin_figdata.py`;
   7. `condensin_figures.py`.

The figures in this repository regenerate pixel-identical from the included data. The tree was inferred once with IQ-TREE 3.1.3; its report is included.

## Sources of data

- **Proteomes and annotations:** NCBI (Yoshida et al. 2017, *PLoS Biology*; Hashimoto et al. 2016, *Nature Communications*; Hara et al. 2021, *Open Biology*).
- **Expression:**
  - GEO GSE94295 (Yoshida et al. 2017);
  - GEO GSE70185 and TSA GBZR01 (Levin et al. 2016, *Nature*).
- **Pfam models:** InterPro.
- **Seed sequences:** UniProt.
- **Condensin II status and architecture types:** Hoencamp et al. 2021, Supplementary Tables S3 and S6.

## Licence and citation

Code: MIT. Data and figures produced here: CC BY 4.0. Third-party data keep their own terms. Please cite the preprint (see `CITATION.cff`).

The analyses, figures and first draft were produced with an AI assistant (Claude, Anthropic) working under the author's direction; the author checked the results and takes responsibility for the content.
