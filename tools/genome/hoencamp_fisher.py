"""Recompute Hoencamp et al. 2021's condensin II status x architecture type test with the tardigrade reclassified.

  python3 tools/genome/hoencamp_fisher.py

Transcribed from Hoencamp et al. 2021 Supplementary Materials, Table S6 (PDF p. 50: subunit accessions, blank = not
found; last column manual_type from Table S3). C. elegans (holocentric, n/a) is left out, as in the source. The
source states only 'dependent variables (p<0.05)' for its Fisher's exact test.
"""
from scipy.stats import fisher_exact

# species: (condensin II complete in Table S6, manual architecture type)
S6 = {
    'Macropus eugenii': (True, 'II'), 'Gallus gallus': (True, 'II'), 'Python bivittatus': (True, 'II'),
    'Xenopus laevis': (True, 'I'), 'Pygocentrus nattereri': (True, 'II'), 'Chiloscyllium punctatum': (True, 'II'),
    'Lethenteron camtschaticum': (True, 'II'), 'Ciona intestinalis/robusta': (True, 'I'),
    'Branchiostoma lanceolatum': (True, 'II'), 'Strongylocentrotus purpuratus': (True, 'II'),
    'Aedes aegypti': (False, 'I'), 'Culex quinquefasciatus': (False, 'I'), 'Drosophila melanogaster': (False, 'I'),
    'Hypsibius dujardini': (False, 'I'), 'Clonorchis sinensis': (False, 'I'), 'Aplysia californica': (True, 'II'),
    'Cristatella mucedo': (False, 'I'), 'Acropora millepora': (True, 'I'), 'Pleurobrachia bachei': (True, 'I'),
    'Agaricus bisporus': (False, 'I'), 'Saccharomyces cerevisiae': (False, 'I'), 'Arachis hypogaea': (True, 'I'),
    'Triticum aestivum': (True, 'I'),
}

def table(d):
    return [[sum(1 for c, t in d.values() if c == comp and t == typ) for typ in ('II', 'I')] for comp in (True, False)]

def main():
    out = []
    for label, d in [('as published', S6), ('tardigrade condensin II complete (this study)', {**S6, 'Hypsibius dujardini': (True, 'I')})]:
        t = table(d)
        out.append(f"{label}: complete [II, I] = {t[0]}, incomplete [II, I] = {t[1]}; Fisher two-sided p = {fisher_exact(t)[1]:.4f}, "
                   f"one-sided p = {fisher_exact(t, alternative='greater')[1]:.4f}")
    print('\n'.join(out))
    return out

if __name__ == '__main__':
    main()
