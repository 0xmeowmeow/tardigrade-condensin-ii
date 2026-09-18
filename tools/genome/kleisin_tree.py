"""Maximum-likelihood kleisin tree placing the tardigrade CAP-H2 candidates.

  ~/miniconda3/envs/tardigrade/bin/python tools/genome/kleisin_tree.py     (after caph2_profile.py)

68 sequences: the 33 CAP-H2 and 24 CAP-H profile seeds, the six tardigrade CAP-H2 candidates (one Paramacrobiotus
isoform), the three tardigrade CAP-H, and the cohesin kleisins human RAD21 (O60216) and yeast Scc1 (Q12158) as
outgroup. mafft --auto; columns with more than 50% gaps removed (631 of the alignment kept); IQ-TREE 3 with
ModelFinder (BIC) and 1,000 ultrafast bootstraps. Output in data/derived/condensin_tree/ (kleisins.treefile,
kleisins.iqtree). First run 19 September 2026 (as an interactive script; this file reproduces it).
"""
import pathlib, re, subprocess, urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[2]
BIN = pathlib.Path.home() / 'miniconda3/envs/tardigrade/bin'
DER = ROOT / 'data/derived'
OUT = DER / 'condensin_tree'
TARD = {'OWA52468.1': 'TARD__Hypsibius_exemplaris__OWA52468_CAPH2a', 'OWA51095.1': 'TARD__Hypsibius_exemplaris__OWA51095_CAPH2b',
        'GAU96458.1': 'TARD__Ramazzottius_varieornatus__GAU96458_CAPH2', 'GAU88883.1': 'TARD__Ramazzottius_varieornatus__GAU88883_CAPH2',
        'XP_055348578.1': 'TARD__Paramacrobiotus_metropolitanus__XP_055348578_CAPH2',
        'XP_055355123.1': 'TARD__Paramacrobiotus_metropolitanus__XP_055355123_CAPH2',
        'OQV23316.1': 'TARD__Hypsibius_exemplaris__OQV23316_CAPH', 'GAU89914.1': 'TARD__Ramazzottius_varieornatus__GAU89914_CAPH',
        'XP_055347164.1': 'TARD__Paramacrobiotus_metropolitanus__XP_055347164_CAPH'}
OUTGROUP = [('O60216', 'OUT__Homo_sapiens__RAD21'), ('Q12158', 'OUT__Saccharomyces_cerevisiae__SCC1')]

def fasta(p):
    d, k = {}, None
    for l in open(p):
        l = l.rstrip()
        if l.startswith('>'): k = l[1:]; d[k] = ''
        else: d[k] += l
    return d

def short(h, tag):
    m = re.search(r'\|([A-Z0-9]+)\|(\S+)', h); sp = re.search(r'OS=(\S+ \S+)', h)
    return f"{tag}__{sp.group(1).replace(' ', '_') if sp else h.split()[0]}__{m.group(1) if m else h.split()[0]}"

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    seqs = [(short(h, 'CAPH2'), s) for h, s in fasta(DER / 'condensin_caph2_seeds.fa').items()]
    seqs += [(short(h, 'CAPH'), s) for h, s in fasta(DER / 'condensin_caph_seeds.fa').items()]
    prot = {}
    for f in ('h.faa', 'r.faa', 'paramacrobiotus.faa'):
        prot.update({h.split()[0]: s for h, s in fasta(DER / 'condensin_check' / f).items()})
    seqs += [(n, prot[a]) for a, n in TARD.items()]
    for acc, n in OUTGROUP:
        seqs.append((n, ''.join(urllib.request.urlopen(f'https://rest.uniprot.org/uniprotkb/{acc}.fasta').read().decode().split('\n')[1:])))
    (OUT / 'kleisins.fa').write_text(''.join(f'>{n}\n{s}\n' for n, s in seqs))
    (OUT / 'kleisins.aln').write_text(subprocess.run([str(BIN / 'mafft'), '--auto', '--thread', '16', '--quiet', str(OUT / 'kleisins.fa')],
                                                     check=True, capture_output=True, text=True).stdout)
    aln = fasta(OUT / 'kleisins.aln')
    L = len(next(iter(aln.values())))
    keep = [i for i in range(L) if sum(s[i] == '-' for s in aln.values()) / len(aln) <= 0.5]
    (OUT / 'kleisins.trim.aln').write_text(''.join(f'>{k}\n{"".join(s[i] for i in keep)}\n' for k, s in aln.items()))
    print(len(seqs), 'sequences;', L, 'columns ->', len(keep))
    subprocess.run([str(BIN / 'iqtree3'), '-s', 'kleisins.trim.aln', '-m', 'MFP', '-B', '1000', '-T', '12', '--prefix', 'kleisins', '-redo'],
                   cwd=OUT, check=True)

if __name__ == '__main__':
    main()
