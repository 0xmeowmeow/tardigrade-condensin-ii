"""Profile-HMM search for the condensin II kleisin CAP-H2 in three tardigrade proteomes.

  ~/miniconda3/envs/tardigrade/bin/python tools/genome/caph2_profile.py   (run condensin_check.py first)

condensin_check.py found every condensin subunit except CAP-H2 by blastp and PSI-BLAST. Kleisins evolve fast, so
here: a profile built from 33 CAP-H2 sequences of other animals (UniProt, one per species, no tardigrades;
data/derived/condensin_caph2_seeds.fa), aligned with mafft, is searched against the three proteomes. A CAP-H
profile from 24 sequences (condensin_caph_seeds.fa) is the positive control and the discriminator: a true CAP-H2
should score on the CAP-H2 profile and not the CAP-H one. Each hit is then scanned with the Pfam CNDH2_N
(PF06278), CNDH2_M (PF16858), CNDH2_C (PF16869) and Cnd2 (PF05786, CAP-H) models. Null: the H. exemplaris proteome
with every sequence shuffled, three times. Then co-expression of the H. exemplaris genes with the other subunits
across Yoshida 2017's 48 samples (Spearman on TPM, ranked against every expressed gene). Output:
data/derived/condensin_caph2_profile.tsv.
"""
import gzip, pathlib, random, subprocess, urllib.request
import pyhmmer
from pyhmmer.easel import Alphabet, MSAFile, SequenceFile, TextSequence, DigitalSequenceBlock
from pyhmmer.plan7 import Builder, Background, HMMFile

ROOT = pathlib.Path(__file__).resolve().parents[2]
BIN = pathlib.Path.home() / 'miniconda3/envs/tardigrade/bin'
WORK = ROOT / 'data/derived/condensin_check'
DER = ROOT / 'data/derived'
PROTEOMES = {'H. exemplaris': 'h.faa', 'R. varieornatus': 'r.faa', 'Paramacrobiotus (GCF_019649055.1)': 'paramacrobiotus.faa'}
PFAM = {'PF06278': 'CNDH2_N', 'PF16858': 'CNDH2_M', 'PF16869': 'CNDH2_C', 'PF05786': 'Cnd2'}
ab = Alphabet.amino(); bg = Background(ab)
nm = lambda x: x if isinstance(x, str) else x.decode()

def profile(name):
    aln = WORK / f'{name}.aln'
    aln.write_text(subprocess.run([str(BIN / 'mafft'), '--auto', '--quiet', str(DER / f'condensin_{name}_seeds.fa')], check=True, capture_output=True, text=True).stdout)
    with MSAFile(str(aln), digital=True, alphabet=ab, format='afa') as f: msa = f.read()
    msa.name = name.encode()
    return Builder(ab).build_msa(msa, bg)[0]

def load(path):
    with SequenceFile(str(path), digital=True, alphabet=ab) as sf: return sf.read_block()

def main():
    H = {'caph2': profile('caph2'), 'caph': profile('caph')}
    pf = []
    for acc in PFAM:
        p = WORK / f'{acc}.hmm'
        if not p.exists(): p.write_bytes(gzip.decompress(urllib.request.urlopen(f'https://www.ebi.ac.uk/interpro/api/entry/pfam/{acc}?annotation=hmm').read()))
        with HMMFile(str(p)) as f: pf.append(f.read())
    rows = []
    for sp, fn in PROTEOMES.items():
        db = load(WORK / fn)
        byname = {nm(s.name): s for s in db}
        top = {n: {nm(h.name): h.evalue for h in hits} for n, hits in zip(H, pyhmmer.hmmsearch(list(H.values()), db, E=1e-3))}
        for acc in sorted(set(top['caph2']) - set(top['caph']), key=lambda a: top['caph2'][a]):
            s = byname[acc]
            doms = {}
            for scan in pyhmmer.hmmscan([s], pf, E=1e-3, domE=1e-3):
                for h in scan:
                    for d in h.domains: doms.setdefault(nm(h.name), f'{d.i_evalue:.0e} [{d.env_from}-{d.env_to}]')
            rows.append([sp, acc, str(len(s.sequence)), f"{top['caph2'][acc]:.1e}", f"{top['caph'].get(acc, 1):.1e}"] + [doms.get(PFAM[a], '') for a in PFAM])
        for acc in top['caph']:
            rows.append([sp, acc, str(len(byname[acc].sequence)), f"{top['caph2'].get(acc, 1):.1e}", f"{top['caph'][acc]:.1e}", '', '', '', 'CAP-H (control)'])
    text = [s for s in SequenceFile(str(WORK / 'h.faa'), digital=False)]
    null = []
    for rep in range(3):
        rnd = random.Random(rep); blk = []
        for s in text:
            q = list(s.sequence); rnd.shuffle(q); blk.append(TextSequence(name=s.name, sequence=''.join(q)).digitize(ab))
        hits = next(iter(pyhmmer.hmmsearch([H['caph2']], DigitalSequenceBlock(ab, blk), E=10)))
        null.append(min([h.evalue for h in hits], default=None))
    out = DER / 'condensin_caph2_profile.tsv'
    out.write_text('species\tprotein\tlength\tE_caph2_profile\tE_caph_profile\t' + '\t'.join(f'{PFAM[a]} ({a})' for a in PFAM) + '\n'
                   + ''.join('\t'.join(r) + '\n' for r in rows) + f'# null: shuffled H. exemplaris proteome, best CAP-H2 profile E per replicate: {null}\n')
    print(out.read_text())
    coexpression()

def coexpression():
    import sqlite3, collections, numpy as np
    db = sqlite3.connect(ROOT / 'data/derived/genome.sqlite')
    tg = collections.defaultdict(list)
    for t, g in db.execute("select target, gene from target_gene where dataset='yoshida2017'"): tg[t].append(g)
    samples = [s for (s,) in db.execute("select sample from sample where dataset='yoshida2017' order by sample")]
    rows = collections.defaultdict(lambda: collections.defaultdict(float))
    for t, s, v in db.execute("select target, sample, value from expression where dataset='yoshida2017' and unit='tpm'"):
        for g in tg.get(t, []): rows[g][s] += v
    genes = [g for g in rows if sum(rows[g][s] > 1 for s in samples) >= 6]
    R = np.array([[rows[g][s] for s in samples] for g in genes]).argsort(1).argsort(1).astype(float)
    R -= R.mean(1, keepdims=True); R /= np.linalg.norm(R, axis=1, keepdims=True)
    idx = {g: i for i, g in enumerate(genes)}
    G = {'CAP-H2 BV898_16921': 'gene-BV898_16921', 'CAP-H2 BV898_15594': 'gene-BV898_15594', 'CAP-D3': 'gene-BV898_04123',
         'CAP-G2': 'gene-BV898_12630', 'SMC2': 'gene-BV898_10436', 'SMC4': 'gene-BV898_17418', 'CAP-H': 'gene-BV898_02764'}
    lines = [f'# co-expression, Yoshida 2017, {len(samples)} samples, {len(genes)} expressed genes: Spearman (share of all genes ranking higher)']
    for x in list(G)[:2]:
        r = R @ R[idx[G[x]]]
        lines.append(f'# {x}: ' + '; '.join(f'{y} {r[idx[G[y]]]:.2f} ({100 * (r > r[idx[G[y]]]).mean():.1f}%)' for y in G if y != x))
    out = DER / 'condensin_caph2_profile.tsv'
    out.write_text(out.read_text() + '\n'.join(lines) + '\n'); print('\n'.join(lines))

if __name__ == '__main__':
    main()
