"""Data for the condensin II preprint figures: Pfam domain positions and expression profiles.

  ~/miniconda3/envs/tardigrade/bin/python tools/genome/condensin_figdata.py   (after caph2_profile.py)

Writes data/derived/condensin_figdata.json:
  domains:    Pfam CNDH2_N / CNDH2_M / CNDH2_C / Cnd2 envelopes (i-E <= 1e-3) for human NCAPH2 and NCAPH, the six
              tardigrade CAP-H2 candidates (one Paramacrobiotus isoform) and the three tardigrade CAP-H controls
  yoshida:    mean TPM per condition (Yoshida 2017, 16 conditions x 3 reps) for the H. exemplaris condensin genes
  levin:      counts per million per single embryo against minutes (Levin 2016, 62 embryos)
Plot with tools/genome/condensin_figures.py (system python, matplotlib).
"""
import collections, json, pathlib, sqlite3
import pyhmmer
from pyhmmer.easel import Alphabet, SequenceFile
from pyhmmer.plan7 import HMMFile

ROOT = pathlib.Path(__file__).resolve().parents[2]
WORK = ROOT / 'data/derived/condensin_check'
DER = ROOT / 'data/derived'
PFAM = {'PF06278': 'CNDH2_N', 'PF16858': 'CNDH2_M', 'PF16869': 'CNDH2_C', 'PF05786': 'Cnd2'}
SEQS = [  # label, file, accession
    ('Human NCAPH2 (CAP-H2)', 'condensin_caph2_seeds.fa', 'Q6IBW4'),
    ('H. exemplaris OWA52468.1', 'h.faa', 'OWA52468.1'),
    ('H. exemplaris OWA51095.1', 'h.faa', 'OWA51095.1'),
    ('R. varieornatus GAU96458.1', 'r.faa', 'GAU96458.1'),
    ('R. varieornatus GAU88883.1', 'r.faa', 'GAU88883.1'),
    ('P. metropolitanus XP_055348578.1', 'paramacrobiotus.faa', 'XP_055348578.1'),
    ('P. metropolitanus XP_055355132.1', 'paramacrobiotus.faa', 'XP_055355132.1'),
    ('Human NCAPH (CAP-H)', 'condensin_caph_seeds.fa', 'Q15003'),
    ('H. exemplaris OQV23316.1 (CAP-H)', 'h.faa', 'OQV23316.1'),
    ('R. varieornatus GAU89914.1 (CAP-H)', 'r.faa', 'GAU89914.1'),
    ('P. metropolitanus XP_055347164.1 (CAP-H)', 'paramacrobiotus.faa', 'XP_055347164.1'),
]
GENES = {'CAP-H2 BV898_16921': 'gene-BV898_16921', 'CAP-H2 BV898_15594': 'gene-BV898_15594', 'CAP-D3': 'gene-BV898_04123',
         'CAP-G2': 'gene-BV898_12630', 'SMC2': 'gene-BV898_10436', 'SMC4': 'gene-BV898_17418', 'CAP-H': 'gene-BV898_02764'}
nm = lambda x: x if isinstance(x, str) else x.decode()

def domains():
    ab = Alphabet.amino()
    pf = []
    for acc in PFAM:
        with HMMFile(str(WORK / f'{acc}.hmm')) as f: pf.append(f.read())
    out = []
    for label, fn, acc in SEQS:
        path = (DER if fn.startswith('condensin_') else WORK) / fn
        with SequenceFile(str(path), digital=True, alphabet=ab) as sf:
            s = next(x for x in sf.read_block() if acc in nm(x.name).split('|') or nm(x.name) == acc)
        doms = []
        for scan in pyhmmer.hmmscan([s], pf, E=1e-3, domE=1e-3):
            for h in scan:
                for d in h.domains:
                    if d.i_evalue <= 1e-3: doms.append({'pfam': nm(h.name), 'from': d.env_from, 'to': d.env_to, 'E': d.i_evalue})
        out.append({'label': label, 'acc': acc, 'length': len(s.sequence), 'domains': doms})
    return out

def expression():
    db = sqlite3.connect(ROOT / 'data/derived/genome.sqlite')
    ids = tuple(GENES.values())
    y = collections.defaultdict(lambda: collections.defaultdict(list))
    order = []
    for g, samp, stage, cond, v in db.execute(f"""select t.gene, e.sample, s.stage, s.condition, sum(e.value) from expression e
            join target_gene t on t.dataset=e.dataset and t.target=e.target join sample s on s.dataset=e.dataset and s.sample=e.sample
            where e.dataset='yoshida2017' and e.unit='tpm' and t.gene in {ids} group by t.gene, e.sample order by e.sample"""):
        key = f'{stage}: {cond}'
        if key not in order: order.append(key)
        y[g][key].append(v)
    yosh = {'conditions': order, 'mean_tpm': {n: [sum(y[g][c]) / len(y[g][c]) for c in order] for n, g in GENES.items() if g in y}}
    tot = dict(db.execute("select sample, sum(value) from expression where dataset='levin2016' group by sample"))
    tmin = dict(db.execute("select sample, time_min from sample where dataset='levin2016'"))
    lv = collections.defaultdict(lambda: collections.defaultdict(float))
    for g, samp, v in db.execute(f"""select t.gene, e.sample, e.value from expression e join target_gene t
            on t.dataset=e.dataset and t.target=e.target where e.dataset='levin2016' and t.gene in {ids}"""):
        lv[g][samp] += v
    samples = sorted(tmin, key=lambda s: tmin[s])
    levin = {'minutes': [tmin[s] for s in samples],
             'cpm': {n: [1e6 * lv[g][s] / tot[s] for s in samples] for n, g in GENES.items() if g in lv}}
    return yosh, levin

if __name__ == '__main__':
    yosh, levin = expression()
    out = DER / 'condensin_figdata.json'
    out.write_text(json.dumps({'domains': domains(), 'yoshida': yosh, 'levin': levin}, indent=1))
    print('wrote', out)
