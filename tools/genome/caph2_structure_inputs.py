"""Inputs for structure predictions of the tardigrade CAP-H2 kleisins (Boltz-2), alone and bound to the condensin ring.

  ~/miniconda3/envs/tardigrade/bin/python tools/genome/caph2_structure_inputs.py

Kleisins bridge the two SMC ATPase heads: in condensin the kleisin N-terminal domain binds the SMC2 neck (the coiled coil
leaving the head) and the C-terminal winged-helix domain binds the SMC4 head. Predictions per species (human as the
reference; H. exemplaris a and b; R. varieornatus a and b; Paramacrobiotus a and b):
  mono_<sp>      full-length CAP-H2 alone (compare domain confidence with AlphaFold DB)
  n_<sp>         CAP-H2 N-terminal region + that species' SMC2 head-and-neck construct
  c_<sp>         CAP-H2 C-terminal region + that species' SMC4 head-and-neck construct
Controls: swapped partners (C region with SMC2, N region with SMC4), the same regions with residues shuffled, and the
condensin I kleisin CAP-H C region with SMC4 (a real kleisin, the other complex).
SMC head constructs: the N-terminal lobe up to ~50 residues into the coiled coil, a GGSGGSGGSGGS linker, then from ~50
residues before the C-terminal lobe to the end. Boundaries come from human AlphaFold models (long helical runs) and are
mapped to each tardigrade protein through a MAFFT alignment. Kleisin regions: the Pfam CNDH2_N / CNDH2_C envelopes
(condensin_figdata.json) extended by 40 residues each side.
MSAs: each chain's MSA from the ColabFold MMseqs2 server (Boltz --use_msa_server), all inputs public sequences.
Output: data/derived/caph2_structure/inputs/*.yaml
"""
import json, pathlib, random, subprocess, gzip

ROOT = pathlib.Path(__file__).resolve().parents[2]
BIN = pathlib.Path.home() / 'miniconda3/envs/tardigrade/bin'
OUT = ROOT / 'data/derived/caph2_structure'; INP = OUT / 'inputs'; INP.mkdir(parents=True, exist_ok=True)
WORK = ROOT / 'data/derived/condensin_check'
LINK = 'GGSGGSGGSGGS'
HUMAN_SMC = {'SMC2': ('O95347', (1, 195), (1000, 1197)), 'SMC4': ('Q9NTJ3', (56, 313), (1105, 1288))}
SMC = {'hsa': {'SMC2': 'O95347', 'SMC4': 'Q9NTJ3'},
       'hex': {'SMC2': 'OQV15428.1', 'SMC4': 'OWA52980.1'},
       'rva': {'SMC2': 'GAV04166.1', 'SMC4': 'GAV06634.1'},
       'pme': {'SMC2': 'XP_055343015.1', 'SMC4': 'XP_055342529.1'}}
KLEISIN = {'hsa': ('Q6IBW4', 'hsa'), 'hex_a': ('OWA52468.1', 'hex'), 'hex_b': ('OWA51095.1', 'hex'),
           'rva_a': ('GAU96458.1', 'rva'), 'rva_b': ('GAU88883.1', 'rva'),
           'pme_a': ('XP_055355132.1', 'pme'), 'pme_b': ('XP_055348578.1', 'pme')}
CAPH_HEX = 'OQV23316.1'

def fasta(p):
    d, k = {}, None
    op = gzip.open if str(p).endswith('.gz') else open
    for l in op(p, 'rt'):
        l = l.strip()
        if l.startswith('>'): k = l[1:].split()[0]; d[k] = ''
        else: d[k] += l
    return d

def uniprot(acc):
    import urllib.request
    f = OUT / f'{acc}.fasta'
    if not f.exists(): f.write_bytes(urllib.request.urlopen(f'https://rest.uniprot.org/uniprotkb/{acc}.fasta').read())
    return ''.join(f.read_text().split('\n')[1:])

def seqs():
    S = {}
    for f in ('h.faa', 'r.faa', 'paramacrobiotus.faa'): S.update(fasta(WORK / f))
    for acc in ('O95347', 'Q9NTJ3', 'Q6IBW4'): S[acc] = uniprot(acc)
    return S

def map_range(hseq, tseq, a, b):
    """map human residue range [a,b] (1-based) onto tseq via a pairwise mafft alignment"""
    fa = OUT / 'pair.fa'; fa.write_text(f'>h\n{hseq}\n>t\n{tseq}\n')
    aln = subprocess.run([str(BIN / 'mafft'), '--auto', '--quiet', str(fa)], capture_output=True, text=True, check=True).stdout
    al = {}; k = None
    for l in aln.splitlines():
        if l.startswith('>'): k = l[1:]; al[k] = ''
        else: al[k] += l.strip()
    hi = ti = 0; m = {}
    for x, y in zip(al['h'], al['t']):
        if x != '-': hi += 1
        if y != '-': ti += 1
        if x != '-' and y != '-': m[hi] = ti
    near = lambda r, step: next(m[q] for q in range(r, r + 200 * step, step) if q in m)
    return near(a, 1), near(b, -1)

def construct(S, sp, smc):
    acc, (n0, n1), (c0, c1) = HUMAN_SMC[smc]
    h = S[acc]; t = S[SMC[sp][smc]]
    if sp == 'hsa': a, b, c, d = n0, n1, c0, c1
    else:
        a, b = map_range(h, t, n0, n1); c, d = map_range(h, t, c0, c1)
    return t[a - 1:b] + LINK + t[c - 1:d], f'{a}-{b}+{c}-{d}'

def region(seq, dom, pad=40):
    return seq[max(0, dom['from'] - 1 - pad):dom['to'] + pad], f"{max(1, dom['from'] - pad)}-{min(len(seq), dom['to'] + pad)}"

def yaml(name, chains):
    lines = ['version: 1', 'sequences:']
    for cid, s in chains:
        lines += ['  - protein:', f'      id: {cid}', f'      sequence: {s}']
    (INP / f'{name}.yaml').write_text('\n'.join(lines) + '\n')

def main():
    S = seqs()
    fd = json.loads((ROOT / 'data/derived/condensin_figdata.json').read_text())
    doms = {d['acc']: {x['pfam']: x for x in d['domains']} for d in fd['domains']}
    meta = {}
    for sp in SMC:
        for smc in ('SMC2', 'SMC4'):
            meta[f'{sp}_{smc}'] = construct(S, sp, smc)
    rnd = random.Random(0)
    shuf = lambda s: ''.join(rnd.sample(s, len(s)))
    for key, (acc, sp) in KLEISIN.items():
        k = S[acc]; d = doms[acc]
        N, nr = region(k, d['CNDH2_N']); C, cr = region(k, d['CNDH2_C'])
        yaml(f'mono_{key}', [('A', k)])
        yaml(f'n_{key}', [('A', N), ('B', meta[f'{sp}_SMC2'][0])])
        yaml(f'c_{key}', [('A', C), ('B', meta[f'{sp}_SMC4'][0])])
        if key in ('hsa', 'hex_a', 'hex_b'):
            yaml(f'ctl_swapN_{key}', [('A', N), ('B', meta[f'{sp}_SMC4'][0])])
            yaml(f'ctl_swapC_{key}', [('A', C), ('B', meta[f'{sp}_SMC2'][0])])
            yaml(f'ctl_shufN_{key}', [('A', shuf(N)), ('B', meta[f'{sp}_SMC2'][0])])
            yaml(f'ctl_shufC_{key}', [('A', shuf(C)), ('B', meta[f'{sp}_SMC4'][0])])
        meta[key] = {'acc': acc, 'N': nr, 'C': cr, 'len': len(k)}
    capH = S[CAPH_HEX]
    # CAP-H C-terminal region: the last 200 residues, 837-1036 (the Cnd2 hit near the C terminus)
    yaml('ctl_caph_c_hex', [('A', capH[-200:]), ('B', meta['hex_SMC4'][0])])
    (OUT / 'constructs.json').write_text(json.dumps(meta, indent=1))
    print(len(list(INP.glob('*.yaml'))), 'inputs;', {k: v[1] if isinstance(v, tuple) else v for k, v in meta.items()})

def n_long():
    """Round 2 for the N-terminal region: SMC2 constructs with ~150 residues of each coiled-coil strand (human 1-300 + 950-1197)
    in case the neck site was cut short. Output: data/derived/caph2_structure/inputs_n2/"""
    global INP
    INP = OUT / 'inputs_n2'; INP.mkdir(parents=True, exist_ok=True)
    HUMAN_SMC['SMC2'] = ('O95347', (1, 300), (950, 1197))
    S = seqs(); fd = json.loads((ROOT / 'data/derived/condensin_figdata.json').read_text())
    doms = {d['acc']: {x['pfam']: x for x in d['domains']} for d in fd['domains']}
    rnd = random.Random(1); shuf = lambda s: ''.join(rnd.sample(s, len(s))); meta = {}
    for key in ('hsa', 'hex_a', 'hex_b', 'rva_b', 'pme_b'):
        acc, sp = KLEISIN[key]; N, nr = region(S[acc], doms[acc]['CNDH2_N'])
        smc2, rng = construct(S, sp, 'SMC2'); meta[key] = {'N': nr, 'SMC2': rng}
        yaml(f'n2_{key}', [('A', N), ('B', smc2)])
        if key in ('hsa', 'hex_b'): yaml(f'ctl_n2shuf_{key}', [('A', shuf(N)), ('B', smc2)])
    (OUT / 'constructs_n2.json').write_text(json.dumps(meta, indent=1)); print(meta)

if __name__ == '__main__':
    import sys
    n_long() if '--n-long' in sys.argv else main()
