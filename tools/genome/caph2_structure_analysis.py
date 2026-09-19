"""Summarise the Boltz-2 predictions of the tardigrade CAP-H2 kleisins (inputs from caph2_structure_inputs.py).

  ~/miniconda3/envs/boltz/bin/python tools/genome/caph2_structure_analysis.py [--n2]
  --n2: the N-terminal retry with the longer SMC2 neck (runs_n2, inputs from --n-long) -> summary_n2.tsv

Per prediction, the best of the diffusion samples by Boltz confidence score, and the spread over samples:
  iptm            interface pTM (0-1): how confident the model is in the relative placement of the two chains
  plddt_A         mean pLDDT of chain A (the kleisin region, or the full kleisin for mono_*)
  ipae            local interface error (A): for each well-folded kleisin residue (pLDDT >= 70) in contact, mean PAE to the SMC residues
                  within 8 A of it, averaged; core_contacts counts those residues
  contacts_A      kleisin residues with any heavy atom within 5 A of the SMC chain
For mono_*: mean pLDDT over the Pfam CNDH2_N and CNDH2_C envelopes, against the AlphaFold DB model where one exists.
Output: data/derived/caph2_structure/summary.tsv
"""
import glob, json, pathlib, sys
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
D = ROOT / 'data/derived/caph2_structure'
N2 = '--n2' in sys.argv
PRED = D / ('runs_n2/boltz_results_inputs_n2/predictions' if N2 else 'runs/boltz_results_inputs/predictions')
AF = {'hsa': 'Q6IBW4', 'hex_a': 'A0A9X6NG78', 'hex_b': 'A0A9X6RKL3', 'rva_a': 'A0A1D1VD77', 'rva_b': 'A0A1D1UR96'}

def atoms(pdb):
    out = []
    for l in open(pdb):
        if l.startswith(('ATOM', 'HETATM')) and l[76:78].strip() != 'H':
            out.append((l[21], int(l[22:26]), l[12:16].strip(), np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]), float(l[60:66])))
    return out

def contacts(at, cut=5.0):
    A = [(r, x) for c, r, n, x, b in at if c == 'A']; B = np.array([x for c, r, n, x, b in at if c == 'B'])
    if not len(B): return set(), set()
    ra, rb = set(), set()
    for r, x in A:
        d = np.linalg.norm(B - x, axis=1)
        if d.min() < cut: ra.add(r)
    Aa = np.array([x for r, x in A]); Br = [r for c, r, n, x, b in at if c == 'B']
    for i, x in enumerate(B):
        if np.linalg.norm(Aa - x, axis=1).min() < cut: rb.add(Br[i])
    return ra, rb

def plddt_by_res(at, chain):
    d = {}
    for c, r, n, x, b in at:
        if c == chain: d.setdefault(r, []).append(b)
    return {r: np.mean(v) for r, v in d.items()}

def af_domain_plddt(acc, doms):
    f = ROOT / f'data/raw/alphafold/AF-{acc}-F1-model_v6.pdb'
    if not f.exists(): return None
    p = {}
    for l in open(f):
        if l.startswith('ATOM') and l[12:16].strip() == 'CA': p[int(l[22:26])] = float(l[60:66])
    return [np.mean([p[r] for r in range(a, b + 1) if r in p]) for a, b in doms]

def main():
    cons = json.loads((D / 'constructs.json').read_text())
    fd = json.loads((ROOT / 'data/derived/condensin_figdata.json').read_text())
    doms = {d['acc']: {x['pfam']: (x['from'], x['to']) for x in d['domains']} for d in fd['domains']}
    rows = []
    for pdir in sorted(PRED.iterdir()):
        name = pdir.name
        confs = sorted(glob.glob(str(pdir / 'confidence_*_model_*.json')))
        if not confs: continue
        cs = [json.loads(open(c).read()) for c in confs]
        best = max(range(len(cs)), key=lambda i: cs[i]['confidence_score'])
        at = atoms(pdir / f'{name}_model_{best}.pdb')
        pl = plddt_by_res(at, 'A')
        plA = np.mean(list(pl.values()))
        plA = plA * 100 if plA <= 1 else plA
        row = {'name': name, 'iptm': cs[best].get('iptm', float('nan')), 'iptm_range': f"{min(c.get('iptm', 0) for c in cs):.2f}-{max(c.get('iptm', 0) for c in cs):.2f}",
               'plddt_A': plA, 'contacts_A': '', 'ipae': '', 'core_contacts': '', 'domains': ''}
        if not name.startswith('mono_'):
            ra, rb = contacts(at)
            row['contacts_A'] = len(ra)
            pae = np.load(pdir / f'pae_{name}_model_{best}.npz'); pae = pae[pae.files[0]]
            nA = len(pl)
            # local interface error: for each well-folded kleisin residue (pLDDT >= 70) touching the SMC chain, the mean PAE to the
            # SMC residues within 8 A of it; the flexible flanks of the construct are left out
            resA = {}; resB = {}
            for c, r, n, x, b in at: (resA if c == 'A' else resB).setdefault(r, []).append(x)
            Bk = sorted(resB); Bx = [np.array(resB[r]) for r in Bk]
            vals = []
            for r in sorted(ra):
                if (pl[r] * 100 if pl[r] <= 1 else pl[r]) < 70: continue
                xa = np.array(resA[r]); near = [nA + q - 1 for q, xb in zip(Bk, Bx) if np.min(np.linalg.norm(xa[:, None] - xb[None], axis=2)) < 8]
                if near: vals.append(pae[r - 1][near].mean())
            row['ipae'] = float(np.mean(vals)) if vals else ''
            row['core_contacts'] = len(vals)
        else:
            key = name[5:]; acc = cons[key]['acc']; dd = doms.get(acc, {})
            want = [dd[k] for k in ('CNDH2_N', 'CNDH2_C') if k in dd]
            ours = [np.mean([pl[r] for r in range(a, b + 1) if r in pl]) for a, b in want]
            ours = [v * 100 if v <= 1 else v for v in ours]
            afv = af_domain_plddt(AF[key], want) if key in AF else None
            row['domains'] = 'N {:.0f}, C {:.0f}'.format(*ours) + (' (AlphaFold DB: N {:.0f}, C {:.0f})'.format(*afv) if afv else '')
        rows.append(row)
    out = D / ('summary_n2.tsv' if N2 else 'summary.tsv')
    cols = ['name', 'iptm', 'iptm_range', 'plddt_A', 'ipae', 'contacts_A', 'core_contacts', 'domains']
    fmt = lambda v: f'{v:.3f}' if isinstance(v, float) else str(v)
    out.write_text('\t'.join(cols) + '\n' + ''.join('\t'.join(fmt(r[c]) for c in cols) + '\n' for r in rows))
    print(out.read_text())

if __name__ == '__main__':
    main()
