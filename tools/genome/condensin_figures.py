"""Figures for the condensin II preprint (docs/papers/condensin-ii/figs/).

  python3 tools/genome/condensin_figures.py        (system python with matplotlib; run condensin_figdata.py first)

fig1-domains:     Pfam domain maps of human CAP-H2 and CAP-H and the tardigrade kleisins
fig2-tree:        maximum-likelihood kleisin tree (data/derived/condensin_tree/kleisins.treefile), if present
fig3-expression:  (a) mean TPM per Yoshida 2017 condition, (b) Levin 2016 single embryos, (c) co-expression ranks
"""
import collections, json, pathlib, re, sqlite3
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / 'docs/papers/condensin-ii/figs'
D = json.loads((ROOT / 'data/derived/condensin_figdata.json').read_text())
COL = {'CNDH2_N': '#2b6cb0', 'CNDH2_M': '#63b3ed', 'CNDH2_C': '#c05621', 'Cnd2': '#6b8e23'}
plt.rcParams.update({'font.size': 8, 'font.family': 'DejaVu Sans', 'axes.spines.top': False, 'axes.spines.right': False})

def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ('pdf', 'png'): fig.savefig(OUT / f'{name}.{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)

def fig_domains():
    rows = D['domains']
    fig, ax = plt.subplots(figsize=(6.8, 3.6))
    for i, r in enumerate(rows):
        y = len(rows) - i - (0.6 if i >= 7 else 0)
        ax.plot([0, r['length']], [y, y], color='#555', lw=1.2, solid_capstyle='butt')
        for d in r['domains']:
            ax.add_patch(plt.Rectangle((d['from'], y - 0.3), d['to'] - d['from'], 0.6, color=COL[d['pfam']], lw=0))
        ax.text(-20, y, r['label'], ha='right', va='center', fontsize=7, style='italic' if not r['label'].startswith('Human') else 'normal')
        ax.text(r['length'] + 15, y, f"{r['length']} aa", va='center', fontsize=6, color='#555')
    ax.set_xlim(-10, 1200); ax.set_ylim(0, len(rows) + 1); ax.set_yticks([]); ax.spines['left'].set_visible(False)
    ax.set_xlabel('residue')
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=COL[k]) for k in COL], labels=list(COL), ncol=4, frameon=False,
              loc='upper center', bbox_to_anchor=(0.45, 1.1), fontsize=7)
    save(fig, 'fig1-domains')

def fig_tree():
    tf = ROOT / 'data/derived/condensin_tree/kleisins.treefile'
    if not tf.exists(): print('no tree yet'); return
    from Bio import Phylo
    t = Phylo.read(str(tf), 'newick')
    t.root_with_outgroup(*[c for c in t.get_terminals() if c.name.startswith('OUT__')])
    t.ladderize()
    colour = {'CAPH2': '#2b6cb0', 'CAPH': '#6b8e23', 'TARD': '#c05621', 'OUT': '#777777'}
    labels = {}
    for c in t.get_terminals():
        kind, sp, acc = c.name.split('__')
        labels[c] = f"{sp.replace('_', ' ')} {acc.replace('_CAPH2a', ' (a)').replace('_CAPH2b', ' (b)').replace('_CAPH2', '').replace('_CAPH', '')}"
    for c in t.get_nonterminals():
        c.name = None
        if c.confidence is not None and c.confidence < 70: c.confidence = None
    fig, ax = plt.subplots(figsize=(6.8, 9.5))
    Phylo.draw(t, axes=ax, do_show=False, label_func=lambda c: labels.get(c, ''),
               label_colors={labels[c]: colour[c.name.split('__')[0]] for c in labels},
               branch_labels=lambda c: f'{c.confidence:.0f}' if c.confidence is not None and not c.is_terminal() else '')
    for txt in ax.texts: txt.set_fontsize(5.5)
    ax.set_ylabel(''); ax.set_yticks([]); ax.spines['left'].set_visible(False)
    ax.set_xlabel('substitutions per site')
    save(fig, 'fig2-tree')

def coexpression_all():
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
    return genes, R

def fig_expression():
    G = {'CAP-H2 BV898_16921': '#c05621', 'CAP-H2 BV898_15594': '#dd6b20', 'CAP-D3': '#2b6cb0', 'CAP-G2': '#4299e1',
         'SMC2': '#555', 'SMC4': '#999', 'CAP-H': '#6b8e23'}
    fig = plt.figure(figsize=(6.8, 5.6))
    a = fig.add_subplot(2, 2, (1, 2))
    Y = D['yoshida']; x = np.arange(len(Y['conditions']))
    for n, c in G.items():
        v = np.array(Y['mean_tpm'][n]); a.plot(x, v / v.max(), '-o', ms=2.5, lw=1.2 if 'H2' in n else 0.8, color=c, label=n)
    short = [re.sub(r'day (\d) after laying', r'E\1', re.sub(r'juvenile day (\d)', r'J\1', c.split(': ')[1])) for c in Y['conditions']]
    short = [s.replace('active', 'act').replace('anhydrobiotic tun', 'tun') for s in short]
    a.set_xticks(x); a.set_xticklabels(short, rotation=45, ha='right', fontsize=6.5)
    a.set_ylabel('TPM / gene maximum'); a.legend(fontsize=6, ncol=4, frameon=False, loc='upper center', bbox_to_anchor=(0.5, 1.16))
    a.set_title('a', loc='left', fontweight='bold', x=-0.07)
    b = fig.add_subplot(2, 2, 3)
    L = D['levin']; h = np.array(L['minutes']) * 2 / 60  # the workbook's minutes are half of real time (matched to Yoshida 2017 daily samples)
    for n, c in G.items():
        if n not in L['cpm']: continue
        v = np.array(L['cpm'][n]); k = np.ones(5) / 5
        b.plot(h, v, '.', ms=2, color=c, alpha=0.4); b.plot(h[2:-2], np.convolve(v, k, 'valid'), lw=1, color=c)
    b.set_xlabel('hours after laying (single embryos)'); b.set_ylabel('counts per million'); b.set_yscale('symlog', linthresh=10)
    b.set_title('b', loc='left', fontweight='bold', x=-0.15)
    c = fig.add_subplot(2, 2, 4)
    genes, R = coexpression_all()
    idx = {g: i for i, g in enumerate(genes)}
    r = R @ R[idx['gene-BV898_16921']]
    c.hist(np.delete(r, idx['gene-BV898_16921']), bins=80, color='#ccc')
    ids = {'CAP-D3': 'gene-BV898_04123', 'CAP-G2': 'gene-BV898_12630', 'SMC2': 'gene-BV898_10436', 'SMC4': 'gene-BV898_17418',
           'CAP-H': 'gene-BV898_02764', 'CAP-H2 BV898_15594': 'gene-BV898_15594'}
    for i, (n, g) in enumerate(ids.items()):
        c.axvline(r[idx[g]], color=G[n], lw=1); c.text(r[idx[g]], c.get_ylim()[1] * (0.92 - 0.1 * i), n.replace('CAP-H2 ', ''), fontsize=5.5, color=G[n], ha='right')
    c.set_xlabel('Spearman ρ with CAP-H2 BV898_16921'); c.set_ylabel(f'genes (n = {len(genes) - 1:,})')
    c.set_title('c', loc='left', fontweight='bold', x=-0.15)
    fig.tight_layout()
    save(fig, 'fig3-expression')

if __name__ == '__main__':
    fig_domains(); fig_expression(); fig_tree()
    print('figures in', OUT)
