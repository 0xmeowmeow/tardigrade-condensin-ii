"""Microsynteny around the tardigrade CAP-H2 paralogs: are the 'a' copies in three species at orthologous loci, and the 'b' copies?

  ~/miniconda3/envs/tardigrade/bin/python tools/genome/caph2_synteny.py

Gene order from the GenBank/RefSeq GFFs (one representative protein per gene, the longest). Orthologs between each pair of
species are reciprocal best hits from DIAMOND blastp (--more-sensitive, E <= 1e-5) over whole proteomes. For every pair
of CAP-H2 loci (one per species), count anchors: genes within W genes of the first locus whose ortholog lies within W
genes of the second. Null: the same count for every orthologous gene pair between the two species (how much
neighbourhood an ortholog pair keeps in general), and for 20,000 random gene pairs (p_random_pairs: share with at least as
many anchors).
Output: data/derived/condensin_synteny.tsv and a summary on stdout.
"""
import collections, gzip, pathlib, random, re, subprocess

ROOT = pathlib.Path(__file__).resolve().parents[2]
BIN = pathlib.Path.home() / 'miniconda3/envs/tardigrade/bin'
WORK = ROOT / 'data/derived/condensin_synteny'; WORK.mkdir(parents=True, exist_ok=True)
SP = {'Hex': ('data/raw/GCA_002082055.1/GCA_002082055.1_nHd_3.1_genomic.gff.gz', 'data/raw/GCA_002082055.1/GCA_002082055.1_nHd_3.1_protein.faa.gz'),
      'Rva': ('data/raw/GCA_001949185.1/GCA_001949185.1_Rvar_4.0_genomic.gff.gz', 'data/raw/GCA_001949185.1/GCA_001949185.1_Rvar_4.0_protein.faa.gz'),
      'Pme': ('data/raw/GCF_019649055.1/GCF_019649055.1_Prichtersi_v1.0_genomic.gff.gz', 'data/raw/GCF_019649055.1/GCF_019649055.1_Prichtersi_v1.0_protein.faa.gz')}
LOCI = {'Hex a OWA52468.1': 'OWA52468.1', 'Hex b OWA51095.1': 'OWA51095.1', 'Rva a GAU96458.1': 'GAU96458.1', 'Rva b GAU88883.1': 'GAU88883.1',
        'Pme a XP_055355132.1': 'XP_055355132.1', 'Pme ? XP_055348578.1': 'XP_055348578.1'}
W = 20

def fasta(path):
    d, k = {}, None
    for l in gzip.open(ROOT / path, 'rt'):
        if l.startswith('>'): k = l[1:].split()[0]; d[k] = []
        else: d[k].append(l.strip())
    return {k: ''.join(v) for k, v in d.items()}

def genes(sp):
    gff, faa = SP[sp]; seqs = fasta(faa)
    loc = {}
    for l in gzip.open(ROOT / gff, 'rt'):
        if l.startswith('#'): continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'CDS': continue
        pid = re.search(r'protein_id=([^;]+)', f[8]); g = re.search(r'(?:locus_tag|gene)=([^;]+)', f[8])
        if not pid: continue
        pid = pid.group(1); gid = g.group(1).split('-')[0] if g else pid
        s, e = int(f[3]), int(f[4])
        if gid not in loc: loc[gid] = [f[0], s, e, set()]
        loc[gid][1] = min(loc[gid][1], s); loc[gid][2] = max(loc[gid][2], e); loc[gid][3].add(pid)
    rep, order, pos = {}, collections.defaultdict(list), {}
    for gid, (sc, s, e, pids) in loc.items():
        pids = [p for p in pids if p in seqs]
        if not pids: continue
        best = max(sorted(pids), key=lambda p: len(seqs[p]))
        for p in pids: rep[p] = best
        order[sc].append((s, best))
    for sc, lst in order.items():
        lst.sort()
        for i, (s, p) in enumerate(lst): pos[p] = (sc, i)
    reps = {p: seqs[p] for p in set(rep.values())}
    return rep, pos, reps

def rbh(a, b, A, B):
    out = {}
    for x, X, y, Y in ((a, A, b, B), (b, B, a, A)):
        fa = WORK / f'{x}.faa'; db = WORK / f'{y}.dmnd'
        if not fa.exists(): fa.write_text(''.join(f'>{k}\n{v}\n' for k, v in X[2].items()))
        fb = WORK / f'{y}.faa'
        if not fb.exists(): fb.write_text(''.join(f'>{k}\n{v}\n' for k, v in Y[2].items()))
        if not db.exists(): subprocess.run([str(BIN / 'diamond'), 'makedb', '--in', str(fb), '-d', str(db), '--quiet'], check=True)
        tsv = WORK / f'{x}_vs_{y}.tsv'
        if not tsv.exists():
            subprocess.run([str(BIN / 'diamond'), 'blastp', '-q', str(fa), '-d', str(db), '-o', str(tsv), '--more-sensitive', '-e', '1e-5',
                            '--max-target-seqs', '1', '--threads', '16', '--quiet'], check=True)
        best = {}
        for l in open(tsv):
            q, s = l.split('\t')[:2]
            best.setdefault(q, s)
        out[(x, y)] = best
    ab, ba = out[(a, b)], out[(b, a)]
    return {q: s for q, s in ab.items() if ba.get(s) == q}

def window(pos, byidx, p):
    sc, i = pos[p]
    return [byidx[sc][j] for j in range(max(0, i - W), i + W + 1) if j != i and j < len(byidx[sc])]

def main():
    S = {sp: genes(sp) for sp in SP}
    byidx = {}
    for sp, (rep, pos, reps) in S.items():
        d = collections.defaultdict(dict)
        for p, (sc, i) in pos.items(): d[sc][i] = p
        byidx[sp] = {sc: [d[sc][i] for i in range(len(d[sc]))] for sc in d}
    rows, summary = [], []
    for a, b in (('Hex', 'Rva'), ('Hex', 'Pme'), ('Rva', 'Pme')):
        R = rbh(a, b, S[a], S[b])
        posA, posB = S[a][1], S[b][1]
        def anchors(pa, pb):
            wb = set(window(posB, byidx[b], pb))
            return sum(1 for g in window(posA, byidx[a], pa) if R.get(g) in wb)
        # null: every ortholog pair
        null = [anchors(x, y) for x, y in R.items() if x in posA and y in posB]
        rnd = random.Random(0); keys = [k for k in R if k in posA]; vals = [R[k] for k in keys if R[k] in posB]
        rnull = [anchors(rnd.choice(keys), rnd.choice(vals)) for _ in range(20000)]
        summary.append(f'{a}-{b}: {len(R):,} reciprocal best hits; ortholog pairs keep a median of {sorted(null)[len(null) // 2]} anchors within ±{W} genes '
                       f'(mean {sum(null) / len(null):.2f}; {100 * sum(n >= 1 for n in null) / len(null):.0f}% keep at least one); random pairs mean {sum(rnull) / len(rnull):.3f}')
        for la, pa in LOCI.items():
            if not la.startswith(a): continue
            for lb, pb in LOCI.items():
                if not lb.startswith(b): continue
                ra, rb = S[a][0].get(pa, pa), S[b][0].get(pb, pb)
                n = anchors(ra, rb)
                pct = 100 * sum(x >= n for x in null) / len(null)
                k = sum(x >= n for x in rnull); prand = f'{k / len(rnull):.2g}' if k else f'<{1 / len(rnull):.0e}'
                rows.append([la, lb, str(n), f'{pct:.1f}', prand, posA[ra][0], posB[rb][0], 'yes' if R.get(ra) == rb else 'no'])
    out = ROOT / 'data/derived/condensin_synteny.tsv'
    out.write_text('locus_1\tlocus_2\tanchors_within_20_genes\tpct_of_ortholog_pairs_with_at_least_as_many\tp_random_pairs\tscaffold_1\tscaffold_2\treciprocal_best_hit\n'
                   + ''.join('\t'.join(r) + '\n' for r in rows) + ''.join(f'# {s}\n' for s in summary))
    print(out.read_text())

if __name__ == '__main__':
    main()
