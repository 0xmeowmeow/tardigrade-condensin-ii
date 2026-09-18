"""Check Hoencamp et al. 2021's statement that the tardigrade lacks every condensin accessory subunit.

  ~/miniconda3/envs/tardigrade/bin/python tools/genome/condensin_check.py

Hoencamp et al. 2021 (Science 372:984, Supplementary Materials p. 10–11) found SMC2 and SMC4 in "Hypsibius
dujardini" but none of the condensin I (CAP-D2, CAP-G, CAP-H) or condensin II (CAP-D3, CAP-G2, CAP-H2) accessory
subunits, searching with tblastn against nucleotide databases, and left it to be investigated further.

Here: human reviewed UniProt sequences of all eight condensin subunits, blastp against three tardigrade
proteomes (H. exemplaris nHd_3.1, R. varieornatus Rvar_4.0, Paramacrobiotus GCF_019649055.1), best hit kept only
if it is a reciprocal best hit against the whole human reviewed proteome. The kleisins (CAP-H, CAP-H2) evolve
fast, so they also get 5 rounds of PSI-BLAST, with the same reciprocal check. Output: data/derived/condensin_check.tsv.
"""
import gzip, pathlib, subprocess, urllib.request, collections

ROOT = pathlib.Path(__file__).resolve().parents[2]
BIN = pathlib.Path.home() / 'miniconda3/envs/tardigrade/bin'
WORK = ROOT / 'data/derived/condensin_check'; WORK.mkdir(parents=True, exist_ok=True)
UNI = ROOT / 'data/raw/uniprot'; UNI.mkdir(parents=True, exist_ok=True)
SUBUNITS = {'O95347': 'SMC2', 'Q9NTJ3': 'SMC4', 'Q15021': 'NCAPD2', 'Q9BPX3': 'NCAPG', 'Q15003': 'NCAPH',
            'P42695': 'NCAPD3', 'Q86XI2': 'NCAPG2', 'Q6IBW4': 'NCAPH2'}
PROTEOMES = {'H. exemplaris': 'data/raw/GCA_002082055.1/GCA_002082055.1_nHd_3.1_protein.faa.gz',
             'R. varieornatus': 'data/raw/GCA_001949185.1/GCA_001949185.1_Rvar_4.0_protein.faa.gz',
             'Paramacrobiotus (GCF_019649055.1)': 'data/raw/GCF_019649055.1/GCF_019649055.1_Prichtersi_v1.0_protein.faa.gz'}

def run(*a): return subprocess.run([str(BIN / a[0]), *a[1:]], check=True, capture_output=True, text=True).stdout

def fasta(path):
    seqs, k = collections.OrderedDict(), None
    op = gzip.open if str(path).endswith('.gz') else open
    for l in op(path, 'rt'):
        if l.startswith('>'): k = l[1:].split()[0]; seqs[k] = [l.strip(), '']
        elif k: seqs[k][1] += l.strip()
    return seqs

def main():
    human = UNI / 'human_sprot.fasta'
    if not human.exists():
        urllib.request.urlretrieve('https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28reviewed%3Atrue%29+AND+%28organism_id%3A9606%29', human)
    hs = fasta(human)
    q = WORK / 'condensin_human.fa'
    q.write_text(''.join(f'>{SUBUNITS[k.split("|")[1]]}\n{v[1]}\n' for k, v in hs.items() if k.split('|')[1] in SUBUNITS))
    run('makeblastdb', '-in', str(human), '-dbtype', 'prot', '-out', str(WORK / 'human'))
    hname = {k: v[0].split(' OS=')[0].split(' ', 1)[1] for k, v in hs.items()}
    rows = []
    for sp, path in PROTEOMES.items():
        tag = sp.split()[0].strip('.').lower()
        seqs = fasta(ROOT / path); faa = WORK / f'{tag}.faa'
        faa.write_text(''.join(f'{v[0]}\n{v[1]}\n' for v in seqs.values()))
        run('makeblastdb', '-in', str(faa), '-dbtype', 'prot', '-out', str(WORK / tag))
        hits = {}
        for l in run('blastp', '-query', str(q), '-db', str(WORK / tag), '-evalue', '1e-5', '-max_target_seqs', '5', '-num_threads', '8', '-outfmt', '6 qseqid sseqid evalue bitscore length').splitlines():
            g, s, e, b, ln = l.split('\t')
            if g not in hits or float(b) > hits[g][2]: hits[g] = (s, float(e), float(b), int(ln), 'blastp')
        for g in ('NCAPH', 'NCAPH2'):
            if g in hits: continue
            one = WORK / f'{g}.fa'; one.write_text(f'>{g}\n' + [v[1] for k, v in hs.items() if SUBUNITS.get(k.split("|")[1]) == g][0] + '\n')
            out = run('psiblast', '-query', str(one), '-db', str(WORK / tag), '-num_iterations', '5', '-evalue', '1e-3', '-max_target_seqs', '5', '-num_threads', '8', '-outfmt', '6 sseqid evalue bitscore length')
            best = sorted((l.split('\t') for l in out.splitlines() if l and not l.startswith('Search') and '\t' in l), key=lambda r: float(r[1]))
            if best: hits[g] = (best[0][0], float(best[0][1]), float(best[0][2]), int(best[0][3]), 'psiblast')
        sub = WORK / f'{tag}.hits.faa'
        sub.write_text(''.join(f'{seqs[h[0]][0]}\n{seqs[h[0]][1]}\n' for h in hits.values()))
        back = {}
        for l in run('blastp', '-query', str(sub), '-db', str(WORK / 'human'), '-evalue', '1e-3', '-max_target_seqs', '5', '-num_threads', '8', '-outfmt', '6 qseqid sseqid bitscore').splitlines():
            s, h, b = l.split('\t')
            if s not in back or float(b) > back[s][1]: back[s] = (h, float(b))
        for acc, g in SUBUNITS.items():
            if g not in hits: rows.append([sp, g, '', '', '', '', '', 'not found']); continue
            s, e, b, ln, how = hits[g]; bh = back.get(s, ('', 0))[0]
            rbh = bh.split('|')[1] == acc if bh else False
            rows.append([sp, g, s, seqs[s][0].split(' ', 1)[1], f'{e:.1e}', how, hname.get(bh, 'no human hit'), 'reciprocal best hit' if rbh else 'not reciprocal'])
    out = ROOT / 'data/derived/condensin_check.tsv'
    out.write_text('species\tsubunit\ttardigrade_protein\tannotation\tevalue\tmethod\tbest_human_hit_back\tverdict\n' + ''.join('\t'.join(r) + '\n' for r in rows))
    for r in rows: print(f'{r[0]:36} {r[1]:7} {r[2]:16} {r[4]:8} {r[5]:8} {r[7]:20} back: {r[6]}')

if __name__ == '__main__':
    main()
