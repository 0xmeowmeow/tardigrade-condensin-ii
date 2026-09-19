"""Turn a predicted complex (PDB) into a U3D mesh for embedding in the preprint PDF (media9), plus the view settings.

  ~/miniconda3/envs/trace/bin/python tools/genome/structure_u3d.py web/genomics/structures/c_hex_a.pdb docs/papers/condensin-ii/figs/c_hex_a

Each chain is drawn as a smooth tube through its C-alpha atoms (Catmull-Rom spline, 8 points per residue). Chain A (the
kleisin region) is coloured by the model's per-residue confidence (pLDDT, B-factor column): red below 50, white at 70,
blue above 90. Chain B (the SMC head construct) is grey, and its artificial GGS linker is left out.
Writes <out>.u3d and <out>.views.json (centre and radius for the default camera).
"""
import json, sys
import numpy as np
import pymeshlab

def read(pdb):
    ch = {}
    for l in open(pdb):
        if l.startswith('ATOM') and l[12:16].strip() == 'CA':
            ch.setdefault(l[21], []).append((int(l[22:26]), l[17:20], np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]), float(l[60:66])))
    return ch

def spline(P, n=8):
    P = np.vstack([P[0], P, P[-1]]); out = []
    for i in range(1, len(P) - 2):
        p0, p1, p2, p3 = P[i - 1], P[i], P[i + 1], P[i + 2]
        for t in np.linspace(0, 1, n, endpoint=False):
            t2, t3 = t * t, t * t * t
            out.append(0.5 * ((2 * p1) + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t2 + (-p0 + 3 * p1 - 3 * p2 + p3) * t3))
    out.append(P[-2]); return np.array(out)

def tube(path, cols, r, sides=10):
    T = np.gradient(path, axis=0); T /= np.linalg.norm(T, axis=1, keepdims=True) + 1e-9
    a = np.array([0, 0, 1.0]) if abs(T[0][2]) < 0.9 else np.array([1.0, 0, 0])
    N = np.cross(T[0], a); N /= np.linalg.norm(N); V, C, F = [], [], []
    for i, (p, t) in enumerate(zip(path, T)):
        N = N - np.dot(N, t) * t; N /= np.linalg.norm(N) + 1e-9; B = np.cross(t, N)
        for k in range(sides):
            th = 2 * np.pi * k / sides; V.append(p + r * (np.cos(th) * N + np.sin(th) * B)); C.append(cols[i])
    for i in range(len(path) - 1):
        for k in range(sides):
            a0, a1 = i * sides + k, i * sides + (k + 1) % sides; b0, b1 = a0 + sides, a1 + sides
            F += [[a0, b0, a1], [a1, b0, b1]]
    return np.array(V), np.array(C), np.array(F)

def plddt_col(b):
    if b <= 70: t = max(0, (b - 50) / 20); return [1, 0.25 + 0.75 * t, 0.25 + 0.75 * t, 1]      # red -> white
    t = min(1, (b - 70) / 20); return [1 - 0.8 * t, 1 - 0.55 * t, 1, 1]                        # white -> blue

def main(pdb, out):
    ch = read(pdb); Vs, Cs, Fs, off = [], [], [], 0
    for cid, res in ch.items():
        if cid == 'B':   # drop the GGSGGS... linker: a run of G/S residues of length >= 8
            seq = ''.join('G' if r[1] == 'GLY' else 'S' if r[1] == 'SER' else 'x' for r in res); i = seq.find('GGSGGSGG')
            segs = [res[:i], res[i + 12:]] if i >= 0 else [res]
        else: segs = [res]
        for seg in segs:
            if len(seg) < 2: continue
            P = np.array([r[2] for r in seg]); path = spline(P)
            b = np.repeat([r[3] for r in seg], 8)[:len(path)]; b = np.concatenate([b, [b[-1]] * (len(path) - len(b))])
            cols = [plddt_col(x) for x in b] if cid == 'A' else [[0.62, 0.66, 0.64, 1]] * len(path)
            V, C, F = tube(path, cols, 1.25 if cid == 'A' else 0.9)
            Vs.append(V); Cs.append(C); Fs.append(F + off); off += len(V)
    V, C, F = np.vstack(Vs), np.vstack(Cs), np.vstack(Fs)
    ms = pymeshlab.MeshSet(); ms.add_mesh(pymeshlab.Mesh(vertex_matrix=V, face_matrix=F.astype(np.int32), v_color_matrix=C))
    ms.save_current_mesh(out + '.u3d', save_vertex_color=True)
    c = V.mean(0); r = float(np.linalg.norm(V - c, axis=1).max())
    json.dump({'centre': c.round(2).tolist(), 'radius': round(r, 1), 'vertices': len(V)}, open(out + '.views.json', 'w'))
    print(out + '.u3d', len(V), 'vertices; radius', round(r, 1))

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
