#!/usr/bin/env python3
"""Window-sensitivity analysis for the preregistered primary endpoint.

The native system has not equilibrated at 300 ns: group-mean Gly61 occupancy
decays from ~95% in the first 50 ns to ~48% in the last, while the glycated
group is flat near 25%. The measured effect size therefore depends on which
window is analysed. This sweeps the window across the trajectory and reports
the endpoint in each, so the reader can see whether the CONCLUSION is robust
even though the MAGNITUDE is not.
"""
import os, itertools, math
import numpy as np
import mdtraj as md

D = "/scratch/a2148a01"
SITE, HIS48, GLY61, CYS57 = 142, 47, 60, 56
PREREG_CUT = 0.45
PREREG_ALPHA = 0.0042
FRAME_PS = 200.0
N_BOOT = 10000
RNG = np.random.default_rng(20260904)

GROUPS = {"native": [f"native_mono_rep{i}" for i in range(1, 11)],
          "glyc":   [f"glyc_mono_rep{i}"   for i in range(1, 11)]}

WINDOWS = [(0, 150), (50, 200), (100, 250), (150, 300), (0, 300)]


def perm_p(a, b):
    obs = abs(np.mean(b) - np.mean(a))
    pool = np.concatenate([a, b]); k = len(a); hits = tot = 0
    for idx in itertools.combinations(range(len(pool)), k):
        m = np.zeros(len(pool), bool); m[list(idx)] = True
        if abs(pool[~m].mean() - pool[m].mean()) >= obs - 1e-12:
            hits += 1
        tot += 1
    return hits / tot


def boot_ci(a, b):
    d = [np.mean(RNG.choice(b, len(b), True)) - np.mean(RNG.choice(a, len(a), True))
         for _ in range(N_BOOT)]
    return np.percentile(d, 2.5), np.percentile(d, 97.5)


# load once, slice many
dist = {g: [] for g in GROUPS}
for g, tags in GROUPS.items():
    for tag in tags:
        base = tag.rsplit("_rep", 1)[0]
        top = md.load_prmtop(os.path.join(D, base + ".prmtop"))
        dcd = os.path.join(D, "prod", tag, "full.dcd")
        if not os.path.exists(dcd):
            print("  !! missing " + dcd); continue
        assert top.residue(SITE).name in ("ARG", "MGH")
        assert top.residue(HIS48).name in ("HIS", "HIE", "HID", "HIP")
        assert top.residue(GLY61).name == "GLY"
        assert top.residue(CYS57).name == "CYS"
        t = md.load(dcd, top=top)
        pairs = [[SITE, HIS48], [SITE, CYS57], [SITE, GLY61]]
        dist[g].append(md.compute_contacts(t, pairs, scheme="closest-heavy")[0])
        print(f"  loaded {tag:20s} frames={t.n_frames}", flush=True)

print(__doc__)
for label, col in [("PREREGISTERED CORE composite (His48/Cys57/Gly61)", None),
                   ("component Gly61 alone", 2),
                   ("component His48 alone", 0)]:
    print("\n" + "=" * 96)
    print(label + f"   |  contact < {PREREG_CUT*10:.1f} A  |  prereg alpha = {PREREG_ALPHA}")
    print("=" * 96)
    print(f"{'window (ns)':>14s} {'native':>16s} {'glycated':>16s} {'delta':>9s} "
          f"{'95% CI':>20s} {'p':>9s}")
    for lo_ns, hi_ns in WINDOWS:
        i0, i1 = int(lo_ns * 1000 / FRAME_PS), int(hi_ns * 1000 / FRAME_PS)
        vals = {}
        for g in GROUPS:
            v = []
            for d in dist[g]:
                w = d[i0:i1]
                occ = (w < PREREG_CUT).mean(axis=0) * 100
                v.append(occ.mean() if col is None else occ[col])
            vals[g] = np.array(v)
        a, b = vals["native"], vals["glyc"]
        lo, hi = boot_ci(a, b); p = perm_p(a, b)
        flag = " *" if p <= PREREG_ALPHA else ""
        print(f"{lo_ns:>6d}-{hi_ns:<7d} {a.mean():8.2f}+-{a.std(ddof=1):5.2f} "
              f"{b.mean():8.2f}+-{b.std(ddof=1):5.2f} {b.mean()-a.mean():+9.2f} "
              f"[{lo:+7.2f},{hi:+7.2f}] {p:9.4f}{flag}", flush=True)
print("\n* = passes the preregistered Bonferroni alpha of 0.0042")
print("=" * 96)
