#!/usr/bin/env python3
"""CEL at Lys122 / Lys128 vs unmodified native, implicit GB-OBC2, n = 3, 100 ns.

Arai 1987 identified Lys122 and Lys128 as the glycation sites whose modification
inactivates SOD1; Monteiro Neto 2023 names them with Arg143 as the three
residues whose modification most affects SOD1 structure and function. All three
help stabilise the electrostatic loop.

Arg143 is UNMODIFIED in all three conditions here (confirmed by assertion), so
this asks a specific question: does glycation at Lys122 or Lys128 perturb the
site-143 core contact network REMOTELY? If it does, the effect reported for
MG-H1 at Arg143 is part of a general electrostatic-loop destabilisation. If it
does not, that effect is specific to modification of Arg143 itself.

MATCHED WINDOW. The CEL runs are 100 ns; the native apo-SH runs are 300 ns. To
compare like with like we take the FINAL 50 ns OF THE FIRST 100 ns from every
condition, so all three are analysed over the same interval of simulated time.
This is a shorter window than the 150 ns used for the explicit-solvent work and
the values are NOT comparable with it.

n = 3 per group: the exact permutation floor is 2/C(6,3) = 0.10, so NO ENDPOINT
HERE CAN REACH SIGNIFICANCE. Read the deltas and confidence intervals.

Global metrics use the 151 protein residues common to all three topologies
(excluding both modified positions), so the numbers are strictly comparable.
"""
import os, sys, itertools, math
import numpy as np
import mdtraj as md

D = os.environ.get("MGO_DIR", "/home/jyoti/Projects/MGO")
STRIDE = 20              # traj.dcd is 10 ps/frame -> 200 ps
WIN_LO_NS, WIN_HI_NS = 50.0, 100.0     # final 50 ns of the first 100 ns
FRAME_PS = 10.0 * STRIDE
SITE, HIS48, GLY61, CYS57 = 142, 47, 60, 56
LYS122, LYS128 = 121, 127
PREREG_CUT = 0.45
N_BOOT = 10000
RNG = np.random.default_rng(20260904)

CONDS = {
    "native": (os.path.join(D, "sod1_new", "native_SH.prmtop"),
               [f"nativeSH_rep{i}" for i in range(1, 4)]),
    "cel122": (os.path.join(D, "sod1_cel", "cel122.prmtop"),
               [f"cel122_rep{i}" for i in range(1, 4)]),
    "cel128": (os.path.join(D, "sod1_cel", "cel128.prmtop"),
               [f"cel128_rep{i}" for i in range(1, 4)]),
}


def analyse(top_f, tag):
    top = md.load_prmtop(top_f)
    dcd = os.path.join(D, tag, "traj.dcd")
    if not os.path.exists(dcd):
        print("  !! missing " + dcd); return None

    # the site must be UNMODIFIED arginine in every condition, or the question
    # being asked is not the one described above
    assert top.residue(SITE).name == "ARG", tag + " site is " + top.residue(SITE).name
    assert top.residue(HIS48).name in ("HIS", "HIE", "HID", "HIP")
    assert top.residue(GLY61).name == "GLY"
    assert top.residue(CYS57).name == "CYS"

    t = md.load(dcd, top=top, stride=STRIDE)
    i0, i1 = int(WIN_LO_NS * 1000 / FRAME_PS), int(WIN_HI_NS * 1000 / FRAME_PS)
    assert t.n_frames >= i1, f"{tag}: {t.n_frames} frames, need {i1}"
    t = t[i0:i1]

    pairs = [[SITE, HIS48], [SITE, CYS57], [SITE, GLY61]]
    d = md.compute_contacts(t, pairs, scheme="closest-heavy")[0]
    occ = (d < PREREG_CUT).mean(axis=0) * 100

    common = t.topology.select(f"protein and not resid {LYS122} and not resid {LYS128}")
    prot = t.atom_slice(common)
    bb = prot.topology.select("backbone")
    ss = md.compute_dssp(prot, simplified=True)

    return dict(core_his=float(occ[0]), core_cys=float(occ[1]), core_gly=float(occ[2]),
                core_occ=float(occ.mean()),
                rmsd=float(md.rmsd(prot, prot, 0, atom_indices=bb).mean() * 10),
                rg=float(md.compute_rg(prot).mean() * 10),
                beta=float((ss == "E").mean() * 100),
                n_frames=t.n_frames, n_common=prot.n_residues)


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


print(__doc__)
res = {}
for cond, (top_f, tags) in CONDS.items():
    res[cond] = []
    for tag in tags:
        r = analyse(top_f, tag)
        if r:
            res[cond].append(r)
            print(f"  {tag:16s} frames={r['n_frames']:5d} common_res={r['n_common']}  "
                  f"core={r['core_occ']:5.1f}%  his48={r['core_his']:5.1f}%  "
                  f"cys57={r['core_cys']:5.1f}%  gly61={r['core_gly']:5.1f}%", flush=True)

METRICS = [("core_occ", "Site-143 CORE composite (%)"),
           ("core_his", "  His48 (%)"),
           ("core_cys", "  Cys57 (%)"),
           ("core_gly", "  Gly61 (%)"),
           ("rmsd",     "Backbone RMSD (A)"),
           ("rg",       "Radius of gyration (A)"),
           ("beta",     "Beta-sheet content (%)")]

for other in ("cel122", "cel128"):
    n_a, n_b = len(res["native"]), len(res[other])
    if n_a < 2 or n_b < 2:
        print(f"\n!! {other}: too few replicates"); continue
    floor = 2 / math.comb(n_a + n_b, n_a)
    print("\n" + "=" * 100)
    print(f"native (n={n_a})  vs  {other} (n={n_b})   |   window {WIN_LO_NS:.0f}-{WIN_HI_NS:.0f} ns"
          f"   |   permutation floor p = {floor:.3g}")
    print("=" * 100)
    print(f"{'endpoint':30s} {'native':>15s} {other:>15s} {'delta':>9s} {'95% CI':>20s} {'p':>8s} {'sep':>4s}")
    for key, name in METRICS:
        a = np.array([r[key] for r in res["native"]])
        b = np.array([r[key] for r in res[other]])
        lo, hi = boot_ci(a, b); p = perm_p(a, b)
        sep = "yes" if (max(a) < min(b) or max(b) < min(a)) else "no"
        print(f"{name:30s} {a.mean():7.2f}+-{a.std(ddof=1):6.2f} "
              f"{b.mean():7.2f}+-{b.std(ddof=1):6.2f} {b.mean()-a.mean():+9.2f} "
              f"[{lo:+7.2f},{hi:+7.2f}] {p:8.4f} {sep:>4s}")

print("\n" + "=" * 100)
print("Arg143 is UNMODIFIED in all three conditions. Any change in the site-143")
print("core here is a REMOTE effect of modifying Lys122 or Lys128, not a direct")
print("consequence of modifying the site itself.")
print("n = 3 per group: p cannot fall below 0.10 regardless of effect size.")
print("=" * 100)
