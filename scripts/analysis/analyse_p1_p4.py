#!/usr/bin/env python3
"""P1 and P4: the two preregistered predictions whose data was never analysed.

P1 -- apo-SH shows a LARGER glycation effect than apo-SS.
  Measure: Delta(core integrity, native -> glycated) in apo-SH vs the same
  Delta in apo-SS. Refuted if the apo-SH difference is equal to or smaller
  than the apo-SS one.

  Tested here at a MATCHED 300 ns for the first time. The earlier 100 ns test
  was uninformative because the native apo-SH endpoint saturated near 99%,
  leaving no room for a difference to appear. Native core occupancy declines
  substantially over 300 ns, so the saturation may be relieved -- the printed
  saturation check reports whether it is.

  The interaction is tested by permuting the native/glycated labels WITHIN each
  disulfide state independently, giving C(6,3)^2 = 400 arrangements and a
  two-sided floor of 0.005 -- twenty times finer than the 0.10 floor that binds
  each pairwise comparison at n = 3.

P4 -- the three arginine adducts differ in magnitude, not in kind.
  Measure: core integrity and site H-bond count across native, MG-H1,
  argpyrimidine and CEA at Arg143, all apo-SS at a matched 100 ns.
  Refuted if the adducts differ in the DIRECTION of their effect.

All conditions use the same endpoints, cutoff and statistics as the rest of the
study. n = 3 per group throughout: pairwise permutation floor is 0.10.

DISULFIDE STATE is asserted from the BOND LIST, not residue names: these
topologies name the bonded pair CYS rather than CYX, so a name check silently
passes on the wrong system.
"""
import os, sys, itertools, math
import numpy as np
import mdtraj as md

D = os.environ.get("MGO_DIR", "/home/jyoti/Projects/MGO")
STRIDE = 20
FRAME_PS = 10.0 * STRIDE
SITE, HIS48, GLY61, CYS57 = 142, 47, 60, 56
PREREG_CUT = 0.45
N_BOOT = 10000
RNG = np.random.default_rng(20260904)

T = {
    "native_SS": f"{D}/TRANSFER_to_JyotiPC/native_SOD1apo.prmtop",
    "MGH_SS":    f"{D}/TRANSFER_to_JyotiPC/glyc_SOD1apo_MGH143.prmtop",
    "ARP_SS":    f"{D}/sod1_adducts/glyc_SOD1apo_ARP143.prmtop",
    "CEA_SS":    f"{D}/sod1_adducts/glyc_SOD1apo_CEA143.prmtop",
    "native_SH": f"{D}/sod1_new/native_SH.prmtop",
    "glyc_SH":   f"{D}/sod1_new/glyc_SH.prmtop",
}


def has_disulfide(top):
    return any(b[0].element.symbol == "S" and b[1].element.symbol == "S"
               for b in top.bonds)


def analyse(top_key, tag, lo_ns, hi_ns, want_ss):
    top = md.load_prmtop(T[top_key])
    dcd = f"{D}/{tag}/traj.dcd"
    if not os.path.exists(dcd):
        print(f"  !! missing {dcd}"); return None

    assert top.residue(SITE).name in ("ARG", "MGH", "ARP", "CEA"), \
        f"{tag} site is {top.residue(SITE).name}"
    assert top.residue(HIS48).name in ("HIS", "HIE", "HID", "HIP")
    assert top.residue(GLY61).name == "GLY"
    assert top.residue(CYS57).name in ("CYS", "CYX")
    assert has_disulfide(top) == want_ss, \
        f"{tag}: disulfide={has_disulfide(top)}, expected {want_ss}"

    t = md.load(dcd, top=top, stride=STRIDE)
    i0, i1 = int(lo_ns * 1000 / FRAME_PS), int(hi_ns * 1000 / FRAME_PS)
    assert t.n_frames >= i1, f"{tag}: {t.n_frames} frames, need {i1}"
    t = t[i0:i1]

    d = md.compute_contacts(t, [[SITE, HIS48], [SITE, CYS57], [SITE, GLY61]],
                            scheme="closest-heavy")[0]
    occ = (d < PREREG_CUT).mean(axis=0) * 100

    at = t.topology.atom
    per_frame = [sum(1 for don, _, acc in fr
                     if at(don).residue.index == SITE
                     or at(acc).residue.index == SITE)
                 for fr in md.wernet_nilsson(t[::5], periodic=False)]

    return dict(core_occ=float(occ.mean()), core_his=float(occ[0]),
                core_cys=float(occ[1]), core_gly=float(occ[2]),
                hbond=float(np.mean(per_frame)), n_frames=t.n_frames)


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


def interaction_p(sh_n, sh_g, ss_n, ss_g):
    """Permute condition labels WITHIN each disulfide state; 400 arrangements."""
    obs = abs((np.mean(sh_g) - np.mean(sh_n)) - (np.mean(ss_g) - np.mean(ss_n)))
    pool_sh = np.concatenate([sh_n, sh_g]); pool_ss = np.concatenate([ss_n, ss_g])
    k = len(sh_n); hits = tot = 0
    for i_sh in itertools.combinations(range(len(pool_sh)), k):
        m1 = np.zeros(len(pool_sh), bool); m1[list(i_sh)] = True
        d_sh = pool_sh[~m1].mean() - pool_sh[m1].mean()
        for i_ss in itertools.combinations(range(len(pool_ss)), k):
            m2 = np.zeros(len(pool_ss), bool); m2[list(i_ss)] = True
            d_ss = pool_ss[~m2].mean() - pool_ss[m2].mean()
            if abs(d_sh - d_ss) >= obs - 1e-12:
                hits += 1
            tot += 1
    return hits / tot, obs


print(__doc__)

# ---------------------------------------------------------------- P1
P1 = [("SS_native", "native_SS",  [f"native_rep{i}"   for i in (1, 2, 3)], True),
      ("SS_glyc",   "MGH_SS",     [f"glycfix_rep{i}"  for i in (1, 2, 3)], True),
      ("SH_native", "native_SH",  [f"nativeSH_rep{i}" for i in (1, 2, 3)], False),
      ("SH_glyc",   "glyc_SH",    [f"glycSH_rep{i}"   for i in (1, 2, 3)], False)]

print("\n" + "=" * 100)
print("P1 -- apo-SH vs apo-SS glycation effect, MATCHED 300 ns, window 150-300 ns")
print("=" * 100)
r1 = {}
for label, key, tags, want_ss in P1:
    r1[label] = []
    for tag in tags:
        r = analyse(key, tag, 150.0, 300.0, want_ss)
        if r:
            r1[label].append(r)
            print(f"  {label:10s} {tag:16s} core={r['core_occ']:5.1f}%  "
                  f"his48={r['core_his']:5.1f}%  cys57={r['core_cys']:5.1f}%  "
                  f"gly61={r['core_gly']:5.1f}%  hb={r['hbond']:4.2f}", flush=True)

for key, name in [("core_occ", "CORE composite (%)"), ("hbond", "Site H-bonds")]:
    ss_n = np.array([r[key] for r in r1["SS_native"]])
    ss_g = np.array([r[key] for r in r1["SS_glyc"]])
    sh_n = np.array([r[key] for r in r1["SH_native"]])
    sh_g = np.array([r[key] for r in r1["SH_glyc"]])
    d_ss, d_sh = ss_g.mean() - ss_n.mean(), sh_g.mean() - sh_n.mean()
    p_int, obs = interaction_p(sh_n, sh_g, ss_n, ss_g)
    print(f"\n  {name}")
    print(f"    apo-SS  native {ss_n.mean():7.2f}+-{ss_n.std(ddof=1):5.2f}   "
          f"glycated {ss_g.mean():7.2f}+-{ss_g.std(ddof=1):5.2f}   delta {d_ss:+7.2f}"
          f"   (pairwise p={perm_p(ss_n, ss_g):.2f})")
    print(f"    apo-SH  native {sh_n.mean():7.2f}+-{sh_n.std(ddof=1):5.2f}   "
          f"glycated {sh_g.mean():7.2f}+-{sh_g.std(ddof=1):5.2f}   delta {d_sh:+7.2f}"
          f"   (pairwise p={perm_p(sh_n, sh_g):.2f})")
    verdict = "SUPPORTED" if abs(d_sh) > abs(d_ss) else "REFUTED"
    print(f"    |delta_SH| = {abs(d_sh):.2f}  vs  |delta_SS| = {abs(d_ss):.2f}"
          f"   -->  P1 {verdict}")
    print(f"    interaction |d_SH - d_SS| = {obs:.2f}, permutation p = {p_int:.4f}"
          f"  (floor 0.005 over 400 arrangements)")
    if key == "core_occ":
        print(f"    SATURATION CHECK: native occupancies  apo-SS "
              f"{[round(x,1) for x in ss_n]}   apo-SH {[round(x,1) for x in sh_n]}")
        print("    (the 100 ns test failed because native apo-SH sat near 99%)")

# ---------------------------------------------------------------- P4
P4 = [("native", "native_SS", [f"native_rep{i}"  for i in (1, 2, 3)]),
      ("MG-H1",  "MGH_SS",    [f"glycfix_rep{i}" for i in (1, 2, 3)]),
      ("ARP",    "ARP_SS",    [f"arp143_rep{i}"  for i in (1, 2, 3)]),
      ("CEA",    "CEA_SS",    [f"cea143_rep{i}"  for i in (1, 2, 3)])]

print("\n" + "=" * 100)
print("P4 -- three arginine adducts at Arg143, apo-SS, MATCHED 100 ns, window 50-100 ns")
print("=" * 100)
r4 = {}
for label, key, tags in P4:
    r4[label] = []
    for tag in tags:
        r = analyse(key, tag, 50.0, 100.0, True)
        if r:
            r4[label].append(r)
            print(f"  {label:7s} {tag:16s} core={r['core_occ']:5.1f}%  "
                  f"hb={r['hbond']:4.2f}", flush=True)

for key, name in [("core_occ", "CORE composite (%)"), ("hbond", "Site H-bonds")]:
    print(f"\n  {name}")
    base = np.array([r[key] for r in r4["native"]])
    print(f"    {'native':8s} {base.mean():7.2f}+-{base.std(ddof=1):5.2f}")
    signs = []
    for label in ("MG-H1", "ARP", "CEA"):
        v = np.array([r[key] for r in r4[label]])
        d = v.mean() - base.mean()
        lo, hi = boot_ci(base, v)
        signs.append(np.sign(d))
        print(f"    {label:8s} {v.mean():7.2f}+-{v.std(ddof=1):5.2f}   delta {d:+7.2f}"
              f"  [{lo:+7.2f},{hi:+7.2f}]   p={perm_p(base, v):.2f}")
    same = len(set(signs)) == 1
    print(f"    directions {'ALL AGREE' if same else 'DISAGREE'}  -->  "
          f"P4 {'supported' if same else 'REFUTED'} on this endpoint")

print("\n" + "=" * 100)
print("n = 3 per group: pairwise p cannot fall below 0.10. The P1 interaction")
print("test permutes within each disulfide state and has a floor of 0.005.")
print("=" * 100)
