#!/usr/bin/env python3
"""
Is the CA stereocentre of residue 143 correct (L-amino acid)?

If the model builder swapped the C-beta and H-alpha directions, the residue
would be D rather than L - a far more serious defect than a bad rotamer.

Chirality test: the signed volume of (N-CA) x (C-CA) . (CB-CA).
All 152 other residues in the same file give the reference sign for L.
"""
import numpy as np
from openmm.app import AmberPrmtopFile, AmberInpcrdFile
from openmm import unit

IN = "TRANSFER_to_JyotiPC"


def load(base, crd):
    prm = AmberPrmtopFile(f"{IN}/{base}.prmtop")
    xyz = np.array(AmberInpcrdFile(f"{IN}/{crd}").positions.value_in_unit(unit.angstrom))
    return prm.topology, xyz


def chirality(xyz, byname, cb):
    try:
        CA, N, C, CB = (xyz[byname[n]] for n in ("CA", "N", "C", cb))
    except KeyError:
        return None
    return float(np.dot(np.cross(N - CA, C - CA), CB - CA))


for tag, base, crd, cbname in (
        ("native",            "native_SOD1apo",      "native_SOD1apo.inpcrd",            "CB"),
        ("glycated as-built", "glyc_SOD1apo_MGH143", "glyc_SOD1apo_MGH143.inpcrd",       "C1"),
        ("glycated repaired", "glyc_SOD1apo_MGH143", "glyc_SOD1apo_MGH143_fixed.inpcrd", "C1")):
    top, xyz = load(base, crd)
    residues = list(top.residues())

    # reference: sign across all standard residues that have a real CB
    signs = []
    for r in residues:
        if r.index == 142:
            continue
        nm = {a.name: a.index for a in r.atoms()}
        if "CB" in nm:
            v = chirality(xyz, nm, "CB")
            if v is not None:
                signs.append(np.sign(v))
    ref = np.sign(np.mean(signs))

    r143 = residues[142]
    nm143 = {a.name: a.index for a in r143.atoms()}
    v143 = chirality(xyz, nm143, cbname)
    ok = np.sign(v143) == ref

    print(f"\n=== {tag} ===")
    print(f"  reference sign over {len(signs)} residues : {ref:+.0f} "
          f"(consistent: {np.mean(np.array(signs)==ref)*100:.0f}%)")
    print(f"  residue143 ({r143.name}) signed volume    : {v143:+.2f}  -> "
          f"{'L (correct)' if ok else 'D  ** INVERTED **'}")

    # where do CB-equivalent and HA sit relative to the native reference?
    if tag != "native":
        ntop, nxyz = load("native_SOD1apo", "native_SOD1apo.inpcrd")
        nres = list(ntop.residues())[142]
        nnm = {a.name: a.index for a in nres.atoms()}
        for lbl, gi, ni in (("C1  vs native CB", cbname, "CB"),
                            ("C1  vs native HA", cbname, "HA"),
                            ("HA  vs native CB", "HA",   "CB"),
                            ("HA  vs native HA", "HA",   "HA")):
            if gi in nm143 and ni in nnm:
                d = np.linalg.norm(xyz[nm143[gi]] - nxyz[nnm[ni]])
                print(f"    {lbl} : {d:6.2f} A")
