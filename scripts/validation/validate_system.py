#!/usr/bin/env python3
"""
Final acceptance check for a built system.
usage: validate_system.py <prmtop> <inpcrd> [expected_atoms]

Gate criteria
  - zero residues confirmed D by BOTH chirality tests
  - every MGH: L, and CA-C1 bond in 1.40-1.65 A
  - integer net charge
  - minimised potential energy negative and physical
"""
import sys
import numpy as np
from openmm.app import AmberPrmtopFile, AmberInpcrdFile
from openmm import unit, app
import openmm

prmtop, inpcrd = sys.argv[1], sys.argv[2]
expect = int(sys.argv[3]) if len(sys.argv) > 3 else None

prm = AmberPrmtopFile(prmtop)
crd = AmberInpcrdFile(inpcrd)
xyz = np.array(crd.positions.value_in_unit(unit.angstrom))
top = prm.topology
residues = list(top.residues())
fail = []


def dih(p0, p1, p2, p3):
    b0, b1, b2 = p0-p1, p2-p1, p3-p2
    b1n = b1/np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n)*b1n
    w = b2 - np.dot(b2, b1n)*b1n
    return np.degrees(np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))


print(f"system      : {prmtop}")
print(f"atoms {top.getNumAtoms()}   residues {len(residues)}")
if expect and top.getNumAtoms() != expect:
    fail.append(f"atom count {top.getNumAtoms()} != expected {expect}")

rows = []
for r in residues:
    nm = {a.name: a.index for a in r.atoms()}
    cb = "CB" if "CB" in nm else ("C1" if "C1" in nm else None)
    if not cb or not all(k in nm for k in ("N", "CA", "C", "HA")):
        continue
    CA = xyz[nm["CA"]]
    vol = float(np.dot(np.cross(xyz[nm["N"]]-CA, xyz[nm["C"]]-CA), xyz[nm[cb]]-CA))
    imp = dih(xyz[nm["HA"]], xyz[nm["N"]], xyz[nm["C"]], xyz[nm[cb]])
    rows.append((r, nm, cb, vol, imp))

rv = np.sign(np.median([v for *_, v, _ in rows]))
ri = np.sign(np.median([i for *_, i in rows]))
badD = [(r.name, r.index+1) for r, nm, cb, v, i in rows
        if np.sign(v) != rv and np.sign(i) != ri]
print(f"D-residues (both methods)  : {len(badD)} {badD[:6]}")
if badD:
    fail.append(f"{len(badD)} D residues")

mgh = [(r, nm, cb, v, i) for r, nm, cb, v, i in rows if r.name == "MGH"]
print(f"MGH residues               : {[r.index+1 for r, *_ in mgh]}")
for r, nm, cb, v, i in mgh:
    b = np.linalg.norm(xyz[nm["CA"]]-xyz[nm[cb]])
    ok = v > 0 and i > 0 and 1.40 < b < 1.65
    print(f"   MGH{r.index+1:<4d} volume {v:+6.2f}  improper {i:+6.1f}  "
          f"CA-C1 {b:.3f} A   {'OK' if ok else '** FAIL **'}")
    if not ok:
        fail.append(f"MGH{r.index+1}")

system = prm.createSystem(nonbondedMethod=app.CutoffNonPeriodic,
                          nonbondedCutoff=1.8*unit.nanometer,
                          constraints=None, implicitSolvent=app.OBC2)
nb = next(system.getForce(k) for k in range(system.getNumForces())
          if isinstance(system.getForce(k), openmm.NonbondedForce))
q = sum(nb.getParticleParameters(k)[0].value_in_unit(unit.elementary_charge)
        for k in range(system.getNumParticles()))
integ = openmm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond,
                                        0.002*unit.picosecond)
ctx = openmm.Context(system, integ, openmm.Platform.getPlatformByName("CPU"))
ctx.setPositions(crd.positions)
e0 = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
openmm.LocalEnergyMinimizer.minimize(ctx, maxIterations=2000)
e1 = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
print(f"net charge                 : {q:+.3f} e")
print(f"PE pre-minimisation        : {e0:,.1f} kJ/mol")
print(f"PE post-minimisation       : {e1:,.1f} kJ/mol")
if abs(q - round(q)) > 0.01:
    fail.append(f"non-integer charge {q:.3f}")
if e1 > 0:
    fail.append(f"minimised energy positive ({e1:.3g})")

print()
print("GATE: " + ("PASS" if not fail else "FAIL -> " + "; ".join(fail)))
sys.exit(1 if fail else 0)
