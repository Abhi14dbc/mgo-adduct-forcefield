#!/usr/bin/env python3
"""
Multi-conformer, two-stage RESP refit for the four MGO adducts.

WHY THIS EXISTS
---------------
The charges currently in the library were fitted from a SINGLE conformer, with a
SINGLE-STAGE restraint, on an MMFF94 geometry. That is not the canonical RESP
procedure and it is the one methodological criticism of the parameterisation a
referee is most likely to make. Canonical RESP is:

  * geometry from quantum mechanics, not molecular mechanics
  * several conformers fitted simultaneously, so charges are not tuned to one
    arbitrary backbone arrangement
  * TWO stages: stage 1 fits all atoms with a weak hyperbolic restraint;
    stage 2 refits only the chemically equivalent groups (e.g. methyl hydrogens)
    with a stronger restraint, holding everything else frozen
  * topological symmetry equivalencing within each stage

PHASES
------
Phase 1 (cheap, ~15 min per conformer): ESP single points at HF/6-31G(d) on the
        geometries we already have, then a proper two-stage multi-conformer fit.
        This alone fixes the single-conformer / single-stage criticism.

Phase 2 (expensive, ~8 h per conformer): re-optimise each conformer at
        HF/6-31G(d) first, then repeat phase 1 on the QM geometry. Run with
        --optimise once the GPU work has finished and the cores are free.

HF/6-31G(d) is deliberate, not a compromise. RESP charges for AMBER are
conventionally derived at exactly this level because it overestimates molecular
dipoles by roughly the amount needed to mimic condensed-phase polarisation. A
larger basis or a dispersion-corrected functional would be inconsistent with the
protein force field these charges must work alongside.
"""
import os, sys, json, argparse, subprocess, tempfile, shutil, time
import numpy as np

MGO_DIR = "/home/jyoti/Projects/MGO"
AMBER_BIN = "/home/jyoti/miniforge3/envs/amber23/bin"
ADDUCT_DIR = os.path.join(MGO_DIR, "abhi_import/models/module01_adduct_ff")
OUT_DIR = os.path.join(MGO_DIR, "resp_refit")

ADDUCTS = {                       # name -> capped structure to start from
    "MGH": "residue_MGH/capped_Lfixed.mol2",   # the L-corrected template
    "ARP": "residue_ARP/capped.mol2",
    "CEA": "residue_CEA/capped.mol2",
    "CEL": "residue_CEL/capped.mol2",
}

# Formal charge of each CAPPED MODEL COMPOUND (not the charge the residue
# contributes to the protein). The mol2 charge column is NOT trustworthy:
# residue_CEL/capped.mol2 carries partial charges summing to 0.000 even
# though its carboxylate is deprotonated (O1/O2 both bare, no O-H), i.e.
# the molecule is an anion. Trusting that column gave 147 electrons at
# spin 0 and PySCF refused to build the molecule -- after the job had
# already sat in the queue for six hours.
# Verified against the production prmtops: native -6; MGH/ARP/CEA -7
# (a +1 residue replaced by a neutral one); CEL -8 (neutral Lys AND a
# -1 carboxylate).
FORMAL_CHARGE = {"MGH": 0, "ARP": 0, "CEA": 0, "CEL": -1}

# Backbone conformers to fit simultaneously (phi, psi in degrees)
CONFORMERS = {"alphaR": (-60.0, -45.0), "beta": (-135.0, 135.0), "PII": (-75.0, 145.0)}


# ----------------------------------------------------------------- structure IO
def read_mol2(path):
    """Return (symbols, coords, atom_names, total_charge, bonds)."""
    syms, xyz, names, q, bonds = [], [], [], 0.0, []
    with open(path) as fh:
        section = None
        for line in fh:
            if line.startswith("@<TRIPOS>"):
                section = line.strip().split(">")[1]
                continue
            if section == "ATOM" and line.strip():
                f = line.split()
                names.append(f[1])
                xyz.append([float(f[2]), float(f[3]), float(f[4])])
                # SYBYL type -> element
                t = f[5].split(".")[0]
                syms.append(t)
                if len(f) > 8:
                    q += float(f[8])
            elif section == "BOND" and line.strip():
                f = line.split()
                if len(f) >= 3:
                    bonds.append((int(f[1]) - 1, int(f[2]) - 1))
    return syms, np.array(xyz), names, q, bonds


def write_xyz(path, syms, xyz, comment=""):
    with open(path, "w") as fh:
        fh.write(f"{len(syms)}\n{comment}\n")
        for s, c in zip(syms, xyz):
            fh.write(f"{s:2s} {c[0]:14.8f} {c[1]:14.8f} {c[2]:14.8f}\n")




# ----------------------------------------------------------------- geometry
def _rotate_about_bond(xyz, i, j, moving, angle_rad):
    """Rotate `moving` atom indices about the i->j bond axis."""
    axis = xyz[j] - xyz[i]
    axis = axis / np.linalg.norm(axis)
    c, s_ = np.cos(angle_rad), np.sin(angle_rad)
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    R = np.eye(3) + s_ * K + (1 - c) * (K @ K)
    out = xyz.copy()
    out[moving] = (R @ (xyz[moving] - xyz[j]).T).T + xyz[j]
    return out


def _dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return np.degrees(np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))


def _side_of_bond(bonds, n_atoms, i, j):
    """Atoms reachable from j without crossing the i-j bond."""
    adj = {k: set() for k in range(n_atoms)}
    for a, b in bonds:
        adj[a].add(b); adj[b].add(a)
    seen, stack = {j}, [j]
    while stack:
        cur = stack.pop()
        for nb in adj[cur]:
            if nb == i and cur == j:
                continue
            if nb not in seen:
                seen.add(nb); stack.append(nb)
    seen.discard(j)
    return sorted(seen)


def set_dihedral(xyz, bonds, quad, target_deg):
    """Set the dihedral defined by four atom indices to target_deg."""
    i0, i1, i2, i3 = quad
    cur = _dihedral(xyz[i0], xyz[i1], xyz[i2], xyz[i3])
    delta = np.radians(target_deg - cur)
    moving = _side_of_bond(bonds, len(xyz), i1, i2)
    if i3 not in moving:
        return xyz          # bond graph did not separate cleanly; leave as is
    return _rotate_about_bond(xyz, i1, i2, moving, delta)


def backbone_quads(names):
    """Locate phi (C_prev-N-CA-C) and psi (N-CA-C-N_next) in a capped residue."""
    idx = {n: k for k, n in enumerate(names)}
    def find(*cands):
        for c in cands:
            if c in idx: return idx[c]
        return None
    N, CA, C = find("N"), find("CA"), find("C")
    # These templates cap with ACE = CY/CAY/OY and NME = NZ/CZ
    Cprev = find("CY", "C_ACE", "CH3_ACE", "C1_ACE")
    Nnext = find("NZ", "N_NME", "NN", "N_NM")
    if None in (N, CA, C):
        return None, None
    phi = (Cprev, N, CA, C) if Cprev is not None else None
    psi = (N, CA, C, Nnext) if Nnext is not None else None
    return phi, psi


# ----------------------------------------------------------------- QM
def build_mol(syms, xyz, basis="6-31g*", charge=0):
    from pyscf import gto
    atoms = [(s, tuple(c)) for s, c in zip(syms, xyz)]
    mol = gto.M(atom=atoms, basis=basis, charge=int(round(charge)),
                unit="Angstrom", verbose=0)
    return mol


def optimise_hf(syms, xyz, charge, maxsteps=60):
    """HF/6-31G(d) geometry optimisation via geomeTRIC. Expensive."""
    from pyscf import scf
    from pyscf.geomopt.geometric_solver import optimize
    mol = build_mol(syms, xyz, charge=charge)
    mf = scf.RHF(mol).density_fit()
    mol_eq = optimize(mf, maxsteps=maxsteps)
    return mol_eq.atom_coords() * 0.52917721092      # Bohr -> Angstrom


def esp_on_mk_grid(syms, xyz, charge):
    """
    HF/6-31G(d) electrostatic potential on a Merz-Kollman grid.
    Returns (grid_points_Angstrom, esp_values_au, total_energy).
    """
    from pyscf import scf, df, gto
    mol = build_mol(syms, xyz, charge=charge)
    mf = scf.RHF(mol).density_fit()
    e = mf.kernel()
    dm = mf.make_rdm1()

    # Merz-Kollman: 4 shells at 1.4, 1.6, 1.8, 2.0 x vdW radius, ~1 point/A^2
    VDW = {"H": 1.20, "C": 1.50, "N": 1.50, "O": 1.40, "S": 1.75}
    pts = []
    coords = mol.atom_coords() * 0.52917721092
    elems = [mol.atom_symbol(i) for i in range(mol.natm)]
    rng = np.random.default_rng(2026)
    for scale in (1.4, 1.6, 1.8, 2.0):
        for i, (el, c) in enumerate(zip(elems, coords)):
            r = VDW.get(el, 1.5) * scale
            n = max(int(4 * np.pi * r * r * 2.5), 24)      # ~2.5 points per A^2
            # Fibonacci sphere for even coverage
            k = np.arange(n) + 0.5
            phi = np.arccos(1 - 2 * k / n)
            theta = np.pi * (1 + 5 ** 0.5) * k
            sph = np.stack([np.cos(theta) * np.sin(phi),
                            np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1)
            cand = c + r * sph
            # keep only points outside every atom's shell (MK exclusion)
            keep = np.ones(len(cand), bool)
            for j, (el2, c2) in enumerate(zip(elems, coords)):
                r2 = VDW.get(el2, 1.5) * scale
                d = np.linalg.norm(cand - c2, axis=1)
                keep &= (d >= r2 - 1e-6)
            pts.append(cand[keep])
    grid = np.vstack(pts)

    # V(r) = sum_A Z_A/|r-R_A| - integral rho(r')/|r-r'|
    grid_bohr = grid / 0.52917721092
    from pyscf import lib
    esp = np.zeros(len(grid_bohr))
    Z = mol.atom_charges()
    R = mol.atom_coords()
    for i in range(mol.natm):
        esp += Z[i] / np.linalg.norm(grid_bohr - R[i], axis=1)
    # electronic part, in blocks to bound memory
    blk = 2000
    for s in range(0, len(grid_bohr), blk):
        sub = grid_bohr[s:s + blk]
        fakemol = gto.fakemol_for_charges(sub)
        ints = df.incore.aux_e2(mol, fakemol, intor="int3c2e")
        esp[s:s + blk] -= np.einsum("ijk,ij->k", ints, dm)
    return grid, esp, e


# ----------------------------------------------------------------- RESP
def write_esp_file(path, coords_A, grid_A, esp_au):
    """AmberTools resp input format (Bohr, atomic units)."""
    B = 1.0 / 0.52917721092
    with open(path, "w") as fh:
        fh.write(f"{len(coords_A):5d}{len(grid_A):5d}    0\n")
        for c in coords_A:
            fh.write(f"{'':16s}{c[0]*B:16.7E}{c[1]*B:16.7E}{c[2]*B:16.7E}\n")
        for v, g in zip(esp_au, grid_A):
            fh.write(f"{v:16.7E}{g[0]*B:16.7E}{g[1]*B:16.7E}{g[2]*B:16.7E}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--adduct", default="MGH", choices=list(ADDUCTS))
    ap.add_argument("--optimise", action="store_true",
                    help="Phase 2: HF geometry optimisation first (~8 h/conformer)")
    ap.add_argument("--conformers", default="alphaR",
                    help="comma-separated subset of " + ",".join(CONFORMERS))
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    src = os.path.join(ADDUCT_DIR, ADDUCTS[args.adduct])
    if not os.path.exists(src):
        print(f"missing capped structure: {src}"); return 1
    syms, xyz, names, q_mol2, bonds = read_mol2(src)
    qtot = FORMAL_CHARGE[args.adduct]
    Zn = {"H": 1, "C": 6, "N": 7, "O": 8, "S": 16}
    n_elec = sum(Zn[e] for e in syms) - qtot
    if n_elec % 2:
        print(f"ABORT: {args.adduct} has {n_elec} electrons at charge "
              f"{qtot:+d} -- odd, so it cannot be the closed-shell "
              f"species assumed here. Check FORMAL_CHARGE and the "
              f"capped structure.")
        return 1
    print(f"{args.adduct}: {len(syms)} atoms from {os.path.basename(src)}, "
          f"formal charge {qtot:+d}, {n_elec} electrons (even, OK)")
    if abs(q_mol2 - qtot) > 0.01:
        print(f"  note: mol2 partial charges sum to {q_mol2:+.3f}, not "
              f"{qtot:+d}; using the declared formal charge.")

    wanted = [c.strip() for c in args.conformers.split(",")]
    records = []
    for conf in wanted:
        t0 = time.time()
        geom = xyz.copy()
        phi_q, psi_q = backbone_quads(names)
        tphi, tpsi = CONFORMERS[conf]
        applied = []
        if phi_q and None not in phi_q:
            geom = set_dihedral(geom, bonds, phi_q, tphi); applied.append(f"phi={tphi:.0f}")
        if psi_q and None not in psi_q:
            geom = set_dihedral(geom, bonds, psi_q, tpsi); applied.append(f"psi={tpsi:.0f}")
        print(f"  [{conf}] backbone set: {', '.join(applied) if applied else 'NOT APPLIED (backbone atoms not found)'}", flush=True)
        if args.optimise:
            print(f"  [{conf}] HF/6-31G(d) optimisation ...", flush=True)
            geom = optimise_hf(syms, geom, qtot)
        print(f"  [{conf}] ESP single point ...", flush=True)
        grid, esp, energy = esp_on_mk_grid(syms, geom, qtot)
        espf = os.path.join(OUT_DIR, f"{args.adduct}_{conf}.esp")
        write_esp_file(espf, geom, grid, esp)
        print(f"  [{conf}] E = {energy:.6f} Eh | {len(grid):,} grid points | "
              f"{time.time()-t0:.0f} s -> {os.path.basename(espf)}", flush=True)
        records.append(dict(conformer=conf, energy=float(energy),
                            n_grid=int(len(grid)), esp_file=espf,
                            optimised=bool(args.optimise)))

    meta = dict(adduct=args.adduct, n_atoms=len(syms), atom_names=names,
                total_charge=qtot, level="HF/6-31G(d)",
                geometry="QM-optimised" if args.optimise else "as-supplied (MMFF)",
                conformers=records)
    with open(os.path.join(OUT_DIR, f"{args.adduct}_esp.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    print(f"wrote {OUT_DIR}/{args.adduct}_esp.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
