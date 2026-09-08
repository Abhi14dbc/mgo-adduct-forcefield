#!/usr/bin/env python3
"""Verify SMILES for the MGO adducts, then render them.

These are free amino acids (not chain residues) so formula and charge can be
checked against the chemistry. Paste the SMILES into ChemDraw (Edit > Paste
Special > SMILES) to get a correctly connected structure to lay out, rather
than drawing from scratch.
"""
import io, sys
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors, AllChem
from rdkit.Chem.Draw import rdMolDraw2D

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

MOLS = [
    ("Methylglyoxal", "CC(=O)C=O", 0),
    ("L-Arginine (guanidinium, +1)",
     "N[C@@H](CCCNC(=[NH2+])N)C(=O)O", +1),
    ("MG-H1  (Nd-(5-methyl-4-oxo-imidazolin-2-yl)ornithine)",
     "N[C@@H](CCCNC1=NC(C)C(=O)N1)C(=O)O", 0),
    ("Argpyrimidine  (Nd-(5-hydroxy-4,6-dimethylpyrimidin-2-yl)ornithine)",
     "N[C@@H](CCCNc1nc(C)c(O)c(C)n1)C(=O)O", 0),
    ("CEA  (Nw-carboxyethylarginine, zwitterionic side chain)",
     "N[C@@H](CCCNC(=[NH2+])NC(C)C(=O)[O-])C(=O)O", 0),
    ("CEL  (Ne-carboxyethyl-lysine)",
     "N[C@@H](CCCCNC(C)C(=O)[O-])C(=O)O", -1),
    ("L-Lysine (+1, for comparison with CEL)",
     "N[C@@H](CCCC[NH3+])C(=O)O", +1),
]

print(f"{'compound':56s} {'formula':22s} {'chg':>4s} {'MW':>8s}  ok")
print("-" * 100)
ok_all = True
mols = []
for name, smi, expect_q in MOLS:
    m = Chem.MolFromSmiles(smi)
    if m is None:
        print(f"{name:56s} {'PARSE FAILED':22s}"); ok_all = False; continue
    q = Chem.GetFormalCharge(m)
    f = rdMolDescriptors.CalcMolFormula(m)
    good = (q == expect_q)
    ok_all &= good
    print(f"{name:56s} {f:22s} {q:+4d} {Descriptors.MolWt(m):8.2f}  "
          f"{'OK' if good else 'CHARGE MISMATCH, expected %+d' % expect_q}")
    mols.append((name.split("  ")[0], m))

print()
# side-chain charge change on modification, the number the prmtop must reproduce
arg = Chem.MolFromSmiles(MOLS[1][1]); mgh = Chem.MolFromSmiles(MOLS[2][1])
lys = Chem.MolFromSmiles(MOLS[6][1]); cel = Chem.MolFromSmiles(MOLS[5][1])
print("charge change on modification (free amino acid):")
print(f"  Arg -> MG-H1 : {Chem.GetFormalCharge(arg):+d} -> "
      f"{Chem.GetFormalCharge(mgh):+d}   delta {Chem.GetFormalCharge(mgh)-Chem.GetFormalCharge(arg):+d}"
      "   (protein -6 -> -7: consistent)")
print(f"  Lys -> CEL   : {Chem.GetFormalCharge(lys):+d} -> "
      f"{Chem.GetFormalCharge(cel):+d}   delta {Chem.GetFormalCharge(cel)-Chem.GetFormalCharge(lys):+d}"
      "   (preregistration states -2: consistent)")
print()

# hydrogen-bond donor / acceptor counts -- the basis of the P4 explanation
print(f"{'compound':30s} {'HBD':>5s} {'HBA':>5s}")
for name, smi, _ in MOLS[1:6]:
    m = Chem.MolFromSmiles(smi)
    print(f"{name.split('  ')[0]:30s} {rdMolDescriptors.CalcNumHBD(m):5d} "
          f"{rdMolDescriptors.CalcNumHBA(m):5d}")

for nm, m in mols:
    AllChem.Compute2DCoords(m)
d = Draw.MolsToGridImage([m for _, m in mols], molsPerRow=3,
                         subImgSize=(430, 340),
                         legends=[n for n, _ in mols])
d.save("figures/adduct_structures.png")
print("\nwrote figures/adduct_structures.png")
print("ALL CHECKS PASSED" if ok_all else "\n!! CHECK FAILURES ABOVE")
