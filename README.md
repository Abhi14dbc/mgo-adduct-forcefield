# Methylglyoxal adduct force field and SOD1 glycation simulations

Force-field libraries for four methylglyoxal-derived advanced glycation end-products, the
structure-validation gate used to build systems from them, and the simulation and analysis code
for a study of site-specific glycation in apo-SH superoxide dismutase 1.

## What is here

```
forcefield/     AMBER-compatible mid-chain residue libraries
  MGH/            MG-H1, the hydroimidazolone formed at arginine
  ARP/            argpyrimidine
  CEA/            N-omega-carboxyethyl-arginine
  CEL/            N-epsilon-carboxyethyl-lysine
scripts/
  build/          tleap inputs and the RESP refit
  validation/     the structure-validation gate (two independent chirality tests)
  simulation/     replicate drivers and queue scripts
  analysis/       every script that produced a reported table or figure
topologies/     prmtop/inpcrd for the simulated systems
figures/        the figures as deposited, in vector form
```

## Using the residue libraries

Load alongside ff19SB and GAFF2 in tleap:

```
source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.opc
loadamberprep   forcefield/MGH/MGH.prepi
loadamberparams forcefield/MGH/MGH.frcmod
```

The residues are built as mid-chain units (HEAD = N, TAIL = C) so peptide bonds resolve in
either direction. Backbone atoms carry ff19SB protein types; side-chain atoms carry GAFF2 types.
`parmchk2` reports no missing bonded parameters for any of the four.

## Two things that will catch you

**Disulfide status cannot be read from residue names.** tleap forms a disulfide by an explicit
bond record and does not rename the participating cysteines, so a topology can carry an intact
disulfide while both residues are still named `CYS` rather than `CYX`. Read the bond list:

```python
ss = [(b[0].residue.resSeq, b[1].residue.resSeq) for b in top.bonds
      if b[0].element.symbol == "S" and b[1].element.symbol == "S"]
```

**MDTraj does not classify MG-H1 as protein.** A `select("protein")` atom selection silently
omits the adduct from a glycated system and shifts every subsequent residue index by one, so the
same index then addresses different residues in the modified and unmodified systems. Analyse
whole-system trajectories, or build the selection explicitly.

## Chirality

Adduct grafting by side-chain replacement can invert the CA stereocentre, and the result survives
energy minimisation and passes conventional bond-length and bond-angle checks. `scripts/validation`
tests every residue by two independent methods, the signed volume (N-CA)x(C-CA).(CB-CA) and the
HA-N-C-CB improper dihedral, and requires them to agree. The signed-volume test alone gives false
positives on strained but correctly configured residues. Automated repair of missing terminal
residues can also produce a D-configured terminus, so the check should be applied to rebuilt
termini and not only to grafted adducts.

## Reproducing the analysis

Trajectories are deposited separately on account of their size; see the record referenced in the
accompanying paper. With that record downloaded:

```
./reproduce.sh /path/to/trajectory/deposit
```

That runs both arms and prints the values to compare against the paper. It reads the
trajectories in place and writes nothing.

Analysis scripts contain the random seeds used, so bootstrap intervals reproduce exactly rather
than approximately. Permutation tests are exact enumerations, not sampled.

### Frame spacing

The deposited trajectories are at 200 ps per frame, the resolution at which every reported value
was computed. `analyse_stage3.py` strides its input by 20 because the working trajectories are
written at 10 ps, so on the deposit it needs to be told the input is already strided:

```
MGO_DIR=<deposit>/explicit_mono python analyse_phase1.py
MGO_DIR=<deposit> MGO_RAW_PS=200 MGO_STRIDE=1 python analyse_stage3.py
```

`MGO_RAW_PS` x `MGO_STRIDE` is the analysis frame spacing and must come to 200 ps. Occupancies
and hydrogen-bond counts are time averages and are insensitive to it, but switch rates are not:
a more finely sampled trajectory reports more transitions for identical physics, which is why
both solvent arms are analysed at a matched 200 ps.

Both scripts accept either layout. They try the working layout first and fall back to the
deposited one, so `MGO_DIR` is the only thing that needs setting.

## Requirements

AmberTools (tleap, antechamber, prepgen, parmchk2), OpenMM 8, MDTraj, NumPy, Matplotlib, RDKit.

## Citation

See `CITATION.cff`. Please cite the accompanying paper and this deposit.

## License

Code and libraries: MIT (see LICENSE). Figures: CC BY 4.0.
