# MGO glycation of apo-SOD1 — final protocol (rev. 2026-08-18)

## Systems
| | native | glycated |
|---|---|---|
| file | `native_SOD1apo` | `glyc_SOD1apo_MGH143` |
| atoms | 2187 | 2192 |
| residues | 153 (1 chain) | 153 (1 chain) |
| residue 143 | ARG (+1 e) | MGH / MG-H1 (0 e) |
| net charge | −6 e | −7 e |
| metals | none (apo) | none (apo) |
| Cys57–Cys146 | **disulfide intact** (2.04 Å) | disulfide intact |

Verified: the pair differs at exactly one residue. Net atom change C+3/H+1/O+1 =
Arg⁺ + MGO(C₃H₄O₂) − H₂O − H⁺ → MG-H1 (neutral). The neutral adduct matches the
literature pKa of the hydroimidazolone group (~3.6–4.6), i.e. uncharged at pH 7.

## Simulation settings — UNCHANGED (must match replicate 1)
- AMBER prmtop/inpcrd; implicit solvent **GB-OBC2**
- `nonbondedMethod=CutoffNonPeriodic`, cutoff **1.8 nm**
- `constraints=HBonds`
- **LangevinMiddleIntegrator**, 300 K, friction 1 ps⁻¹, timestep 2 fs
- energy minimisation 2000 iterations
- `setVelocitiesToTemperature(300 K, seed)`, distinct seed per replicate
- report + DCD every 5000 steps; checkpoint every 25000 steps
- 50,000,000 steps = **100 ns** per replicate

### Platform change (physics unchanged)
Runs on **CUDA / mixed precision** (GTX 1070) rather than CPU. CPU measured
3.1 ns/day (~32 d/replicate); CUDA gives **383 ns/day** (~6.3 h/replicate).
Cross-check: minimised energy CPU −22,901 kJ/mol vs GPU −22,936 kJ/mol.
Revert with `MD_PLATFORM=CPU`.

### Seeds
| replicate | seed |
|---|---|
| native_rep2 | 2002 |
| native_rep3 | 2003 |
| glyc_rep2 | 3002 |
| glyc_rep3 | 3003 |

⚠ Confirm replicate 1 (other machine) did not use these values.

## Analysis protocol — REVISED

### Rejected endpoint: salt-bridge loss
Rabbani et al. (2021) motivate salt-bridge loss as the consequence of MG-H1
formation. **Measured and rejected**: Arg143 shows *0.0%* salt-bridge occupancy
in apo-SS SOD1; the nearest carboxylate oxygen is 10.1 Å from the guanidinium.
There is no salt bridge to lose in this system.

### Rejected endpoint: global Cα-RMSD
Underpowered. Two native replicates from identical coordinates diverge to
1.86 Å Cα, so the between-condition effect would have to exceed the
within-condition noise floor. It will not at 100 ns.

### Primary endpoint: site-143 contact-network integrity
Unbiased contact scan of the Arg143 sidechain (heavy atom < 4.5 Å) gives a
well-defined, mechanistically loaded network:

| partner | occupancy | mean min-dist | relevance |
|---|---|---|---|
| Gly61 | 100.0% | 2.92 Å | dominant H-bond |
| His48 | 99.8% | 3.13 Å | Cu-site ligand |
| Cys57 | 99.7% | 3.17 Å | disulfide partner; dimerisation axis |
| Thr58 | 92.1% | 4.03 Å | |
| Val118 | 80.2% | 4.23 Å | |
| Gly141 | 79.0% | 4.11 Å | |
| Ala60 | 76.8% | 4.26 Å | |
| His120 | 69.8% | 4.33 Å | metal-site ligand |
| Pro62 | 37.1% | 4.61 Å | |

Sequence neighbours 142/144 excluded as trivially in contact.

Supports Wright et al. (2020): Arg143 "interposes between disulfide bond
residues Cys57 and Cys146". Cys57 contact is confirmed (99.7%); Cys146 is not
in direct contact (6.36 Å) in the apo-SS resting state. The literature-proposed
Asn53 H-bond is *absent* here (7.85 Å) — it applies to the hCCS-complexed state.

### Stratification — core vs periphery (decided on full-length native data)
Averaging all nine partners is **not** a usable endpoint: it mixes a
near-invariant core with a highly replicate-dependent periphery, and the
periphery dominates the variance.

| set | members | rep2 | rep3 | spread |
|---|---|---|---|---|
| **CORE** | His48, Cys57, Gly61 | 99.815% | 99.852% | **0.037 pt** |
| periphery | Thr58, Ala60, Pro62, Val118, His120, Gly141 | 72.50% | 55.65% | 16.9 pt |
| whole panel | all nine | 81.60% | 70.38% | 11.2 pt |

Worst individual offender: Gly141, 79.0% vs 6.1% (72.9-point spread).

**Primary endpoint = CORE integrity.** It is simultaneously the lowest-noise
metric (0.037-point replicate spread, vs a 1.9 Å noise floor for global Cα-RMSD)
and the mechanistically loaded one: Cys57 = disulfide/dimerisation axis,
His48 = Cu-site ligand, Gly61 = dominant H-bond.

Core sits near ceiling (~99.8%), so it can only be driven *down* — which matches
the directional hypothesis (charge removal weakens charge-assisted contacts).
Report it as an absolute drop in percentage points, not a ratio.

**Periphery is reported but must not carry inference at n=3.** It would need
substantially more sampling (µs-scale or many more replicates) to be usable.

### Secondary endpoints
- H-bond count from the residue-143 sidechain (D–A < 3.5 Å, D–H···A > 120°)
- SASA of residue 143
- **Cys57 SASA** — solvent exposure of the disulfide, proxy for the
  monomerisation/unfolding axis reported by Polykretis et al. (2020)
- Electrostatic-loop (121–142) Cα RMSF
- Global Cα-RMSD — reference only, not an endpoint

### Statistics
The **replicate** is the statistical unit. Frames within a trajectory are
autocorrelated; pooling them inflates n by ~10⁴ and manufactures significance.
Report effect size and per-replicate values; with n=3 a t-test is descriptive,
not confirmatory.

## Structural repair of the glycated model (2026-08-19)

The supplied `glyc_SOD1apo_MGH143.inpcrd` was **structurally invalid** and all
results derived from it (`glyc_rep2`, `glyc_rep3`, and glycated replicate 1 on
the other machine) must be discarded.

**Diagnosis.** Backbone RMSD to native was 0.000 Å (true matched pair), and all
MG-H1 internal bond lengths were normal — but the side chain was attached to CA
~103° off the true Cβ vector, burying it in the protein core:

| | native | glycated (as supplied) | repaired |
|---|---|---|---|
| single-point PE (kJ/mol) | −16,161 | **9.0×10¹³** | −13,646 |
| non-bonded pairs < 1.5 Å | 3 | **37** | 3 |
| worst overlap | 1.01 Å | **0.20 Å** (His120:CA···MGH143:H5) | 3.51 Å clearance |
| CB↔C1 displacement | — | **2.41 Å** | 0.05 Å |
| CD↔C3 displacement | — | **6.52 Å** | — |

χ-angle rotation cannot fix this: χ1 rotates *about* the CA→C1 bond and so can
never move C1. An exhaustive 13,824-rotamer scan found zero clash-free
placements from the broken geometry.

**Consequence for the earlier result.** The −95% "core network integrity" effect
was an artefact: minimisation ejected the adduct from inside His120, which is
why it ended up pinned against His120 (99.6% contact) and far from the native
pocket. The interim conclusion that MG-H1 geometrically cannot reach the native
site was likewise an artefact of scanning from the broken structure.

**Repair (coordinates only; prmtop and all parameters untouched).**
Rigid-body superposition of the intact side chain via (CA,C1,C2,C3) onto the
native (CA,CB,CG,CD) frame, proper rotation only (det = +1, chirality
preserved), stem fit RMSD 0.054 Å; then a scan of the two distal torsions
(χ3 = −180°, χ4 = +50°) for maximum clearance. Built by
`fix_glycated_structure.py` → `glyc_SOD1apo_MGH143_fixed.inpcrd`.

Repaired structure places the adduct in native-like contact range
(His48 3.76 Å, Cys57 3.51 Å, Gly61 4.24 Å vs native 3.1/3.2/2.9 Å) and
minimises to −23,432 kJ/mol, comparable to native (−22,901 / −22,979).

Production runs on the repaired structure: `glycfix_rep2` (seed 3002),
`glycfix_rep3` (seed 3003), launched 2026-08-19 10:51.

## RESULTS (n=3 per condition, completed 2026-08-20 05:44)

All six replicates verified: 10,000 frames each, final step 50,000,000, zero
duplicate steps, T = 300.1-300.3 K (sd ~5.7), no potential-energy drift.

### Analysis windows
Per-replicate, because a fixed window is unsafe here. `native_rep1` spent ~80 ns
in a displaced state before settling, so it is analysed over its converged
100-159 ns window; all others over 10-100 ns. Using the fixed 10-100 ns window
for rep1 gives 64.7% (transient) instead of 99.64%, and that artefact **inverts
the sign** of the apparent glycation effect (+10.1% instead of -10.2%).

### Four metrics separate completely (no overlap between conditions)

| metric | native (per-replicate) | glycated (per-replicate) | effect |
|---|---|---|---|
| H-bonds at site 143 | 2.56 / 3.17 / 3.82 | **1.82 / 1.09 / 1.51** | **-53.7%** |
| Electrostatic-loop RMSF (Å) | 0.96 / 0.89 / 0.72 | **1.61 / 1.46 / 1.41** | **+75%** |
| SASA residue 143 (nm²) | 0.62 / 0.65 / 0.50 | **0.80 / 0.90 / 0.75** | **+38%** |
| CORE integrity (%) | 99.61 / 99.82 / 99.85 | **87.8 / 96.9 / 83.9** | **-10.2%** |

Mechanistically coherent and consistent with Rabbani et al. (2021): removing the
guanidinium strips H-bonding capacity, the residue becomes more solvent-exposed,
the electrostatic loop loosens, and the core contact network weakens.

### Two metrics show NO effect (report as negatives)
- **Cys57 SASA**: native 0.181-0.209, glycated 0.139-0.196 nm² — overlapping.
  The -13.8% mean shift is driven by one replicate; not reliable.
- **Peripheral contacts**: native 55.6-72.5%, glycated 73.4-81.1% — overlapping
  and high-variance in both conditions.

### Effect DIRECTION is solid; MAGNITUDE is not
Native replicates agree to 0.25 points (99.61-99.85). Glycated replicates span
13 points (83.9-96.9) — 50x more scatter. Block traces show why:

    glycfix_rep2  98.6  98.1  98.3  97.8  93.0   still declining at 100 ns
    glycfix_rep3  97.7  46.8  87.1  98.8  96.9   large excursion, recovers

The glycated site is genuinely dynamic and has not settled within 100 ns.
Quote the effect as qualitative/lower-bound unless the runs are extended.
Extending the three glycated replicates to 300 ns costs ~1.6 days at 383 ns/day
and would pin the magnitude down; resume is automatic from each state.chk.

## Known limitations to state in the paper
1. **Disulfide state.** Polykretis et al. (2020) report MG reacts preferentially
   with the *disulfide-reduced* demetallated form. These models are apo with the
   Cys57–Cys146 disulfide **intact** (apo-SS). apo-SH is the more MGO-reactive
   species and is the natural third condition — but it is a different system and
   is not comparable to replicate 1.
2. **Site competition.** Monteiro Neto et al. (2023) emphasise Lys122/Lys128.
   Defence for arginine: MG is an arginine-directed glycating agent
   (Rabbani et al., 2021), and Arg143 is a reported preferential site in apo/Zn
   metalloforms (Polykretis et al., 2020; Wright et al., 2020).
3. **Timescale.** 100 ns addresses local network perturbation only. Claims about
   monomerisation, unfolding or aggregation require µs sampling or enhanced
   sampling and are out of scope for this design.
4. **GB radii.** The prmtops trigger OpenMM's "non-optimal GB parameters for GB
   model OBC2" warning (radius set is not mbondi2). Applies identically to both
   systems and to replicate 1, so the comparison is internally consistent.

## References
- Lei, M., Peng, F., & Ye, Y. (2016). *Front. Mol. Biosci.* 3, 55. https://doi.org/10.3389/fmolb.2016.00055
- Monteiro Neto, J. R., et al. (2023). *BBA Mol. Basis Dis.* 1869(8), 166835. https://doi.org/10.1016/j.bbadis.2023.166835
- Polykretis, P., Luchinat, E., & Boscaro, F. (2020). *Redox Biol.* 30, 101421. https://doi.org/10.1016/j.redox.2019.101421
- Rabbani, N., Xue, M., & Thornalley, P. J. (2021). *Glycoconj. J.* 38(3), 331–340. https://doi.org/10.1007/s10719-021-09980-0
- Wright, G. S. A., et al. (2020). *Angew. Chem. Int. Ed.* https://doi.org/10.1002/anie.202000451
- Dalilah, Y., & Simm, A. (2026). *Front. Mol. Biosci.* 13. https://doi.org/10.3389/fmolb.2026.1838713
