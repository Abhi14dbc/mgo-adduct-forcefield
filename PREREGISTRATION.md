---
title: "Pre-registered analysis plan and predictions"
timestamp: "2026-08-22 00:30 KST"
status: "SEALED - written before the data it concerns exists"
---

# Pre-registration

**Purpose.** The endpoints in this study were partly chosen after inspecting
data — a salt-bridge endpoint was tested and rejected, global RMSD was rejected
for lack of power, and the contact core was discovered by an unbiased scan and
then stratified after seeing replicate variance. Each step was documented and
defensible, but collectively that is a garden of forking paths. This document
fixes the analysis plan for data that **does not yet exist**, so the remaining
conditions are tested rather than described.

## 1. What exists at the time of writing (2026-08-22 00:30 KST)

| condition | status |
|---|---|
| glycated apo-SS, n=3 × 300 ns | **COMPLETE and verified** |
| native apo-SS, n=3 | rep2/rep3 at 100 ns, rep1 fresh run not started |
| apo-SH (native + glycated) | built, validated, **not started** |
| adduct series (ARP, CEA at Arg143) | built, validated, **not started** |
| CEL at Lys122 / Lys128 | built, validated, **not started** |
| dimer (native + glycated) | built, validated, **not started** |

Everything below concerns conditions in the **not started** rows.

## 2. Primary endpoint, fixed in advance

**Site-143 CORE contact-network integrity** — mean occupancy (heavy-atom
contact < 4.5 Å) over the three-residue core **His48 / Cys57 / Gly61**.

Secondary, in this order: hydrogen bonds donated by the site side chain;
SASA of the modified residue; electrostatic-loop Cα RMSF; protein-core Cα RMSF.

Reported but **not** used for inference: peripheral contacts (Thr58, Ala60,
Pro62, Val118, His120, Gly141) — replicate spread up to 73 points; β-sheet
content — GB-OBC is documented to reproduce secondary structure poorly relative
to explicit solvent [Nguyen 2013].

## 3. Statistical rules, fixed in advance

- The **replicate** is the statistical unit. n = 3 per condition. Frames are
  autocorrelated and are never pooled as n.
- **Range separation** (do the per-replicate ranges overlap?) is the primary
  evidence. Welch p and Cohen's d are reported alongside as descriptive.
- **Multiple comparisons will be reported.** With ~12 metrics, α = 0.05
  uncorrected is not a claim of significance. Bonferroni α = 0.0042 will be
  stated, and no result will be called significant on an uncorrected p alone.
- **Analysis windows are per-replicate**, set by explicit convergence checking,
  never a fixed equilibration discard. Justification: one native replicate
  required ~80 ns to converge, and a fixed 10–100 ns window inverted the sign of
  the effect.
- Any endpoint added after seeing these data will be labelled **post hoc** in
  the manuscript.

## 4. PREDICTIONS

Falsifiable, made before the data exist. Each states what would refute it.

### P1 — apo-SH shows a LARGER glycation effect than apo-SS
*Basis:* disulfide reduction is the dominant SOD1 destabiliser and apo-SH has
the highest fibrillar-aggregation propensity [Toichi 2013]; it is also the form
MGO reacts with preferentially [Polykretis 2020]. If glycation and disulfide
reduction act on the same axis, their effects should compound.
*Measure:* Δ(core integrity, native → glycated) in apo-SH vs the same Δ in apo-SS.
*Refuted if:* the apo-SH difference is equal to or smaller than the apo-SS one.

### P2 — the glycated dimer is less stable than the native dimer
*Basis:* lowering net charge increases aggregation propensity, and non-interface
mutations destabilise the SOD1 dimer allosterically [Farrokhzad 2024]; an
allosteric network vulnerable to perturbation is documented [Wells 2021].
MG-H1 lowers net charge by −1 per subunit (verified: −12 e → −14 e).
*Measure:* buried interface area, interface hydrogen-bond count, inter-subunit
centre-of-mass distance.
*Refuted if:* interface metrics are indistinguishable between conditions.

### P3 — Cys57 is the conduit
*Basis:* Arg143's conformation is sensitive to residue 57 and thereby governs
homodimer propensity [Wright 2020]; the Arg143–Cys57 contact is measurably
weakened by glycation in the completed apo-SS data.
*Measure:* if P2 holds, dimer destabilisation should correlate with loss of the
site-143–Cys57 contact across replicates.
*Refuted if:* the dimer is destabilised without any change at Cys57, or Cys57
weakens without dimer effect.

### P4 — the three arginine adducts differ in magnitude, not in kind
*Basis:* MG-H1, argpyrimidine and CEA all neutralise the guanidinium identically
(all verified −7 e), so electrostatics is held constant and only shape and
hydrogen-bonding capacity vary. Argpyrimidine is bulkier and aromatic; CEA is
open-chain.
*Measure:* core integrity and site H-bond count across the four systems.
*Refuted if:* the adducts differ in the DIRECTION of their effect, which would
indicate the mechanism is not simply charge neutralisation plus sterics.

### P5 — CEL at Lys122/Lys128 is structurally less disruptive than MG-H1 at Arg143
*Basis:* in S100A12, CEL modification left oligomerisation unaffected and
RAGE-binding only slightly impaired, whereas arginine modification at R21 was
"the major cause of MGO-induced impairment" [S100A12 study]. Both lysines are
2–3× more solvent-exposed than Arg143 (1.28 and 1.66 vs 0.65 nm²), so their
side chains are less structurally engaged to begin with.
*Refuted if:* CEL produces an effect comparable to or larger than MG-H1.
*Note:* CEL changes charge by −2 (Lys⁺ → carboxylate⁻) versus −1 for the
arginine adducts, so a larger effect would implicate charge magnitude over
structural engagement — also an informative outcome.

### P6 — no sequential-damage propagation (NULL prediction)
*Basis:* already tested on the completed apo-SS data. Three distal lysines
exceeded replicate scatter but changed in **opposite directions** (Lys122 −29 %,
Lys136 +26 %, Lys128 +13 %), and nothing survived multiple-comparison
correction across 14 sites. Restated here as a prediction for the new conditions.
*Prediction:* the other conditions will likewise show redistribution of
accessibility rather than a systematic increase.
*Refuted if:* distal Arg/Lys exposure increases systematically and survives
correction.

## 5. Committed reporting

- All metrics computed will be reported, including nulls.
- The endpoint-selection history (salt bridge tested and rejected; global RMSD
  rejected; core/periphery stratified after seeing variance) will appear in the
  methods, not be omitted.
- Predictions that fail will be reported as failed.
- This document will be deposited with the data.

---

*Written 2026-08-22 00:30 KST, after the glycated apo-SS set completed and
before native 300 ns, apo-SH, the adduct series, CEL and dimer data existed.*
