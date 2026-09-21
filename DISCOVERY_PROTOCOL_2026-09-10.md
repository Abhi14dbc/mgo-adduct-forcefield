---
title: "apo-SH discovery run: what will be done and what will be reported"
written: "2026-09-10, before the trajectories were started"
status: "committed in advance of the data"
---

# Purpose

The site-143 contact core (His48, Cys57, Gly61) was defined by an unbiased contact
scan run on **apo-SS** trajectories in August 2026, and the confirmatory experiment
is **apo-SH**. Section 4.5 of the manuscript reports that the same scan applied to
the confirmatory apo-SH trajectories does not reproduce the stratification: Val118,
Thr58, Ala60 and His120 are all more occupied and less variable than any core member.

That scan is descriptive only, because it runs on the data the endpoint was tested
on. This run generates apo-SH explicit-solvent trajectories that are held out from
every confirmatory statistic, so the question can be answered without circularity:

**If the contact set had been defined in apo-SH explicit solvent, would it have been
His48, Cys57 and Gly61?**

# What will be run

Three native apo-SH explicit-solvent replicates, 300 ns each, OPC water, identical
protocol to the confirmatory set:

- `prod.py --system native_mono --rep 101 / 102 / 103 --ns 300`
- seeds by the existing formula, `crc32("native_mono") % 100000 + rep x 7919`,
  so replicate indices 101-103 cannot collide with the confirmatory 1-10
- same equilibrated starting state (`native_mono_equil.xml`) as the confirmatory set
- native only; defining a contact set does not require a glycated arm

Cost: 3 x 11.6 GPU-hours = approximately 35 GPU-hours, about 1.1 SRU.

**These three replicates will not enter any confirmatory statistic in the
manuscript.** They exist solely to define a contact set.

# Analysis, fixed now

1. The unbiased scan of `probe_site143_contacts.py`: every residue contacting the
   site side chain, heavy-atom criterion 4.5 A, final 150 ns window, ranked by mean
   occupancy with replicate standard deviation reported alongside.
2. Residues will be classified by the same two criteria the apo-SS scan used, high
   mean occupancy and low replicate variance. Sequence neighbours (142, 144) are
   excluded as trivially contacted.
3. Whatever set that yields will be formed into a composite and evaluated on the
   confirmatory n = 10 per group, with the same exact permutation test.

# What will be reported, whichever way it falls

- The apo-SH ranking in full, and where His48, Cys57 and Gly61 place in it.
- The effect size and p value of the apo-SH-derived composite on the confirmatory
  data, **including if it is smaller or non-significant** relative to the
  prespecified endpoint.
- If the apo-SH-derived set is dominated by contacts near saturation, that will be
  stated as an explanation, not used to set the result aside. Val118 sits at
  97.8 +/- 1.8% in the confirmatory native data and has almost no dynamic range;
  a contact that cannot fall cannot show a decrease. This is a foreseeable outcome
  and is being written down before the data exist precisely so that it cannot be
  produced afterwards as a rescue.

The prespecified endpoint remains the primary result regardless of what this run
shows. This run bears on how that endpoint should be described, not on whether the
reported difference is real.
