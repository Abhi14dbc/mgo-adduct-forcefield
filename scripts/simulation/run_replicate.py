#!/usr/bin/env python3
"""
Replicate MD for the MGO-glycation study (apo-SOD1, implicit solvent).

Settings are FIXED to match the existing replicate-1 runs - do not change them,
or the replicates will not be comparable.

Usage:
    python3 run_replicate.py <prmtop> <inpcrd> <outdir> <seed> [steps]

Example:
    python3 run_replicate.py native_SOD1apo.prmtop native_SOD1apo.inpcrd  native_rep2  2002  50000000

Auto-resumes if <outdir>/state.chk exists, so it survives interruption.
On resume the DCD is trimmed back to the checkpoint step first, so restarting
never leaves duplicated//orphaned frames past the checkpoint.
"""
import os, sys, time

from openmm import app, unit, LangevinMiddleIntegrator, Platform
from openmm.app import (AmberPrmtopFile, AmberInpcrdFile, Simulation,
                        StateDataReporter, DCDReporter, CheckpointReporter, PDBFile)

# ---- fixed protocol (must match replicate 1) ----
TEMP_K      = 300.0
FRICTION    = 1.0            # /ps
TIMESTEP_PS = 0.002          # 2 fs
CUTOFF_NM   = 1.8
REPORT      = 5000           # steps
CHECKPOINT  = 25000          # steps
MINIMIZE    = 2000           # iterations


def trim_dcd(dcd_path, prmtop_f, n_keep):
    """Cut traj.dcd back to n_keep frames (frames written after the last
    checkpoint are invalid once we rewind to that checkpoint)."""
    if not os.path.exists(dcd_path) or os.path.getsize(dcd_path) == 0:
        return
    try:
        import mdtraj as md
    except ImportError:
        print("WARNING: mdtraj unavailable, cannot trim DCD on resume", flush=True)
        return
    top = md.load_prmtop(prmtop_f)
    traj = md.load_dcd(dcd_path, top=top)
    if len(traj) <= n_keep:
        print(f"DCD has {len(traj)} frames, checkpoint expects {n_keep} - no trim needed",
              flush=True)
        return
    print(f"trimming DCD {len(traj)} -> {n_keep} frames to match checkpoint", flush=True)
    traj[:n_keep].save_dcd(dcd_path)


def main():
    if len(sys.argv) < 5:
        print(__doc__); return 2
    prmtop_f, inpcrd_f, outdir, seed = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    steps = int(sys.argv[5]) if len(sys.argv) > 5 else 50_000_000
    os.makedirs(outdir, exist_ok=True)

    prm = AmberPrmtopFile(prmtop_f)
    crd = AmberInpcrdFile(inpcrd_f)
    system = prm.createSystem(nonbondedMethod=app.CutoffNonPeriodic,
                              nonbondedCutoff=CUTOFF_NM * unit.nanometer,
                              constraints=app.HBonds,
                              implicitSolvent=app.OBC2)
    integ = LangevinMiddleIntegrator(TEMP_K * unit.kelvin,
                                     FRICTION / unit.picosecond,
                                     TIMESTEP_PS * unit.picosecond)
    integ.setRandomNumberSeed(seed)          # reproducible, distinct per replicate
    # Platform choice only - it does not alter the force field or any protocol
    # setting above. MD_PLATFORM=CPU forces the original CPU path.
    want = os.environ.get("MD_PLATFORM", "CUDA")
    sim = None
    for name, props in ((want, {"Precision": "mixed"} if want in ("CUDA", "OpenCL") else {}),
                        ("CPU", {})):
        try:
            sim = Simulation(prm.topology, system, integ,
                             Platform.getPlatformByName(name), props or None)
            break
        except Exception as ex:
            print(f"platform {name} unavailable ({type(ex).__name__}: {ex})", flush=True)
    if sim is None:
        sim = Simulation(prm.topology, system, integ)
    sim.context.setPositions(crd.positions)

    kj = unit.kilojoule_per_mole
    chk = os.path.join(outdir, "state.chk")
    dcd = os.path.join(outdir, "traj.dcd")
    t0 = time.time()
    print(f"{outdir}: {system.getNumParticles()} atoms | seed {seed} | "
          f"platform {sim.context.getPlatform().getName()}", flush=True)

    append = False
    if os.path.exists(chk) and os.path.getsize(chk) > 1000:
        sim.loadCheckpoint(chk); append = True
        print(f"resumed from checkpoint at step {sim.currentStep}", flush=True)
        trim_dcd(dcd, prmtop_f, sim.currentStep // REPORT)
    else:
        e0 = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kj)
        print(f"initial PE = {e0:,.1f} kJ/mol", flush=True)
        sim.minimizeEnergy(maxIterations=MINIMIZE)
        e1 = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kj)
        print(f"minimized  = {e1:,.1f} kJ/mol (dropped {e0-e1:,.1f})", flush=True)
        sim.context.setVelocitiesToTemperature(TEMP_K * unit.kelvin, seed)

    sim.reporters.append(StateDataReporter(os.path.join(outdir, "md.log"), REPORT,
                                           step=True, potentialEnergy=True, temperature=True,
                                           speed=True, progress=True, totalSteps=steps,
                                           elapsedTime=True, append=append))
    sim.reporters.append(DCDReporter(dcd, REPORT, append=append))
    sim.reporters.append(CheckpointReporter(chk, CHECKPOINT))

    remaining = steps - sim.currentStep
    print(f"running {remaining:,} more steps (target {steps:,} = {steps*2/1e6:.0f} ns)...", flush=True)
    sim.step(remaining)

    with open(os.path.join(outdir, "final.pdb"), "w") as f:
        PDBFile.writeFile(sim.topology,
                          sim.context.getState(getPositions=True).getPositions(), f)
    print(f"DONE {outdir} | {steps:,} steps | wall {time.time()-t0:.0f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
