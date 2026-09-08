#!/usr/bin/env python3
"""Phase 1 production: one replicate, 300 ns, explicit OPC solvent.
Starts from the 5 ns equilibrated state, so no production time is spent
re-equilibrating.
REPLICATE INDEPENDENCE: every replicate begins from the same equilibrated
coordinates but gets fresh velocities from a per-replicate seed, and the
integrator gets that seed too. Without this the replicates would be identical
and n=5 would be a fiction. crc32, not hash() -- hash() is salted per process,
so a restart would silently get a different seed.
300 ns is the floor, not a choice: in the implicit-solvent work 3 of 12
endpoints changed between 100 and 300 ns, two reversing sign."""
import argparse, os, time, zlib
from openmm import app, unit, XmlSerializer
from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform
from mdtraj.reporters import DCDReporter as MDTrajDCD
import mdtraj as md

ap = argparse.ArgumentParser()
ap.add_argument("--system", required=True)
ap.add_argument("--rep", type=int, required=True)
ap.add_argument("--ns", type=float, default=300.0)
ap.add_argument("--dt", type=float, default=4.0)
ap.add_argument("--dir", default="/scratch/a2148a01")
a = ap.parse_args()

tag = f"{a.system}_rep{a.rep}"
out = os.path.join(a.dir, "prod", tag)
os.makedirs(out, exist_ok=True)
prm = os.path.join(a.dir, a.system + ".prmtop")
eq = os.path.join(a.dir, a.system + "_equil.xml")
for f in (prm, eq):
    if not os.path.exists(f):
        raise SystemExit("missing " + f)

top = app.AmberPrmtopFile(prm)
system = top.createSystem(nonbondedMethod=app.PME,
                          nonbondedCutoff=1.0 * unit.nanometer,
                          constraints=app.HBonds, rigidWater=True,
                          hydrogenMass=4.0 * unit.amu)
system.addForce(MonteCarloBarostat(1 * unit.bar, 300 * unit.kelvin, 25))

seed = (zlib.crc32(a.system.encode()) % 100000) + a.rep * 7919
integ = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                 a.dt * unit.femtoseconds)
integ.setRandomNumberSeed(seed)
sim = app.Simulation(top.topology, system, integ,
                     Platform.getPlatformByName("CUDA"), {"Precision": "mixed"})
with open(eq) as fh:
    sim.context.setState(XmlSerializer.deserialize(fh.read()))
sim.context.setVelocitiesToTemperature(300 * unit.kelvin, seed)

n = int(a.ns * 1e6 / a.dt)
prot_idx = md.Topology.from_openmm(top.topology).select("protein")
sim.reporters.append(MDTrajDCD(os.path.join(out, "prot.dcd"),
                               int(10000 / a.dt), atomSubset=prot_idx))
sim.reporters.append(app.DCDReporter(os.path.join(out, "full.dcd"),
                                     int(200000 / a.dt)))
sim.reporters.append(app.StateDataReporter(
    os.path.join(out, "md.log"), int(n / 3000), step=True, time=True,
    potentialEnergy=True, temperature=True, density=True, volume=True,
    progress=True, remainingTime=True, speed=True, totalSteps=n, separator=","))
sim.reporters.append(app.CheckpointReporter(
    os.path.join(out, "state.chk"), int(5e6 / a.dt)))

t0 = time.time()
print(f"START {tag}  seed={seed}  {a.ns:.0f} ns  {n} steps", flush=True)
sim.step(n)
sim.saveState(os.path.join(out, "final.xml"))
with open(os.path.join(out, "final.pdb"), "w") as fh:
    app.PDBFile.writeFile(sim.topology,
        sim.context.getState(getPositions=True).getPositions(), fh)
el = time.time() - t0
print(f"DONE  {tag}  {el/3600:.2f} h  {a.ns/(el/86400):.1f} ns/day", flush=True)
