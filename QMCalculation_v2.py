# %%
# .QM\Scripts\Activate.ps1
# eval "$(/home/aglotov/miniconda3/bin/conda shell.bash hook)"

# Repo contains scripts used to calculate enthalpy for gaseous stand alone molecule, bulk randomly sorted liquid phase molecules, and bulk crystal lattice solid phase molecules using a combination of quantum mechanics and molecular dynamics.

from openmm.app import *
from openmm import *
import openmm.unit as unit

import numpy as np
import pandas as pd
import mdtraj as md
from openff.toolkit.typing.engines.smirnoff import ForceField as OpenFFForceField
from openff.toolkit.topology import Molecule, Topology


# %%
# --- 1. Define Simulation Parameters and Input Data ---

# Parameters
TEMPERATURE = 294.14 * unit.kelvin
#TEMPERATURE = 250.0 * unit.kelvin
PRESSURE = 1.0 * unit.bar
SIMULATION_STEPS = 50000              # 1 ns of simulation (500,000 steps * 2 fs/step)
# N_MOLECULES = 1000                      # Number of molecules in the simulation box
# TARGET_DENSITY = 0.73 * unit.gram/unit.milliliter # Approximate liquid NH3 density

# %%
monomer_molecule = Molecule.from_pdb_and_smiles('water.pdb', 'O')

unique_molecules = [monomer_molecule]

pdb = PDBFile('liquid_water_packed.pdb')

openff = OpenFFForceField('openff-2.1.0.offxml') 

topology = Topology.from_openmm(pdb.topology, unique_molecules=unique_molecules)
system = openff.create_openmm_system(topology)

forces = system.getForces()
nb_force = [f for f in forces if isinstance(f, NonbondedForce)][0]

nb_force.setNonbondedMethod(NonbondedForce.PME)
nb_force.setCutoffDistance(1.0 * unit.nanometer)

integrator = LangevinMiddleIntegrator(TEMPERATURE, 1/unit.picosecond, 0.002*unit.picoseconds) # Smaller step (2fs) for safety

simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

# Set Box Vectors based on the solid PDB (Critical for PBC)
box_vectors = pdb.topology.getPeriodicBoxVectors()
simulation.context.setPeriodicBoxVectors(*box_vectors)

# Run Minimization (Crucial for crystals)
print("Minimizing energy...")
simulation.minimizeEnergy()

barostat = MonteCarloBarostat(PRESSURE, TEMPERATURE, 25) 
system.addForce(barostat)

simulation.context.reinitialize(preserveState=True)

simulation.reporters.append(DCDReporter('output.dcd', 1000))
simulation.reporters.append(StateDataReporter('output_liquid.log', 1000, 
    step=True, potentialEnergy=True, temperature=True, 
    volume=True, totalEnergy=True
))

print(f"Starting Production Run for {SIMULATION_STEPS} steps...")
simulation.step(SIMULATION_STEPS)
print("Simulation finished.")
