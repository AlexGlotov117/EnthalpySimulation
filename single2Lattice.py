import numpy as np
from ase.io import write
from ase.build import make_supercell
from ase import Atoms
from openff.units import unit

a_unit = 4.5181  # Angstroms (a and b axis length)
c_unit = 7.3560  # Angstroms (c axis length)

# Fractional coordinates for the 4 Oxygen atoms in the unit cell
# Explicitly using 1/3 and 2/3 for precision
O_frac_coords = np.array([
    [0.0000, 0.0000, 0.5000], 
    [0.0000, 0.0000, 0.0000], 
    [1/3, 2/3, 3/4], 
    [2/3, 1/3, 1/4]
])

# Replication factors to achieve a box size > 20.0 Angstroms (for 1.0 nm cutoff)
REPLICATION_FACTORS = (5, 5, 3) # Provides ~22.6 x 22.6 x 22.1 Angstroms

# --- 2. Build the Unit Cell ---
# Define the lattice vectors for the hexagonal unit cell
cell_vectors = np.array([
    [a_unit, 0.0, 0.0],
    [a_unit * np.cos(np.deg2rad(120)), a_unit * np.sin(np.deg2rad(120)), 0.0],
    [0.0, 0.0, c_unit]
])

# Create the unit cell using scaled_positions (fractional)
unit_cell_atoms = Atoms(
    symbols='O' * len(O_frac_coords),
    scaled_positions=O_frac_coords, 
    cell=cell_vectors,
    pbc=True
)

# --- 3. Generate the Supercell ---
ice_supercell = make_supercell(unit_cell_atoms, REPLICATION_FACTORS)

# Wrap coordinates to be inside the box, crucial for PDB
ice_supercell.wrap(pbc=True) 

# Calculate the final box dimensions
final_cell = ice_supercell.get_cell()
A = final_cell[0, 0]
B = final_cell[1, 1]
C = final_cell[2, 2]

# --- 4. Export to PDB ---
OUTPUT_PDB = 'solid_water_packed.pdb'
write(OUTPUT_PDB, ice_supercell, format='pdb')