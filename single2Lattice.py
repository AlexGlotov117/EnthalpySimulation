# from ase.io import read, write
# from ase.build import make_supercell
# import numpy as np

# target_molecules=100
# input_file="ice_28unitCell.pdb"

# unit_cell = read(input_file)

# initial_oxygen_atoms = sum(1 for atom in unit_cell if atom.symbol == 'O')

# if initial_oxygen_atoms == 0:
#     # Fallback if the input structure is unconventional or incorrect
#     print("Warning: Could not detect Oxygen atoms to count molecules. Assuming 28 molecules.")
#     molecules_per_cell = 28 
# else:
#     molecules_per_cell = initial_oxygen_atoms

# required_cells = target_molecules / molecules_per_cell

# # Calculate the side length N required for cubic replication
# N = int(np.ceil(required_cells**(1/3)))

# # Ensure N is at least 1
# N = max(1, N)
# replications = (1, 1, 2)

# # Calculate the final count and replication matrix
# final_molecule_count = N**3 * molecules_per_cell
# replication_matrix = np.diag(replications) # [[N, 0, 0], [0, N, 0], [0, 0, N]]

# print(f"Initial molecules (per unit cell): {molecules_per_cell}")
# print(f"Required replication factor (N^3): {N**3}")
# print(f"Targeting: {replications[0]}x{replications[1]}x{replications[2]} supercell.")
# print(f"Resulting in approximately {final_molecule_count} water molecules.")

# # 3. Create the Supercell
# # make_supercell automatically applies the translation and boundary conditions correctly.
# supercell = make_supercell(unit_cell, replication_matrix)


# # 4. Write the final structure to a PDB file
# # ASE ensures the output PDB contains the correct CRYST1 record for the new, larger box.
# write('solid_ice_packed.pdb', supercell, format='proteindatabank')

# # 5. Output Summary
# print("-" * 50)
# print(f"SUCCESS: Generated supercell with {len(supercell)} atoms.")
# print(f"Total Molecules: {final_molecule_count}")
# print(f"Output saved")

# # Print the resulting box dimensions (ASE stores them in the .cell attribute)
# cell_matrix = supercell.get_cell()
# a = np.linalg.norm(cell_matrix[0])
# b = np.linalg.norm(cell_matrix[1])
# c = np.linalg.norm(cell_matrix[2])
# print(f"New Simulation Box Dimensions: A={a:.3f}, B={b:.3f}, C={c:.3f} Angstroms")
# print("-" * 50)

import re
import os

def xyz_to_pdb_water(xyz_filename, pdb_filename, cell_a, cell_b, cell_c, n_total_atoms):
    """
    Converts a water ice structure in XYZ format (O-H-H order assumed) 
    to a standard PDB file including CRYST1 and CONECT records.
    """
    
    # --- PDB FORMAT DEFINITIONS ---
    RESIDUE_NAME = "HOH"
    
    # PDB ATOM record format (fixed width)
    # The Element Symbol (%2s) must be at the very end.
    ATOM_FMT = "ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n"
    
    # CRYST1 record format (cell dimensions and angles)
    CRYST1_FMT = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n"
    # --- END FORMAT DEFINITIONS ---

    if not os.path.exists(xyz_filename):
        print(f"ERROR: Input file not found at {xyz_filename}")
        return

    with open(xyz_filename, 'r') as f:
        # Skip the first two lines of XYZ (atom count and comment)
        f.readline()
        f.readline()
        
        atom_index = 1
        residue_index = 1
        pdb_lines = []
        
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Use split() without arguments to handle any mix of whitespace (spaces or tabs)
            parts = line.split()
            
            # Check for minimum number of parts (Atom_Type X Y Z)
            if len(parts) < 4:
                continue
                
            atom_type = parts[0].upper()
            try:
                # 1. Coordinate Extraction Fix: Explicitly take indices 1, 2, 3
                x, y, z = map(float, [parts[1], parts[2], parts[3]])
            except ValueError:
                print(f"Skipping line due to coordinate parsing error: '{line}'")
                continue

            # --- CORE LOGIC: O-H1-H2 GROUPING ---
            
            # Position within the O-H-H triplet (0 for O, 1 for H1, 2 for H2)
            atom_in_mol_pos = (atom_index - 1) % 3

            element_symbol = atom_type # Element is always the first part of the line (O or H)
            
            if atom_in_mol_pos == 0:
                # This is the Oxygen atom. Start a new residue index.
                atom_name = "O"
                residue_index = (atom_index - 1) // 3 + 1
            elif atom_in_mol_pos == 1:
                # This is the first Hydrogen atom.
                atom_name = "H1"
            elif atom_in_mol_pos == 2:
                # This is the second Hydrogen atom.
                atom_name = "H2"
            else:
                atom_name = element_symbol
            
            # Standard PDB values
            occupancy = 1.00
            b_factor = 0.00
            
            # Format the ATOM record
            # Crucial: Ensure atom_name is correctly aligned for the PDB format
            pdb_line = ATOM_FMT % (
                atom_index, 
                atom_name.ljust(4), # Left-justified ATOM NAME field
                RESIDUE_NAME, 
                ' ', # Chain ID
                residue_index, 
                x, y, z, 
                occupancy, 
                b_factor, 
                element_symbol.ljust(2) # ELEMENT field (right-most)
            )
            pdb_lines.append(pdb_line)
            atom_index += 1

    # Ensure the correct number of atoms was processed
    if atom_index - 1 != n_total_atoms:
        print(f"WARNING: Expected {n_total_atoms} atoms but processed {atom_index - 1}.")
    
    # Calculate total number of water molecules (for CONECT records)
    n_molecules = (atom_index - 1) // 3

    # --- Write the PDB file ---
    with open(pdb_filename, 'w') as f:
        # 1. Write CRYST1 record
        f.write(CRYST1_FMT % (cell_a, cell_b, cell_c, 90.00, 90.00, 90.00))
        
        # 2. Write the ATOM records
        f.writelines(pdb_lines)
        
        # 3. Write CONECT records (O connected to its H1 and H2)
        for i in range(1, n_molecules + 1):
            o_atom_idx = (i - 1) * 3 + 1
            h1_atom_idx = o_atom_idx + 1
            h2_atom_idx = o_atom_idx + 2
            
            # CONECT record links the central O atom to its two H atoms
            conect_line = f"CONECT{o_atom_idx:5d}{h1_atom_idx:5d}{h2_atom_idx:5d}\n"
            f.write(conect_line)
            
        f.write("END\n")
        
    print(f"✅ Successfully converted {xyz_filename} to {pdb_filename}")
    print(f"Total molecules written: {n_molecules}")
    print(f"Total atoms written: {atom_index - 1}")

# --- Execution Parameters ---

# 1. Input/Output file paths
input_xyz = "ice_supercell_4x4x4.xyz"
output_pdb = "ice_supercell_4x4x4.pdb"

# 2. Cell Dimensions (Based on GenIce log: a=7.848134, b=7.377351, c=9.065738)
# We multiply by 4 because you requested a 4x4x4 supercell.
A = 7.84813412606925 * 4 
B = 7.37735062301457 * 4 
C = 9.06573834219084 * 4 

# 3. Total Atoms (Based on GenIce log: 1024 water molecules * 3 atoms/mol)
N_ATOMS = 1024 * 3 

print(f"Using estimated final box dimensions: A={A:.3f}, B={B:.3f}, C={C:.3f} Å")
print(f"Expecting {N_ATOMS} total atoms ({N_ATOMS/3} water molecules).")

# Run the converter
xyz_to_pdb_water(input_xyz, output_pdb, A, B, C, N_ATOMS)