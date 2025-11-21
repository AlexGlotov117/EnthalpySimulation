import openmm.unit as unit

N_MOLECULES = 3072                     # Number of molecules in the simulation box
TARGET_DENSITY = 0.998 * unit.gram/unit.milliliter
MOLAR_MASS = 18.01528 * unit.gram/unit.mole

total_mass = N_MOLECULES * (MOLAR_MASS / unit.AVOGADRO_CONSTANT_NA)
required_volume = total_mass / TARGET_DENSITY
box_length = required_volume**(1/3) # Cube root for cubic box

L = box_length.value_in_unit(unit.angstrom)

a_len = f"{L:9.3f}"
b_len = f"{L:9.3f}"
c_len = f"{L:9.3f}"
alpha = f"{90.00:7.2f}"
beta  = f"{90.00:7.2f}"
gamma = f"{90.00:7.2f}"

# Construct the CRYST1 record string
# For liquid systems, the space group is often 'P 1' and Z=1
print(f"CRYST1{a_len}{b_len}{c_len}{alpha}{beta}{gamma} P 1      1")

packmol_input = f"""
# Packmol control file

# Define the output file name and format
output liquid_water_packed.pdb
filetype pdb

# Set tolerance for atomic overlap (0.0001 Angstroms)
tolerance 2.0

# Define the molecule to be packed
structure water.pdb
  number {N_MOLECULES}
  # Place the molecules randomly within the specified cubic box
  inside box 0. 0. 0. {box_length.value_in_unit(unit.angstrom):.2f} {box_length.value_in_unit(unit.angstrom):.2f} {box_length.value_in_unit(unit.angstrom):.2f}
end structure
"""

with open('packmol.in', 'w') as f:
    f.write(packmol_input)

print("Created 'packmol.in'")