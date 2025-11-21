import openmm.unit as unit
import pandas as pd

# --- 1. Define Constants and Target Pressure ---
LOG_FILE = 'output_solid.log'
PROD_START_STEP = 25000 
N_MOLECULES = 3072
N_A = 6.02214076e23 / unit.mole
N_MOLES= N_MOLECULES / N_A
pressure = 1.0*unit.bar # Target pressure (e.g., standard pressure for liquid)

# --- 2. Load and Average Data ---
data = pd.read_csv(LOG_FILE, sep=',', header=0) # Adjust skiprows as needed
prod_data = data[data['#"Step"'] >= PROD_START_STEP]

avg_U = prod_data['Total Energy (kJ/mole)'].mean() * unit.kilojoule_per_mole / N_MOLECULES
avg_V = prod_data['Box Volume (nm^3)'].mean() * (unit.nanometer)**3

# --- 3. Calculate Molar Enthalpy (H_molar) ---

PV_term = (pressure * avg_V) / N_MOLES

H_molar_kJ_per_mole = avg_U + PV_term

print(f"Average Internal Energy <U>: {avg_U}")
print(f"Target Pressure <P>: {pressure} (used for calculation)")
print(f"Average Volume <V>: {avg_V}")
print(f"PV Term: {PV_term.value_in_unit(unit.kilojoule_per_mole)*unit.kilojoule_per_mole}")
print(f"Total System Enthalpy (H): {H_molar_kJ_per_mole}")