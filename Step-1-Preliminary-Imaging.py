import iolite.ui as ui
import numpy as np
import re

# ============================================================================
# STEP 1: PRELIMINARY NORMALIZATION AND RATIO DIAGNOSTICS
# ============================================================================
# This script performs preliminary oxide-sum normalization using standard
# oxide factors to generate internally consistent element/Si ratio 
# diagnostics for downstream phase segmentation.
#
# Reference: Table A4 (Default Silicate/Matrix factors) in Supplementary Text S1
# ============================================================================

# ================= 1. TARGET ISOTOPE CONFIGURATION =================
# Isotope selections from Table "01-Adopted-isotopes.csv"
TARGET_ISOTOPES = {
    'Li': '7',   'B': '11',   'Na': '23',  'Mg': '24',  'Al': '27', 
    'Si': '28',  'P': '31',   'S': '32',   'Cl': '35',  'K': '39', 
    'Ca': '44',  'Sc': '45',  'Ti': '47',  'V': '51',   'Cr': '52', 
    'Mn': '55',  'Fe': '56',  'Co': '59',  'Ni': '62',  'Cu': '63', 
    'Zn': '66',  'Ga': '71',  'As': '75',  'Rb': '85',  'Sr': '88', 
    'Y': '89',   'Zr': '90',  'Nb': '93',  'Mo': '98',  'Cs': '133', 
    'Ba': '138', 'La': '139', 'Ce': '140', 'Pr': '141', 'Nd': '146', 
    'Sm': '147', 'Eu': '151', 'Gd': '158', 'Tb': '159', 'Dy': '163', 
    'Ho': '165', 'Er': '166', 'Tm': '169', 'Yb': '172', 'Lu': '175', 
    'Hf': '180', 'Pb': '208', 'Th': '232', 'U': '238'
}

# ================= 2. STOICHIOMETRIC FACTORS =================

# Elements excluded from oxide-sum normalization
# S and Cl are excluded because their mass contributions are already
# accounted for in compound stoichiometry (e.g., FeS2, SrSO4)
EXCLUDE_FROM_SUM = ['S', 'O', 'Ar', 'He', 'Total', 'Sum', 'Cl'] 

# Standard Oxide Factors (Default Silicate/Matrix Phase)
# These factors are applied during preliminary normalization
# Derived from Table A4, Column "Silicate/Matrix (Default Phase 5)"
FAC_DEFAULT = {
    'Si': 2.139,  # SiO2
    'Al': 1.889,  # Al2O3
    'K': 1.205,   # K2O
    'Na': 1.348,  # Na2O
    'Ti': 1.668,  # TiO2
    'Mg': 1.658,  # MgO
    'Ca': 1.399,  # CaO (silicate form, not carbonate)
    'Fe': 1.430,  # Fe2O3 (assumed Fe(III) in silicate)
    'Mn': 1.291,  # MnO (typical Mn(II) in silicate)
    'Pb': 1.077,  # PbO
    'P': 2.291,   # P2O5
    'Sr': 1.183,  # SrO
    'Ba': 1.117,  # BaO
    'Zn': 1.245,  # ZnO
    'Zr': 1.351,  # ZrO2
    'Cr': 1.462,  # Cr2O3
    'V': 1.785,   # V2O5
    'Ni': 1.273,  # NiO
    'Co': 1.362,  # Co3O4
    'Cu': 1.252,  # CuO
    'Li': 2.153,  # Li2O
    'Rb': 1.094,  # Rb2O
    'Cs': 1.060,  # Cs2O
    'Sc': 1.534,  # Sc2O3
    'Y': 1.270,   # Y2O3
    'La': 1.173,  # La2O3
    'Ce': 1.171,  # CeO2
    'Pr': 1.170,  # Pr2O3
    'Nd': 1.166,  # Nd2O3
    'Sm': 1.160,  # Sm2O3
    'Eu': 1.158,  # Eu2O3
    'Gd': 1.153,  # Gd2O3
    'Tb': 1.151,  # Tb2O3
    'Dy': 1.148,  # Dy2O3
    'Ho': 1.146,  # Ho2O3
    'Er': 1.143,  # Er2O3
    'Tm': 1.142,  # Tm2O3
    'Yb': 1.139,  # Yb2O3
    'Lu': 1.137,  # Lu2O3
    'Hf': 1.179,  # HfO2
    'Th': 1.138,  # ThO2
    'U': 1.134,   # UO2
    'Ga': 1.344,  # Ga2O3
    'As': 1.534,  # As2O3
    'Mo': 1.500   # MoO3
}

# ================= 3. MAIN PROCESSING FUNCTION =================

def run_step1_preliminary():
    """
    Perform preliminary normalization using standard oxide-sum approach.

    This function:
    1. Loads elemental concentration channels from Iolite workspace
    2. Applies standard oxide stoichiometric factors
    3. Normalizes to 100 wt.% (1,000,000 ppm) oxide sum
    4. Outputs preliminary concentrations for ratio diagnostics

    Output channels: [Element][Mass]_Prelim_ppm
    """
    print("--- STEP 1: Preliminary Normalization (Standard Oxide Sum) ---")

    # 1. Load data channels
    all_channels = data.timeSeriesList()
    data_map = {}
    channel_map = {}

    print("Loading elemental channels...")
    for ch in all_channels:
        try: 
            name = ch.name()
        except: 
            name = ch.name

        # Only process channels ending with '_ppm' (raw concentrations)
        if not name.endswith('_ppm'): 
            continue

        # Parse element and mass number (e.g., 'Ca44_ppm' -> Element='Ca', Mass='44')
        match = re.match(r"([A-Za-z]+)(\d+)", name)
        if not match: 
            continue
        el = match.group(1)
        mass = match.group(2)

        # Filter by target isotopes
        if el in TARGET_ISOTOPES and TARGET_ISOTOPES[el] == mass:
            d = np.nan_to_num(ch.data()).astype(np.float64)
            data_map[el] = d
            channel_map[el] = ch

    if not data_map:
        print("ERROR: No valid elemental channels found. Check channel naming.")
        return

    print(f"Loaded {len(data_map)} elements: {list(data_map.keys())}")

    # 2. Add Si floor for numerical stability
    # A small constant offset is added to Si to prevent division by zero
    # in Si-poor phases (e.g., sulfides, carbonates) while remaining
    # negligible for Si-rich silicate matrices
    SI_FLOOR = 1.0  # ppm
    if 'Si' in data_map:
        data_map['Si'] += SI_FLOOR
        print(f"Applied Si floor: +{SI_FLOOR} ppm for numerical stability")

    # 3. Calculate preliminary oxide masses
    print("\nCalculating oxide masses...")
    oxide_mass_map = {}

    for el, conc in data_map.items():
        if el in EXCLUDE_FROM_SUM:
            continue  # Skip S, Cl, etc.

        fac = FAC_DEFAULT.get(el, 1.0)  # Default factor = 1.0 if not specified
        oxide_mass_map[el] = conc * fac
        print(f"  {el}: Factor = {fac:.3f}")

    # 4. Calculate total oxide mass per pixel
    total_mass = np.zeros_like(list(data_map.values())[0])
    for el, mass in oxide_mass_map.items():
        total_mass += mass

    # Prevent division by zero
    total_mass = np.where(total_mass > 0, total_mass, 1.0)

    print(f"\nTotal oxide mass range: {total_mass.min():.1f} - {total_mass.max():.1f} ppm")

    # 5. Normalize to 100 wt.% (1,000,000 ppm)
    TARGET_SUM = 1000000.0  # ppm (equivalent to 100 wt.%)
    norm_factor = TARGET_SUM / total_mass

    print(f"Applying normalization factor (range: {norm_factor.min():.3f} - {norm_factor.max():.3f})...")

    # 6. Create output channels
    print("\nCreating preliminary concentration channels...")
    t = channel_map[list(channel_map.keys())[0]].time()

    for el, raw_conc in data_map.items():
        prelim_conc = raw_conc * norm_factor

        # Output channel naming: [Element][Mass]_Prelim_ppm
        iso = TARGET_ISOTOPES[el]
        out_name = f"{el}{iso}_Prelim_ppm"

        data.createTimeSeries(out_name, data.Output, t, prelim_conc)
        print(f"  Created: {out_name}")

    print("\n>>> STEP 1 COMPLETE: Preliminary normalization finished <<<")
    print("    Next: Run Step 2 for phase segmentation using element/Si ratios")

# ================= 4. EXECUTION =================
if __name__ == "__main__":
    run_step1_preliminary()
