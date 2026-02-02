import iolite.ui as ui
import numpy as np
import re

# ============================================================================
# STEP 3: PHASE-SPECIFIC QUANTIFICATION (FINAL NORMALIZATION)
# ============================================================================
# This script applies phase-dependent stoichiometric constraints to 
# calculate final element concentrations using adaptive normalization.
#
# Each pixel's normalization factor is determined by its assigned mineral
# phase (from Step 2), ensuring that compound masses sum to 100 wt.%.
#
# Reference: Table A4 in Supplementary Text S1
# ============================================================================

# ================= 1. TARGET ISOTOPE CONFIGURATION =================
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

# ================= 2. STOICHIOMETRIC FACTORS BY PHASE =================
# ⚠️ CRITICAL: All factors MUST match Table A4 exactly

# Elements excluded from mass summation (S already counted in sulfides/sulfates)
EXCLUDE_FROM_SUM = ['S', 'O', 'Ar', 'He', 'Total', 'Sum', 'Cl']

# --- Phase 1: Pyrite/Sulfide ---
# Table A4, Section "Pyrite/Sulfide (Phase 1)"
FAC_PYRITE = {
    'Fe': 2.148,  # FeS2 (disulfide backbone)
    'Pb': 1.155,  # PbS (galena, monosulfide)
    'Zn': 1.490,  # ZnS (sphalerite, monosulfide)
    'Cu': 1.504,  # CuS (covellite approximation)
    'Co': 2.088,  # CoS2 (disulfide substitution)
    'Ni': 2.091,  # NiS2 (vaesite, disulfide substitution)
    'As': 1.428   # AsS (arsenopyrite/realgar approximation)
}

# --- Phase 2: Sulfate ---
# Table A4, Section "Sulfate (Phase 2)"
FAC_SULFATE = {
    'Sr': 2.096,  # SrSO4 (celestite, primary sulfate)
    'Ba': 1.699,  # BaSO4 (barite, primary sulfate)
    'Pb': 1.463,  # PbSO4 (anglesite)
    'Ca': 3.395   # CaSO4 (anhydrite/gypsum)
}

# --- Phase 3: Fe-Mn Oxide ---
# Table A4, Section "Fe-Mn oxide (Phase 3)"
FAC_OXIDE = {
    'Mn': 1.582,  # MnO2 (pyrolusite, Mn(IV) anhydrous)
    'Fe': 1.430,  # Fe2O3 (hematite, Fe(III) anhydrous)
    'Co': 1.362,  # Co3O4 (oxide/adsorbed phase)
    'Ni': 1.273   # NiO (oxide/adsorbed phase)
}

# --- Phase 4: Carbonate ---
# Table A4, Section "Carbonate (Phase 4)"
FAC_CARBONATE = {
    'Ca': 2.497,  # CaCO3 (calcite/aragonite, primary carbonate)
    'Mg': 3.468,  # MgCO3 (magnesite/dolomite backbone)
    'Fe': 2.075,  # FeCO3 (siderite, solid solution component)
    'Mn': 2.093,  # MnCO3 (rhodochrosite, solid solution component)
    'Sr': 1.685,  # SrCO3 (strontianite, lattice substitution)
    'Ba': 1.437,  # BaCO3 (witherite, solid solution component)
    'Zn': 1.920,  # ZnCO3 (smithsonite, solid solution component)
    'Na': 2.305   # Na2CO3 (assumed carbonate substitution/defect)
}

# --- Phase 5: Silicate/Matrix (Default Fallback) ---
# Table A4, Section "Silicate/Matrix (Default Phase 5)"
FAC_SILICATE = {
    'Si': 2.139,  # SiO2 (default backbone)
    'Al': 1.889,  # Al2O3
    'K': 1.205,   # K2O
    'Na': 1.348,  # Na2O
    'Ti': 1.668,  # TiO2
    'Mg': 1.658,  # MgO
    'Ca': 1.399,  # CaO (silicate form, NOT carbonate)
    'Fe': 1.430,  # Fe2O3 (assumed Fe(III) in silicate)
    'Mn': 1.291,  # MnO (typical Mn(II) in silicate)
    'Pb': 1.077,  # PbO
    'P': 2.291,   # P2O5 (apatite/phosphate in silicate matrix)
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

# Phase classification codes (must match Step 2 output)
PHASE_CODES = {
    'Pyrite':     1,
    'Sulfate':    2,
    'FeMn_Oxide': 3,
    'Carbonate':  4,
    'Silicate':   5
}

# ================= 3. MAIN QUANTIFICATION FUNCTION =================

def run_step3_final_calculation():
    """
    Apply phase-specific stoichiometric normalization for final quantification.

    Workflow:
    1. Load phase classification map from Step 2 ('Adapt_Phase_Class')
    2. Load preliminary concentrations from Step 1 ('_Prelim_ppm' channels)
    3. For each pixel, select stoichiometric factors based on assigned phase
    4. Calculate total compound mass (e.g., FeS2, CaCO3, SiO2)
    5. Normalize to 100 wt.% (1,000,000 ppm) compound sum
    6. Output final concentrations as '[Element][Mass]_Final_ppm'

    Reference: Section 2 "Calculation Logic" in Supplementary Text S1
    """
    print("--- STEP 3: Phase-Specific Quantification (Final Step) ---")

    # 1. Load phase classification map
    print("\nLoading phase classification map...")
    try:
        phase_channel = data.timeSeries("Adapt_Phase_Class")
        phase_map = np.nan_to_num(phase_channel.data()).astype(int)
        print("  ✓ Phase map loaded: 'Adapt_Phase_Class'")
    except:
        print("  ERROR: Phase map 'Adapt_Phase_Class' not found.")
        print("         Please run Step 2 first to generate phase classification.")
        return

    # 2. Load preliminary concentration channels
    all_channels = data.timeSeriesList()
    data_map = {}
    channel_map = {}

    print("\nLoading preliminary concentration channels...")
    for ch in all_channels:
        try: 
            name = ch.name()
        except: 
            name = ch.name

        # Load '_Prelim_ppm' or fallback to '_ppm'
        if not (name.endswith('_Prelim_ppm') or name.endswith('_ppm')): 
            continue

        match = re.match(r"([A-Za-z]+)(\d+)", name)
        if not match: 
            continue
        el = match.group(1)
        mass = match.group(2)

        # Prioritize '_Prelim_ppm'
        if el in data_map and name.endswith('_ppm') and not name.endswith('_Prelim_ppm'):
            continue

        if el in TARGET_ISOTOPES and TARGET_ISOTOPES[el] == mass:
            d = np.nan_to_num(ch.data()).astype(np.float64)
            data_map[el] = d
            channel_map[el] = ch

    if not data_map:
        print("  ERROR: No preliminary concentration channels found.")
        return

    print(f"  ✓ Loaded {len(data_map)} elements: {list(data_map.keys())}")

    # 3. Calculate phase-adaptive normalization factors
    print("\nCalculating phase-specific normalization factors...")

    n_pixels = len(phase_map)
    total_mass_per_pixel = np.zeros(n_pixels, dtype=np.float64)

    # Map phase codes to factor dictionaries
    FACTOR_MAP = {
        PHASE_CODES['Pyrite']:     FAC_PYRITE,
        PHASE_CODES['Sulfate']:    FAC_SULFATE,
        PHASE_CODES['FeMn_Oxide']: FAC_OXIDE,
        PHASE_CODES['Carbonate']:  FAC_CARBONATE,
        PHASE_CODES['Silicate']:   FAC_SILICATE
    }

    # Calculate total compound mass for each pixel
    for el, conc_array in data_map.items():
        if el in EXCLUDE_FROM_SUM:
            continue  # Skip S, Cl, etc.

        # For each phase, get the corresponding stoichiometric factor
        for phase_code, fac_dict in FACTOR_MAP.items():
            mask = (phase_map == phase_code)
            factor = fac_dict.get(el, FAC_SILICATE.get(el, 1.0))  # Fallback to silicate factor
            total_mass_per_pixel[mask] += conc_array[mask] * factor

    # Prevent division by zero
    total_mass_per_pixel = np.where(total_mass_per_pixel > 0, total_mass_per_pixel, 1.0)

    print(f"  Total compound mass range: {total_mass_per_pixel.min():.1f} - {total_mass_per_pixel.max():.1f} ppm")

    # 4. Calculate normalization factor
    TARGET_SUM = 1000000.0  # ppm (100 wt.%)
    norm_factor = TARGET_SUM / total_mass_per_pixel

    print(f"  Normalization factor range: {norm_factor.min():.3f} - {norm_factor.max():.3f}")

    # 5. Apply normalization and create final output channels
    print("\nCreating final quantitative channels...")
    t = channel_map[list(channel_map.keys())[0]].time()

    for el, prelim_conc in data_map.items():
        final_conc = prelim_conc * norm_factor

        # Output channel naming: [Element][Mass]_Final_ppm
        iso = TARGET_ISOTOPES[el]
        out_name = f"{el}{iso}_Final_ppm"

        data.createTimeSeries(out_name, data.Output, t, final_conc)
        print(f"  Created: {out_name}")

    # 6. Summary statistics by phase
    print("\n" + "="*70)
    print("FINAL QUANTIFICATION SUMMARY")
    print("="*70)

    for phase_name, phase_code in PHASE_CODES.items():
        mask = (phase_map == phase_code)
        n_pix = np.sum(mask)
        if n_pix > 0:
            avg_norm = norm_factor[mask].mean()
            print(f"  {phase_name:15s} (Phase {phase_code}): {n_pix:6d} pixels, avg norm factor = {avg_norm:.3f}")

    print("\n>>> STEP 3 COMPLETE: Phase-specific quantification finished <<<")
    print("    Output channels: [Element][Mass]_Final_ppm")
    print("    All concentrations normalized to 100 wt.% compound sum")

# ================= 4. EXECUTION =================
if __name__ == "__main__":
    run_step3_final_calculation()
