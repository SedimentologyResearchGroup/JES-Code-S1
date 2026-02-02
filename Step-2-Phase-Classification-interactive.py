import iolite.ui as ui
import numpy as np
import re

# ============================================================================
# STEP 2: PHASE SEGMENTATION (INTERACTIVE VERSION)
# ============================================================================
# This script performs pixel-by-pixel phase classification using element/Si
# ratio diagnostics derived from Step 1 preliminary concentrations.
#
# Classification follows a hierarchical decision tree to minimize 
# misclassification (Priority Order: Pyrite → Sulfate → Fe-Mn Oxide → 
# Carbonate → Silicate/Matrix fallback).
#
# ⚠️ USER ADJUSTABLE: Threshold values below are based on Table A3 in
# Supplementary Text S1. You may adjust these thresholds based on your
# specific sample characteristics and mineralogy.
#
# Reference: Table A3 in Supplementary Text S1
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

# ================= 2. PHASE SEGMENTATION THRESHOLDS =================
# ⚠️⚠️⚠️ USER ADJUSTABLE SECTION ⚠️⚠️⚠️
#
# These thresholds are based on Table A3 in Supplementary Text S1 and
# represent optimized values for the reference dataset. You should adjust
# these values based on:
#   - Your specific mineral assemblage
#   - Expected phase abundances
#   - Visual inspection of preliminary element/Si ratio distributions
#
# To adjust thresholds:
#   1. Run Step 1 to generate preliminary concentrations
#   2. Visualize element/Si ratio distributions in Iolite
#   3. Identify appropriate threshold values for phase separation
#   4. Modify the R_TH dictionary below
#   5. Re-run this script and verify results with Step 2 output
#
# All thresholds are based on element/Si ratios calculated from 
# preliminary concentrations (Step 1 output)

# Silicon floor for ratio calculation (must match Step 1)
SI_FLOOR = 1.0  # ppm

# Threshold dictionary (Ratio Name: Threshold Value)
# Reference: Table A3, Column "Segmentation thresholds"
R_TH = {
    # === 1. Pyrite (Priority: High - 1st) ===
    'Fe_Si_Pyrite': 0.015,    # Fe/Si > 0.015 (Table A3)
                              # Adjust if your pyrite has different Fe content

    'S_Si_Pyrite':  0.001,    # S/Si > 0.001 (Table A3)
                              # Lower threshold for low-S pyrite variants

    # === 2. Sulfate (Priority: Medium - 2nd) ===
    'Sr_Si_Sulf':   0.5,      # Sr/Si > 0.5 (Table A3)
                              # Adjust for celestite (SrSO4) identification

    'Ba_Si_Sulf':   0.05,     # Ba/Si > 0.05 (Table A3)
                              # Adjust for barite (BaSO4) identification

    # === 3. Fe-Mn Oxide (Priority: Medium - 3rd) ===
    'Mn_Si_Oxide':  0.01,     # Mn/Si > 0.01 (Table A3)
                              # Adjust for pyrolusite/ferric oxide mixtures

    # === 4. Carbonate (Priority: Low - 4th) ===
    'Ca_Si_Carb':   1.5,      # Ca/Si > 1.5 (Table A3)
                              # Adjust based on carbonate vs. silicate ratio

    'Al_Ca_Limit':  0.005     # Al/Ca < 0.005 (Table A3)
                              # This filter excludes Ca-bearing aluminosilicates
                              # (e.g., feldspars) from carbonate classification
                              # Increase if your carbonates have high Al
}

# ⚠️ END OF USER ADJUSTABLE SECTION

# Phase classification codes (for output channel "Adapt_Phase_Class")
PHASE_CODES = {
    'Pyrite':       1,  # FeS2-dominated
    'Sulfate':      2,  # SrSO4/BaSO4-dominated
    'FeMn_Oxide':   3,  # MnO2/Fe2O3-dominated
    'Carbonate':    4,  # CaCO3-dominated
    'Silicate':     5   # Default silicate/oxide matrix (fallback)
}

# ================= 3. HELPER FUNCTION =================

def get_d(data_map, el):
    """
    Safely retrieve element concentration array.
    Returns zero array if element not found.

    Args:
        data_map (dict): Dictionary of element concentration arrays
        el (str): Element symbol (e.g., 'Ca', 'Fe', 'Si')

    Returns:
        np.ndarray: Concentration array (ppm)
    """
    if el in data_map:
        return data_map[el]
    else:
        # Return zero array with same shape as first element in data_map
        template = list(data_map.values())[0]
        return np.zeros_like(template)

# ================= 4. MAIN CLASSIFICATION FUNCTION =================

def run_step2_classification():
    """
    Perform hierarchical phase segmentation using element/Si ratios.

    Classification logic (in order of priority):
    1. Pyrite:      High Fe/Si AND high S/Si
    2. Sulfate:     High Sr/Si OR high Ba/Si (celestite/barite)
    3. Fe-Mn Oxide: High Mn/Si (pyrolusite/ferric oxide)
    4. Carbonate:   High Ca/Si AND low Al/Ca (exclude aluminosilicates)
    5. Silicate:    All remaining pixels (default fallback)

    Output channel: Adapt_Phase_Class (integer codes 1-5)
    """
    print("--- STEP 2: Phase Classification (Interactive Mode) ---")

    # 1. Load preliminary concentration channels from Step 1
    all_channels = data.timeSeriesList()
    data_map = {}
    channel_map = {}

    print("Loading preliminary channels from Step 1...")
    for ch in all_channels:
        try: 
            name = ch.name()
        except: 
            name = ch.name

        # Load channels ending with '_Prelim_ppm' (Step 1 output)
        # If not found, fallback to '_ppm' (raw data)
        if not (name.endswith('_Prelim_ppm') or name.endswith('_ppm')): 
            continue

        # Parse element and mass
        match = re.match(r"([A-Za-z]+)(\d+)", name)
        if not match: 
            continue
        el = match.group(1)
        mass = match.group(2)

        # Prioritize '_Prelim_ppm' over raw '_ppm'
        if el in data_map and name.endswith('_ppm') and not name.endswith('_Prelim_ppm'):
            continue  # Skip raw data if preliminary data already loaded

        if el in TARGET_ISOTOPES and TARGET_ISOTOPES[el] == mass:
            d = np.nan_to_num(ch.data()).astype(np.float64)
            data_map[el] = d
            channel_map[el] = ch

    if not data_map:
        print("ERROR: No preliminary concentration channels found.")
        print("       Please run Step 1 first to generate '_Prelim_ppm' channels.")
        return

    print(f"Loaded {len(data_map)} elements: {list(data_map.keys())}")

    # 2. Apply Si floor (must match Step 1 to ensure ratio consistency)
    if 'Si' in data_map:
        data_map['Si'] += SI_FLOOR

    # 3. Extract required element arrays
    Fe = get_d(data_map, 'Fe')
    S  = get_d(data_map, 'S')
    Ca = get_d(data_map, 'Ca')
    Al = get_d(data_map, 'Al')
    Sr = get_d(data_map, 'Sr')
    Ba = get_d(data_map, 'Ba')
    Mn = get_d(data_map, 'Mn')
    Si = get_d(data_map, 'Si')

    # 4. Calculate diagnostic element/Si ratios
    print("\nCalculating diagnostic ratios...")

    # Prevent division by zero
    Si_safe = np.where(Si > 0, Si, 1.0)

    Fe_Si = Fe / Si_safe
    S_Si  = S  / Si_safe
    Ca_Si = Ca / Si_safe
    Sr_Si = Sr / Si_safe
    Ba_Si = Ba / Si_safe
    Mn_Si = Mn / Si_safe

    # Al/Ca ratio (for carbonate filter)
    Ca_safe = np.where(Ca > 0, Ca, 1.0)
    Al_Ca = Al / Ca_safe

    print(f"  Fe/Si range: {Fe_Si.min():.4f} - {Fe_Si.max():.4f}")
    print(f"  S/Si range:  {S_Si.min():.5f} - {S_Si.max():.4f}")
    print(f"  Ca/Si range: {Ca_Si.min():.4f} - {Ca_Si.max():.4f}")
    print(f"  Sr/Si range: {Sr_Si.min():.4f} - {Sr_Si.max():.4f}")
    print(f"  Ba/Si range: {Ba_Si.min():.4f} - {Ba_Si.max():.4f}")
    print(f"  Mn/Si range: {Mn_Si.min():.4f} - {Mn_Si.max():.4f}")

    # 5. Initialize phase classification map
    n_pixels = len(Fe)
    phase_class_map = np.full(n_pixels, PHASE_CODES['Silicate'], dtype=np.int32)

    # 6. Hierarchical phase assignment
    print("\nApplying hierarchical phase segmentation...")
    print("(Using thresholds from Table A3 - see R_TH dictionary for values)")

    # === Priority 1: Pyrite (FeS2) ===
    # Both Fe/Si AND S/Si must exceed thresholds
    mask_pyrite = (Fe_Si > R_TH['Fe_Si_Pyrite']) & (S_Si > R_TH['S_Si_Pyrite'])
    phase_class_map[mask_pyrite] = PHASE_CODES['Pyrite']
    n_pyrite = np.sum(mask_pyrite)
    print(f"  [1] Pyrite pixels: {n_pyrite} ({100*n_pyrite/n_pixels:.2f}%)")
    print(f"      Criteria: Fe/Si > {R_TH['Fe_Si_Pyrite']} AND S/Si > {R_TH['S_Si_Pyrite']}")

    # === Priority 2: Sulfate (SrSO4, BaSO4) ===
    # Sr/Si OR Ba/Si exceeds threshold (exclude pixels already classified as pyrite)
    mask_sulfate = (~mask_pyrite) & ((Sr_Si > R_TH['Sr_Si_Sulf']) | (Ba_Si > R_TH['Ba_Si_Sulf']))
    phase_class_map[mask_sulfate] = PHASE_CODES['Sulfate']
    n_sulfate = np.sum(mask_sulfate)
    print(f"  [2] Sulfate pixels: {n_sulfate} ({100*n_sulfate/n_pixels:.2f}%)")
    print(f"      Criteria: Sr/Si > {R_TH['Sr_Si_Sulf']} OR Ba/Si > {R_TH['Ba_Si_Sulf']}")

    # === Priority 3: Fe-Mn Oxide (MnO2, Fe2O3) ===
    mask_oxide = (~mask_pyrite) & (~mask_sulfate) & (Mn_Si > R_TH['Mn_Si_Oxide'])
    phase_class_map[mask_oxide] = PHASE_CODES['FeMn_Oxide']
    n_oxide = np.sum(mask_oxide)
    print(f"  [3] Fe-Mn Oxide pixels: {n_oxide} ({100*n_oxide/n_pixels:.2f}%)")
    print(f"      Criteria: Mn/Si > {R_TH['Mn_Si_Oxide']}")

    # === Priority 4: Carbonate (CaCO3) ===
    # High Ca/Si AND low Al/Ca (to exclude Ca-bearing aluminosilicates like feldspars)
    mask_carbonate = (~mask_pyrite) & (~mask_sulfate) & (~mask_oxide) & \
                     (Ca_Si > R_TH['Ca_Si_Carb']) & (Al_Ca < R_TH['Al_Ca_Limit'])
    phase_class_map[mask_carbonate] = PHASE_CODES['Carbonate']
    n_carbonate = np.sum(mask_carbonate)
    print(f"  [4] Carbonate pixels: {n_carbonate} ({100*n_carbonate/n_pixels:.2f}%)")
    print(f"      Criteria: Ca/Si > {R_TH['Ca_Si_Carb']} AND Al/Ca < {R_TH['Al_Ca_Limit']}")

    # === Priority 5: Silicate/Matrix (Default fallback) ===
    mask_silicate = ~(mask_pyrite | mask_sulfate | mask_oxide | mask_carbonate)
    n_silicate = np.sum(mask_silicate)
    print(f"  [5] Silicate/Matrix pixels: {n_silicate} ({100*n_silicate/n_pixels:.2f}%)")
    print(f"      (Remaining pixels assigned to default phase)")

    # 7. Create output channel
    print("\nCreating phase classification channel...")
    t = channel_map[list(channel_map.keys())[0]].time()

    output_channel_name = "Adapt_Phase_Class"
    data.createTimeSeries(output_channel_name, data.Output, t, phase_class_map)

    print(f"  Created: {output_channel_name}")
    print("  Phase codes: 1=Pyrite, 2=Sulfate, 3=Fe-Mn Oxide, 4=Carbonate, 5=Silicate")

    print("\n>>> STEP 2 COMPLETE: Phase segmentation finished <<<")
    print("    Next: Run Step 3 for phase-specific quantification")

    # 8. Interactive user confirmation
    response = ui.questionBox(
        "Phase Classification Complete",
        f"Phase segmentation results:\n\n" +
        f"  Pyrite:       {n_pyrite:6d} pixels ({100*n_pyrite/n_pixels:5.2f}%)\n" +
        f"  Sulfate:      {n_sulfate:6d} pixels ({100*n_sulfate/n_pixels:5.2f}%)\n" +
        f"  Fe-Mn Oxide:  {n_oxide:6d} pixels ({100*n_oxide/n_pixels:5.2f}%)\n" +
        f"  Carbonate:    {n_carbonate:6d} pixels ({100*n_carbonate/n_pixels:5.2f}%)\n" +
        f"  Silicate:     {n_silicate:6d} pixels ({100*n_silicate/n_pixels:5.2f}%)\n\n" +
        f"Channel '{output_channel_name}' created.\n\n" +
        f"⚠️ If classification results are unsatisfactory, you can:\n" +
        f"   1. Adjust threshold values in R_TH dictionary (see Table A3)\n" +
        f"   2. Re-run this script\n" +
        f"   3. Visualize phase map to verify results\n\n" +
        f"Proceed to Step 3?",
        ui.Yes | ui.No
    )

    if response == ui.Yes:
        print("User confirmed: Proceeding to Step 3...")
    else:
        print("User chose to review results. Run Step 3 manually when ready.")

# ================= 5. EXECUTION =================
if __name__ == "__main__":
    run_step2_classification()
