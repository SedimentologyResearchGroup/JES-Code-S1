import iolite.ui as ui
import numpy as np
import re

"""
Supplementary Software Code 1
-----------------------------
Title: Adaptive Phase-Specific Normalization Script for LA-ICP-TOF-MS
Description: This script performs pixel-by-pixel mineral phase identification 
             based on element/Si ratios and applies phase-specific 
             stoichiometric factors for quantitative normalization.
Platform: Iolite 4 (Python)
"""

# ==========================================
# 1. USER CONFIGURATION (THRESHOLDS)
# ==========================================
# Diagnostic thresholds determined from preliminary ratio maps.
USER_TH = {
    # Fallback strategy for unclassified pixels
    'Default_Phase': 'Silicate/Oxide Matrix (Standard Oxide Norm)',

    # [Phase 1] Pyrite (Highest Priority)
    # Criteria: High Fe/Si and presence of S
    'Fe_Si_Min': 0.015,
    'S_Si_Min':  0.001,
    
    # [Phase 2] Sulfate (Second Priority)
    # Criteria: High Sr (Celestine) or Ba (Barite) relative to Si
    'Sr_Si_Min': 0.5,
    'Ba_Si_Min': 0.05,
    
    # [Phase 5] Fe-Mn Oxide
    # Criteria: High Mn relative to Si
    'Mn_Si_Min': 0.01,

    # [Phase 4] Carbonate (Lowest Priority among minerals)
    # Criteria: High Ca/Si, low Al/Ca (to exclude Ca-rich silicates)
    'Ca_Si_Min': 1.5,
    'Al_Ca_Max': 0.005
}

# ==========================================
# 2. SYSTEM PARAMETERS & FACTORS
# ==========================================

# Target Isotopes (Element + Mass) validation list
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

# Elements excluded from Total Mass calculation
# S is excluded as it is accounted for in sulfide factors or negligible in oxides
EXCLUDE_FROM_SUM = ['S', 'O', 'Ar', 'He', 'Total', 'Sum', 'Cl'] 

# Stoichiometric Factors (Mass conversion factors)
# Default: Standard Oxides (for Silicate/Matrix)
FAC_DEFAULT = {
    'Si': 2.139, 'Al': 1.889, 'K': 1.205, 'Na': 1.348, 'Ti': 1.668,
    'Mg': 1.658, 'Ca': 1.399, 'Fe': 1.430, 'Mn': 1.291, 'Pb': 1.077,
    'P': 2.291,  'Sr': 1.183, 'Ba': 1.117, 'Zn': 1.245, 'Zr': 1.351,
    'Cr': 1.462, 'V': 1.785,  'Ni': 1.273, 'Co': 1.362, 'Cu': 1.252,
    'Li': 2.153, 'Rb': 1.094, 'Cs': 1.060, 'Sc': 1.534, 'Y': 1.270,
    'La': 1.173, 'Ce': 1.171, 'Pr': 1.170, 'Nd': 1.166, 'Sm': 1.160,
    'Eu': 1.158, 'Gd': 1.153, 'Tb': 1.151, 'Dy': 1.148, 'Ho': 1.146,
    'Er': 1.143, 'Tm': 1.142, 'Yb': 1.139, 'Lu': 1.137, 'Hf': 1.179,
    'Th': 1.138, 'U': 1.134,  'Ga': 1.344, 'As': 1.534, 'Mo': 1.500
}

# Phase-Specific Factors
FAC_PYRITE  = {'Fe': 2.148, 'Pb': 1.155, 'Zn': 1.490, 'Cu': 1.504, 'Co': 2.088, 'Ni': 2.091, 'As': 1.428}
FAC_SULFATE = {'Sr': 2.096, 'Ba': 1.699, 'Pb': 1.463, 'Ca': 3.395}
FAC_CARB    = {'Ca': 2.497, 'Mg': 3.468, 'Fe': 2.075, 'Mn': 2.093, 'Sr': 1.685, 'Ba': 1.437, 'Zn': 1.920, 'Na': 2.305}
FAC_OXIDE   = {'Fe': 1.430, 'Mn': 1.582, 'Co': 1.362, 'Ni': 1.273}

# ==========================================
# 3. CORE PROCESSING LOGIC
# ==========================================

def get_channel_data(ch_list, data_dict, ch_map_dict, mode='ppm'):
    """Helper function to load data with strict isotope name matching."""
    for ch in ch_list:
        try: name = ch.name()
        except: name = ch.name
        
        # Filter channel names
        if mode == 'ppm' and not name.endswith('_ppm'): continue
        if mode == 'prelim' and not (name.endswith('_Prelim_ppm') or name.endswith('_ppm')): continue
        
        # Regex to extract Element and Mass
        match = re.match(r"([A-Za-z]+)(\d+)", name)
        if not match: continue
        el, mass = match.group(1), match.group(2)
        
        # Avoid duplicates
        if el in data_dict and mode == 'prelim' and name.endswith('_ppm') and not name.endswith('_Prelim_ppm'):
            continue
        
        # Validate against target list
        if el in TARGET_ISOTOPES and TARGET_ISOTOPES[el] == mass:
            d = np.nan_to_num(ch.data()).astype(np.float64)
            data_dict[el] = d
            ch_map_dict[el] = ch

def run_workflow():
    print(">>> Starting: Adaptive Phase-Specific Normalization Workflow")
    
    # ----------------------------------------------------
    # Step 1: Preliminary Normalization
    # ----------------------------------------------------
    print("--- Step 1: Calculating Preliminary Concentrations ---")
    all_channels = data.timeSeriesList()
    data_map = {}
    channel_map = {}
    get_channel_data(all_channels, data_map, channel_map, mode='ppm')
    
    if not data_map:
        print("Error: No _ppm data found. Ensure semi-quantitative processing is complete.")
        return

    n_points = len(next(iter(data_map.values())))
    total_mass_prelim = np.zeros(n_points, dtype=np.float64)
    
    # Calculate initial total mass assuming default oxide/silicate composition
    for el, d_array in data_map.items():
        if el in EXCLUDE_FROM_SUM: continue
        total_mass_prelim += (d_array * FAC_DEFAULT.get(el, 1.0))
    
    total_mass_prelim[total_mass_prelim < 1e-6] = 1.0
    norm_factor_prelim = 1000000.0 / total_mass_prelim
    
    # Store preliminary data for ratio calculation
    prelim_data = {}
    for el, d in data_map.items():
        prelim_data[el] = d * norm_factor_prelim

    # ----------------------------------------------------
    # Step 2: Diagnostic Ratios
    # ----------------------------------------------------
    print("--- Step 2: Generating Diagnostic Ratio Channels ---")
    def get_p(k): return prelim_data.get(k, np.zeros(n_points))
    
    # Key elements with small offset to prevent division by zero
    Si = get_p('Si') + 0.001
    Ca = get_p('Ca'); Fe = get_p('Fe'); S = get_p('S')
    Al = get_p('Al'); Mn = get_p('Mn'); Sr = get_p('Sr'); Ba = get_p('Ba')

    # Calculate ratios
    r_Fe_Si = Fe / Si
    r_S_Si  = S / Si
    r_Ca_Si = Ca / Si
    r_Al_Ca = Al / (Ca + 0.001)
    r_Sr_Si = Sr / Si
    r_Ba_Si = Ba / Si
    r_Mn_Si = Mn / Si
    
    # Output diagnostic channels to Iolite UI for visual inspection
    try: t_axis = list(channel_map.values())[0].time()
    except: t_axis = None
    
    data.createTimeSeries("Diag_Fe_Si", data.Output, t_axis, r_Fe_Si)
    data.createTimeSeries("Diag_S_Si",  data.Output, t_axis, r_S_Si)
    data.createTimeSeries("Diag_Ca_Si", data.Output, t_axis, r_Ca_Si)
    data.createTimeSeries("Diag_Al_Ca", data.Output, t_axis, r_Al_Ca)
    data.createTimeSeries("Diag_Sr_Si", data.Output, t_axis, r_Sr_Si)
    data.createTimeSeries("Diag_Ba_Si", data.Output, t_axis, r_Ba_Si)
    data.createTimeSeries("Diag_Mn_Si", data.Output, t_axis, r_Mn_Si)

    # ----------------------------------------------------
    # Step 3: Phase Classification
    # ----------------------------------------------------
    print("--- Step 3: Applying Phase Classification Logic ---")

    # Initialize Phase Map (0 = Default/Matrix)
    phase_map = np.zeros(n_points, dtype=np.float64) 
    
    # Priority Logic: Distinct phases override the default (0)
    
    # A. Pyrite (High Fe & S)
    is_pyrite = (r_Fe_Si > USER_TH['Fe_Si_Min']) & (r_S_Si > USER_TH['S_Si_Min'])
    phase_map[is_pyrite] = 1.0
    
    # B. Sulfate (High Sr or Ba, if not Pyrite)
    is_sulfate = ((r_Sr_Si > USER_TH['Sr_Si_Min']) | (r_Ba_Si > USER_TH['Ba_Si_Min'])) & (phase_map == 0)
    phase_map[is_sulfate] = 2.0
    
    # C. Fe-Mn Oxide (High Mn, if not Pyrite/Sulfate)
    is_oxide = (r_Mn_Si > USER_TH['Mn_Si_Min']) & (phase_map == 0)
    phase_map[is_oxide] = 5.0
    
    # D. Carbonate (High Ca, Low Al, if not others)
    is_carb = (r_Ca_Si > USER_TH['Ca_Si_Min']) & \
              (r_Al_Ca < USER_TH['Al_Ca_Max']) & \
              (phase_map == 0)
    phase_map[is_carb] = 4.0
    
    # Output Phase Map
    data.createTimeSeries("Adapt_Phase_Class", data.Output, t_axis, phase_map)

    # ----------------------------------------------------
    # Step 4: Final Quantification
    # ----------------------------------------------------
    print("--- Step 4: Applying Lattice-Specific Factors ---")
    total_mass_final = np.zeros(n_points, dtype=np.float64)
    
    for el, d_array in data_map.items():
        if el in EXCLUDE_FROM_SUM: continue
        
        # 1. Apply Default Factor (Fallback)
        current_factor = np.full(n_points, FAC_DEFAULT.get(el, 1.0), dtype=np.float64)
        
        # 2. Override with Phase-Specific Factors where applicable
        if el in FAC_PYRITE: current_factor[phase_map == 1.0] = FAC_PYRITE[el]
        if el in FAC_SULFATE: current_factor[phase_map == 2.0] = FAC_SULFATE[el]
        if el in FAC_CARB: current_factor[phase_map == 4.0] = FAC_CARB[el]
        if el in FAC_OXIDE: current_factor[phase_map == 5.0] = FAC_OXIDE[el]
            
        total_mass_final += (d_array * current_factor)

    # Final Normalization
    total_mass_final[total_mass_final < 1e-6] = 1.0
    norm_factor_final = 1000000.0 / total_mass_final
    
    data.createTimeSeries("Step3_Final_Factor", data.Output, t_axis, norm_factor_final)
    
    # Write Final Output Channels
    count = 0
    print("Writing _Final_ppm channels...")
    for el, ch in channel_map.items():
        raw_d = np.nan_to_num(ch.data())
        final_d = np.ascontiguousarray(raw_d * norm_factor_final, dtype=np.float64)
        
        try: old_name = ch.name()
        except: old_name = ch.name
        base_name = old_name.split('_')[0] 
        new_name = base_name + "_Final_ppm"
        
        data.createTimeSeries(new_name, data.Output, t_axis, final_d)
        count += 1
        
    print(f"--- Workflow Complete ---")
    print(f"Generated {count} final channels.")
    print(">>> Check 'Adapt_Phase_Class' to verify mineral segmentation.")
    print(">>> Legend: 0=Matrix, 1=Pyrite, 2=Sulfate, 4=Carbonate, 5=Oxide")

# Execute
run_workflow()