# Pixel-by-Pixel Adaptive Phase-Specific Normalization (PP-APSN)
# Supplementary Code for LA-ICP-TOF-MS Data Processing

## Overview

This repository contains Python scripts for implementing the Pixel-by-Pixel Adaptive Phase-Specific Normalization (PP-APSN) protocol for LA-ICP-MS elemental mapping data. The protocol applies phase-dependent stoichiometric constraints to accurately quantify multi-phase geological samples.

## Requirements

### Software
- **Iolite** (v4.x or later): https://iolite-software.com
- **Python**: 3.7+ (included with Iolite)

### Dependencies
- `numpy` (standard with Iolite)
- `re` (Python standard library)

### Input Data
- Raw LA-ICP-MS concentration data channels (format: `[Element][Mass]_ppm`)
- Two reference CSV files (should be in the same directory):
  - `01-Adopted-isotopes.csv`
  - `02-Precise_Calibration_Table.csv`

## Workflow Overview

The PP-APSN protocol consists of three sequential steps:

```
Step 1: Preliminary Normalization
         ↓
Step 2: Phase Classification (Interactive)
         ↓
Step 3: Final Quantification
```

## Execution Instructions

### Step 1: Preliminary Imaging
**File**: `Step-1-Preliminary-Imaging.py`

**Purpose**: Perform preliminary oxide-sum normalization to generate internally consistent element/Si ratio diagnostics.

**Input**: 
- Raw concentration channels: `[Element][Mass]_ppm`

**Output**: 
- Preliminary concentrations: `[Element][Mass]_Prelim_ppm`

**How to run**:
1. Load your LA-ICP-MS data in Iolite
2. Open the Python console in Iolite
3. Execute: `exec(open('Step-1-Preliminary-Imaging.py').read())`
4. Wait for completion message

**Expected runtime**: 10-30 seconds (depending on dataset size)

---

### Step 2: Phase Classification (Interactive)
**File**: `Step-2-Phase-Classification-interactive.py`

**Purpose**: Perform hierarchical phase segmentation using element/Si ratios.

**Input**: 
- Preliminary concentrations from Step 1: `[Element][Mass]_Prelim_ppm`

**Output**: 
- Phase classification map: `Adapt_Phase_Class` (integer codes 1-5)
  - 1 = Pyrite/Sulfide
  - 2 = Sulfate (celestite/barite)
  - 3 = Fe-Mn Oxide
  - 4 = Carbonate
  - 5 = Silicate/Matrix (default)

**⚠️ IMPORTANT - Threshold Adjustment**:

The script uses default threshold values from **Table A3** in Supplementary Text S1. These values may need adjustment based on your specific sample mineralogy.

**To adjust thresholds**:
1. Run Step 1 first
2. Visualize element/Si ratio distributions in Iolite (e.g., Fe/Si, Ca/Si, Mn/Si)
3. Open `Step-2-Phase-Classification-interactive.py` in a text editor
4. Locate the `R_TH` dictionary (around line 50-80)
5. Modify threshold values based on your observations:

```python
R_TH = {
    'Fe_Si_Pyrite': 0.015,    # Adjust for pyrite identification
    'S_Si_Pyrite':  0.001,    # Adjust for sulfide phases
    'Sr_Si_Sulf':   0.5,      # Adjust for celestite (SrSO4)
    'Ba_Si_Sulf':   0.05,     # Adjust for barite (BaSO4)
    'Mn_Si_Oxide':  0.01,     # Adjust for Fe-Mn oxides
    'Ca_Si_Carb':   1.5,      # Adjust for carbonate phases
    'Al_Ca_Limit':  0.005     # Adjust to exclude aluminosilicates
}
```

**How to run**:
1. Ensure Step 1 is complete
2. Execute: `exec(open('Step-2-Phase-Classification-interactive.py').read())`
3. Review the console output showing pixel counts for each phase
4. A dialog box will appear showing classification results
5. Click "Yes" to proceed or "No" to review/adjust thresholds

**Expected runtime**: 10-30 seconds

---

### Step 3: Final Quantitative Imaging
**File**: `Step-3-Final-Quantitative-Imaging.py`

**Purpose**: Apply phase-specific stoichiometric normalization for final quantification.

**Input**: 
- Preliminary concentrations: `[Element][Mass]_Prelim_ppm`
- Phase classification map: `Adapt_Phase_Class`

**Output**: 
- Final quantified concentrations: `[Element][Mass]_Final_ppm`
- All elements normalized to 100 wt.% compound sum per phase

**How to run**:
1. Ensure Steps 1 and 2 are complete
2. Execute: `exec(open('Step-3-Final-Quantitative-Imaging.py').read())`
3. Check console for summary statistics

**Expected runtime**: 10-30 seconds

---

## Phase Classification Hierarchy

The classification follows a strict priority order to minimize misclassification:

1. **Pyrite/Sulfide** (Highest priority)
   - Criteria: High Fe/Si AND high S/Si
   - Target minerals: FeS₂, PbS, ZnS, CuS

2. **Sulfate** (Medium-high priority)
   - Criteria: High Sr/Si OR high Ba/Si
   - Target minerals: SrSO₄ (celestite), BaSO₄ (barite)

3. **Fe-Mn Oxide** (Medium priority)
   - Criteria: High Mn/Si
   - Target minerals: MnO₂, Fe₂O₃

4. **Carbonate** (Low priority)
   - Criteria: High Ca/Si AND low Al/Ca
   - Target minerals: CaCO₃, MgCO₃, FeCO₃
   - Note: Low Al/Ca filter excludes Ca-feldspars

5. **Silicate/Matrix** (Fallback)
   - All remaining pixels
   - Default normalization using silicate oxide factors

---

## Stoichiometric Factors

All stoichiometric factors are defined in **Table A4** of Supplementary Text S1. The script automatically applies phase-appropriate factors:

- **Pyrite**: FeS₂ (factor = 2.148), PbS, ZnS, etc.
- **Sulfate**: SrSO₄ (factor = 2.096), BaSO₄ (1.699), etc.
- **Fe-Mn Oxide**: MnO₂ (factor = 1.582), Fe₂O₃ (1.430), etc.
- **Carbonate**: CaCO₃ (factor = 2.497), MgCO₃ (3.468), etc.
- **Silicate**: SiO₂ (factor = 2.139), Al₂O₃ (1.889), etc.

---

## Troubleshooting

### "No valid elemental channels found"
- Check that your raw data channels are named: `[Element][Mass]_ppm`
- Example: `Ca44_ppm`, `Fe56_ppm`, `Si28_ppm`

### "Phase map 'Adapt_Phase_Class' not found"
- You must run Step 2 before Step 3
- Check that Step 2 completed successfully

### "No preliminary concentration channels found"
- You must run Step 1 before Steps 2 or 3
- Check for `[Element][Mass]_Prelim_ppm` channels in Iolite

### Unreasonable phase classification results
- Adjust threshold values in `R_TH` dictionary (Step 2)
- Visualize element/Si ratios to identify appropriate cutoffs
- Refer to Table A3 for recommended starting values

---

## Output Data

After completing all three steps, you will have:

1. **Preliminary concentrations**: `[Element][Mass]_Prelim_ppm`
   - Use for ratio diagnostics only

2. **Phase classification**: `Adapt_Phase_Class`
   - Integer map (1-5) showing mineral phase per pixel

3. **Final concentrations**: `[Element][Mass]_Final_ppm`
   - Use for quantitative analysis and interpretation
   - Normalized to 100 wt.% compound sum per phase

---

## Citation

If you use this code in your research, please cite:

[Manuscript citation will be added upon publication]

---

## Reference Tables

- **Table A3**: Phase segmentation thresholds (element/Si ratios)
- **Table A4**: Stoichiometric factors by mineral phase

See Supplementary Text S1 for complete tables and methodology.

---

## Contact & Support

For questions or issues related to this code, please contact:

- **Corresponding Author**: [Fei Li] ([lifei@swpu.edu.cn])
- **Code Developer**: [Fei Li; Code development assisted by large language models Claude and Gemini
for optimization and documentation] ([lifei@swpu.edu.cn])

For questions or issues related to this code:
- Check the troubleshooting section above
- Refer to Supplementary Text S1 for detailed methodology
- Contact the corresponding author

---

## License

This code is provided as supplementary material for academic research. Upon publication, it will be released under an open-source license (to be determined).

---

**Last updated**: February 2026
**Version**: 1.0
