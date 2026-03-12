# Pixel-by-Pixel Adaptive Phase-Specific Normalization (PP-APSN)
# Supplementary Code for LA-ICP-TOF-MS Data Processing

## Overview

This repository contains a Python script for implementing the Pixel-by-Pixel Adaptive Phase-Specific Normalization (PP-APSN) protocol for LA-ICP-MS elemental mapping data. The protocol applies pixel-level mineral phase identification and phase-dependent stoichiometric constraints to accurately quantify multi-phase geological samples within the Iolite 4 platform.

---

## Requirements

### Software
- **Iolite** (v4.x or later): https://iolite-software.com
- **Python**: 3.7+ (included with Iolite)

### Dependencies
- `numpy` (standard with Iolite)
- `re` (Python standard library)

### Input Data
- Semi-quantitative LA-ICP-MS concentration channels in Iolite (format: `[Element][Mass]_ppm`)
- Example: `Ca44_ppm`, `Fe56_ppm`, `Si28_ppm`

---

## File

| File | Description |
|------|-------------|
| `Pixel-by-Pixel Adaptive Phase-Specific Normalization.py` | Complete PP-APSN workflow (all steps in one script) |

---

## Workflow Overview

The PP-APSN protocol executes four sequential steps automatically within a single script run:

```
Step 1: Preliminary Normalization (oxide-sum baseline)
         ↓
Step 2: Diagnostic Ratio Generation (element/Si channels)
         ↓
Step 3: Phase Classification (pixel-by-pixel segmentation)
         ↓
Step 4: Final Phase-Specific Quantification
```

---

## Execution Instructions

### How to Run

1. Load your LA-ICP-MS data in Iolite and ensure semi-quantitative processing is complete
2. Open the Python console in Iolite
3. Execute:
   ```python
   exec(open('Pixel-by-Pixel Adaptive Phase-Specific Normalization.py').read())
   ```
4. Monitor the console output for progress and completion messages

**Expected runtime**: 10–60 seconds depending on dataset size

---

## Step-by-Step Description

### Step 1: Preliminary Normalization

**Purpose**: Establish pixel-wise baseline concentrations using default oxide-sum normalization.

**Input**: Raw `[Element][Mass]_ppm` channels in Iolite

**Process**: All elements are summed using standard silicate oxide stoichiometric factors (`FAC_DEFAULT`) and normalized to 1,000,000 ppm (100 wt.%). This produces preliminary concentrations used solely for ratio diagnostics.

---

### Step 2: Diagnostic Ratio Generation

**Purpose**: Compute element/Si ratio maps for phase discrimination.

**Output channels** (visible in Iolite for inspection):

| Channel Name | Ratio | Target Phase |
|---|---|---|
| `Diag_Fe_Si` | Fe/Si | Pyrite |
| `Diag_S_Si` | S/Si | Pyrite/Sulfide |
| `Diag_Ca_Si` | Ca/Si | Carbonate |
| `Diag_Al_Ca` | Al/Ca | Silicate filter |
| `Diag_Sr_Si` | Sr/Si | Celestite |
| `Diag_Ba_Si` | Ba/Si | Barite |
| `Diag_Mn_Si` | Mn/Si | Fe-Mn Oxide |

> **Tip**: Visualize these channels in Iolite before running to verify that the classification thresholds in `USER_TH` are appropriate for your sample.

---

### Step 3: Phase Classification

**Purpose**: Assign each pixel to a mineral phase using a hierarchical priority scheme.

**Output channel**: `Adapt_Phase_Class`

#### Phase Codes

| Code | Phase | Priority |
|---|---|---|
| `1` | Pyrite / Sulfide | Highest |
| `2` | Sulfate (Celestite / Barite) | High |
| `5` | Fe-Mn Oxide | Medium |
| `4` | Carbonate | Low |
| `0` | Silicate / Oxide Matrix (Default) | Fallback |

#### Classification Criteria

**Phase 1 — Pyrite** (highest priority):  
`Fe/Si > Fe_Si_Min` **AND** `S/Si > S_Si_Min`  
Target minerals: FeS₂, PbS, ZnS, CuFeS₂

**Phase 2 — Sulfate**:  
`Sr/Si > Sr_Si_Min` **OR** `Ba/Si > Ba_Si_Min` *(if not Pyrite)*  
Target minerals: SrSO₄ (celestite), BaSO₄ (barite)

**Phase 5 — Fe-Mn Oxide**:  
`Mn/Si > Mn_Si_Min` *(if not Pyrite or Sulfate)*  
Target minerals: MnO₂, Fe₂O₃

**Phase 4 — Carbonate**:  
`Ca/Si > Ca_Si_Min` **AND** `Al/Ca < Al_Ca_Max` *(if not others)*  
Target minerals: CaCO₃, MgCO₃, FeCO₃  
*(Low Al/Ca filter excludes Ca-bearing silicates)*

**Phase 0 — Silicate/Matrix**:  
All remaining pixels (fallback)

---

### Step 4: Final Phase-Specific Quantification

**Purpose**: Apply phase-appropriate stoichiometric factors per pixel and renormalize to 100 wt.%.

**Output channels**: `[Element]_Final_ppm` for all input elements  
**Additional output**: `Step3_Final_Factor` — the normalization factor map per pixel

Each pixel receives element-specific stoichiometric factors according to its assigned phase (`Adapt_Phase_Class`). Pixels in the default phase use the standard oxide factors.

---

## User Configuration — Threshold Adjustment

Thresholds are defined in the `USER_TH` dictionary near the top of the script (Section 1). These values are initial estimates derived from preliminary ratio maps and **should be adjusted** based on your sample's mineralogy.

```python
USER_TH = {
    'Default_Phase': 'Silicate/Oxide Matrix (Standard Oxide Norm)',
    
    # Phase 1 — Pyrite
    'Fe_Si_Min': 0.015,   # Minimum Fe/Si for pyrite identification
    'S_Si_Min':  0.001,   # Minimum S/Si for sulfide presence
    
    # Phase 2 — Sulfate
    'Sr_Si_Min': 0.5,     # Minimum Sr/Si for celestite
    'Ba_Si_Min': 0.05,    # Minimum Ba/Si for barite
    
    # Phase 5 — Fe-Mn Oxide
    'Mn_Si_Min': 0.01,    # Minimum Mn/Si for oxide phases
    
    # Phase 4 — Carbonate
    'Ca_Si_Min': 1.5,     # Minimum Ca/Si for carbonate
    'Al_Ca_Max': 0.005    # Maximum Al/Ca to exclude aluminosilicates
}
```

**Recommended workflow for threshold tuning**:
1. Run the script once with default values
2. Inspect the `Diag_*` ratio channels in Iolite
3. Review `Adapt_Phase_Class` to assess segmentation quality
4. Adjust values in `USER_TH` and re-run as needed
5. Refer to **Table A3** in Supplementary Text S1 for recommended starting values

---

## Stoichiometric Factors

All stoichiometric factors are defined in **Table A4** of Supplementary Text S1. The script applies the following phase-specific factor sets:

| Phase | Example Factor | Key Elements |
|---|---|---|
| Default (Silicate) | SiO₂ = 2.139, Al₂O₃ = 1.889 | All elements (fallback) |
| Pyrite | FeS₂ = 2.148 | Fe, Pb, Zn, Cu, Co, Ni, As |
| Sulfate | SrSO₄ = 2.096, BaSO₄ = 1.699 | Sr, Ba, Pb, Ca |
| Carbonate | CaCO₃ = 2.497, MgCO₃ = 3.468 | Ca, Mg, Fe, Mn, Sr, Ba, Zn, Na |
| Fe-Mn Oxide | MnO₂ = 1.582, Fe₂O₃ = 1.430 | Fe, Mn, Co, Ni |

---

## Output Summary

After a successful run, the following channels are available in Iolite:

| Channel | Description | Recommended Use |
|---|---|---|
| `Diag_*` | Element/Si ratio maps | Phase boundary visualization |
| `Adapt_Phase_Class` | Phase segmentation map (codes 0–5) | Verify mineral identification |
| `Step3_Final_Factor` | Per-pixel normalization factor | Quality control |
| `[Element]_Final_ppm` | Final quantified concentrations | Quantitative analysis |

---

## Troubleshooting

**"No valid elemental channels found" / `data_map` is empty**  
→ Confirm that semi-quantitative processing is complete in Iolite  
→ Check that channels are named `[Element][Mass]_ppm` (e.g., `Fe56_ppm`, `Si28_ppm`)

**Unreasonable or unexpected phase classification**  
→ Inspect `Diag_*` channels to assess ratio distributions  
→ Adjust thresholds in `USER_TH` and re-run  
→ Refer to Table A3 for recommended starting values

**`_Final_ppm` channels not generated**  
→ Check the console for errors during Step 4  
→ Ensure the input `_ppm` channels exist and contain valid data

---

## Citation

If you use this code in your research, please cite:

*[Manuscript citation will be added upon publication]*

See Supplementary Text S1 for complete methodology, Table A3 (phase segmentation thresholds), and Table A4 (stoichiometric factors).

---

## Contact & Support

- **Corresponding Author**: Fei Li ([lifei@swpu.edu.cn])
- **Code Developer**: Fei Li; code development assisted by large language models Claude and Gemini for optimization and documentation

For issues, consult the troubleshooting section above or refer to Supplementary Text S1.

---

## License

This code is provided as supplementary material for academic research. Upon publication, it will be released under an open-source license (to be determined).

---

**Last updated**: March 2026  
**Version**: 2.0
