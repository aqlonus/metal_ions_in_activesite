# Metal Ion Analysis Pipeline

A Python-based analysis pipeline for studying metal ion coordination in protein structures from PDB files.

## Features
- Identifies metal ions (Mg, Co, Ni, Zn, Mn, Ca, Fe, Cu, Cd) in PDB structures
- Analyzes coordination geometry and distances
- Categorizes results based on residue numbers
- Generates statistical analysis of metal-ligand distances

## Requirements
- Python 3.x
- BioPython
- NumPy
- Pandas

## Installation
```bash
pip install -r requirements.txt
```

## Usage
1. Place your PDB files in the `pdb_files` directory
2. Run the analysis:
```bash
python metal_ion_analysis_complete.py
```

3. Results will be organized in:
   - `results/` - Initial analysis results
   - `results/metal_1/` - First category of metal coordination
   - `results/metal_2/` - Second category of metal coordination

## Output
- CSV files containing metal coordination details
- Statistical analysis text files for each metal
- Categorized results based on residue numbers

## Directory Structure
