# membrane_permeability

## Permeability Analysis Script

The "calc_permeabilities.py" script (inside folder Pcalc) analyzes the molecular permeability using the ISDM method. 
## Features
📂 Extracts `.tpr` and `.xvg` files from umbrella sampling simulations.
⚙️ Runs WHAM analysis using GROMACS.
📊 Computes permeability using the ISDM method.
📈 Calculates jackknife error and standard error for reliability assessment.

### Requirements
- Python 3
- NumPy
- GROMACS (version 2021.2 or later)
- Umbrella sampling data files (`.tpr`, `.xvg`)
- isdm.py script (must be in the same directory for permeability calculations to work)

### Installation
Clone the repository:
```bash
git clone https://github.com/yourusername/permeability-analysis.git
cd permeability-analysis
```
Ensure dependencies are installed:
```bash
pip install numpy
```

### Usage
Run the script with the compound name as an argument:
```bash
python script.py COMPOUND_NAME
```
where `COMPOUND_NAME` is the name of the molecule being analyzed.

### Script Breakdown
#### 📌 Extracting .tpr Files
```python
get_tprs(base_dir, zzmin, zzmax, zzstep, output_file)
```
Collects `.tpr` files required for WHAM analysis.

#### 📌 Extracting Position & Force Data
```python
get_pulls_all(...)
```
Extracts `.xvg` files for all replicates.

```python
get_pulls_rep(...)
```
Extracts `.xvg` files for a specific replicate.

#### 📌 Running WHAM Analysis
```python
g_wham(...)
```
Runs WHAM to generate a potential of mean force (PMF) profile.

### Computing Permeability
```python
isdm_all(...)
isdm_rep(...)
```
Runs the ISDM method for permeability.

```python
extract_permeability(...)
```
Calculates and stores permeability values along with error estimates.

### Output
The script generates:
📝 `permeability_results.txt` with computed permeability values and errors.
📂 Intermediate `.xvg` and `.dat` files for WHAM and ISDM processing.

 

