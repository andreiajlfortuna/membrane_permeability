#  Membrane Permeability of Halogenated Drugs

<p align="center">
  <img width="700" alt="Graphical abstract: halogen and hydrogen bonds between a halogenated drug and a POPC bilayer" src="https://github.com/user-attachments/assets/325b3f45-6fa8-4509-ad42-24f00183b2af" />
</p>




This repository has part of the workflow used in our work *Impact of Halogen Anisotropy on the Membrane Permeability of Halogenated Drugs* (DOI: [10.1021/acs.jmedchem.6c01696](https://doi.org/10.1021/acs.jmedchem.6c01696)).

Bromazepam was chosen as example but the same methodology can be applied to other halogenated molecules, adjusting the EP parameters accordingly.

## Methods overview

### Ligand parametrization
- **Geometry optimization:** B3LYP/6-311G(d,p) (Gaussian 09).
- **Charges:** RESP charges fitted to ESP calculated at HF/6-31G(d) (6-311G(d) for iodine, with an MK radius of 2.3 Å).  σ-hole emulated with an extra point (EP) of charge, placed at a distance R<sub>min</sub> from the halogen along the C–X axis (180°) and implemented as a GROMACS type 2 virtual site. A "no EP" set of RESP charges is also computed for comparison.
- **Force field:** GAFF
- **Multiconformational RESP** for molecules where an intramolecular HB with the halogen can occur and changes the σ-hole (furosemide, metolazone).

### Simulations (GROMACS 2021.2)
- **System:** one solute in a POPC bilayer (128 lipids, 5652 TIP3P waters), Lipid14 force field, 310 K.
- **Minimization and Equilibration:** energy minimization (steepest descent), then three NVT steps with decreasing restraints.
- **Unbiased MD:** 5 independent replicates of 250 ns (NPT; v-rescale thermostat, semi-isotropic Berendsen barostat, PME, P-LINCS, 2 fs time step).
- **Steered MD:** solute pulled from the aqueous phase to the bilayer center to generate the starting configurations.
- **Umbrella sampling:** 19 windows spaced 2 Å apart (0–36 Å from the bilayer center), 3 replicates of 150 ns per window (first 50 ns discarded).
- **PMF:** computed with WHAM.

### Analysis
- **Diffusivity D(z):** position-dependent, from the autocorrelation of z-fluctuations in each window.
- **Permeability (ISDM):** inhomogeneous solubility-diffusion model combining the PMF and D(z).
- **Permeability ranking:** P<sub>ΔG<sub>ranking</sub></sub> = k<sub>B</sub>T / (ΔG<sub>max</sub> − ΔG<sub>min</sub>).
- **Insertion depth:** calculated with MembIT, accounting for local membrane deformation.
- **XB and HB detection (MDAnalysis):** geometric criteria, with XB defined by distance below the sum of van der Waals radii and R–X···Y angle above 140°, and HB by D···A distance below 3 Å and D–H···A angle above 150°.
- **Interaction energies:** electrostatic and van der Waals ligand–membrane contributions along the insertion coordinate, for the whole molecule and by atom or group.
- **Statistics:** Jackknife error estimates; EP vs no EP compared using non-overlapping error bars or a composite error criterion.
- **Descriptors:** consensus logP and TPSA from SwissADME, compared with experimental Caco-2 P<sub>app</sub> values.

