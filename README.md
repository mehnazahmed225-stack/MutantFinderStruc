MutantFinderStruc is a desktop application designed for researchers investigating the structural and energetic effects of missense mutations in proteins. The application enables users to load a protein structure in PDB format and automatically generate descriptive statistics, including amino‑acid composition, chain‑level distributions, and residue charge properties. Its core functionality performs an exhaustive in‑silico missense mutation scan using PyRosetta, in which each residue is systematically substituted with all other amino acids, followed by structural relaxation. For each mutant, the tool computes Rosetta energy values, ΔΔG relative to the wild type, and Cα RMSD to quantify both energetic and structural perturbations. Results are exported in tabular form and complemented by interactive visualizations, including RMSD–ΔΔG plots, RMSD distributions, and Ramachandran analysis for structural quality assessment. By integrating mutation generation, scoring, and visualization into a single graphical interface, MutantFinderStruc streamlines otherwise manual Rosetta workflows and enables faster identification of mutations with potentially significant structural or stability impacts.

Installation instructions:

Important note: has to be run on linux-64.

1. download these files and extract as needed
2. Open in VS code
3. Open terminal in VS code(make sure conda in enabled)
4. Enter command: conda env create -f environment.yml
5. for pyrosetta seperately, use command:pip install pyrosetta --find-links https://west.rosettacommons.org/pyrosetta/quarterly/release
6. Enter command: conda activate MutantFinderStruc
7. Select MutantFinderStruc as interpreter for VS code.
8. Run code to confirm if all libraries are installed
9. To make code an app: pyinstaller MutantFinderStruc.py
