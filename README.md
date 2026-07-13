# MutantFinder

MutantFinder is a desktop application for protein structure analysis and computational mutagenesis built using PyRosetta, BioPython, PySide6, and Matplotlib. The software enables users to analyze amino acid composition, perform saturation mutagenesis, estimate mutation-induced energetic changes, generate mutant structures, and visualize structural properties through interactive plots.

---

## Overview

The program combines protein structure analysis with large-scale single-point mutation screening.

Given a protein structure in PDB format, MutantFinderStruc:

1. Extracts sequence and composition information.
2. Computes amino acid statistics.
3. Performs systematic saturation mutagenesis.
4. Relaxes wild-type and mutant structures using Rosetta FastRelax.
5. Calculates energetic and structural metrics for each mutant.
6. Saves mutant structures as PDB files.
7. Generates graphical visualizations of mutation effects.
8. Provides a Ramachandran plot viewer for backbone conformation analysis.

---

## Features

### Structural Analysis

- Counts amino acids within each chain.
- Determines amino acid composition of the entire protein.
- Calculates the number of:
  - Polar residues
  - Nonpolar residues
  - Positively charged residues
  - Negatively charged residues

### Saturation Mutagenesis

- Performs systematic single-point mutagenesis.
- Generates all 19 non-native substitutions for each residue.
- Excludes residue 1 from mutation scanning by design.
- Saves all generated mutant structures as individual PDB files.

### Stability Analysis

- Relaxes structures using Rosetta FastRelax.
- Calculates Rosetta full-atom energy scores.
- Computes mutation-induced energy changes (ΔΔG).
- Calculates backbone Cα RMSD values.

### Visualization

- Interactive RMSD vs ΔΔG scatter plot.
- RMSD distribution histogram.
- Ramachandran plot viewer with secondary-structure coloring.

---

## Methods

### Protein Structure Parsing

Protein structures are loaded from PDB files using BioPython's `PDBParser`.

For each chain:

- Standard amino acid residues are identified.
- Three-letter residue codes are converted to one-letter codes.
- Amino acid sequences are reconstructed.
- Residue counts are determined.

Outputs include:

- Residue count per chain
- Amino acid composition per chain

---

### Amino Acid Composition Analysis

The complete protein sequence is extracted using PyRosetta.

Counts are calculated for all twenty standard amino acids and reported for the entire protein.

Residues are additionally classified into four physicochemical groups:

#### Polar

```
S T N Q Y C
```

#### Nonpolar

```
A V L I M F W P G
```

#### Positively Charged

```
K R H
```

#### Negatively Charged

```
D E
```

---

### Wild-Type Structure Preparation

The input structure is loaded into PyRosetta and relaxed using Rosetta FastRelax.

The Rosetta all-atom score function is used:

```python
scorefxn = get_fa_scorefxn()
```

FastRelax optimizes side-chain conformations and reduces unfavorable steric interactions to produce a minimized wild-type reference structure.

The energy of the relaxed wild-type structure is recorded as:

```text
WT_dG
```

---

### Saturation Mutagenesis

Single-point saturation mutagenesis is performed for every residue beginning at position 2.

For each position:

- The native residue is identified.
- All nineteen alternative amino acids are generated.
- The original amino acid is excluded.

Mutation names follow standard notation:

```text
A25V
L42R
F98Y
```

where:

```text
NativeResidue + Position + MutantResidue
```

---

### Mutant Structure Generation

For every mutation:

1. A copy of the relaxed wild-type structure is created.
2. The target residue is replaced using Rosetta's `MutateResidue` mover.
3. The mutant structure is relaxed using FastRelax.
4. Energetic and structural measurements are calculated.
5. The mutant structure is saved as a PDB file.

Workflow:

```text
Wild Type
    ↓
Mutation
    ↓
FastRelax
    ↓
Energy Evaluation
    ↓
Save Mutant Structure
```

---

### Energy Calculations

The Rosetta score function is used to calculate mutant structure energies.

For each relaxed mutant:

```text
Mut_dG
```

is calculated.

Mutation effects are estimated as:

```text
ΔΔG = Mut_dG − WT_dG
```

where:

- WT_dG = relaxed wild-type energy
- Mut_dG = relaxed mutant energy

#### Interpretation

```text
Negative ΔΔG
    More energetically favorable than wild type

Near-zero ΔΔG
    Similar energetic behavior to wild type

Positive ΔΔG
    Less energetically favorable than wild type
```

### Important Note

Reported ΔΔG values are differences in Rosetta Energy Units (REU) and should be interpreted as relative energetic changes predicted by Rosetta rather than experimentally measured free energies.

---

### Structural Deviation Analysis

Structural changes induced by mutation are quantified using Cα RMSD.

RMSD is calculated using:

```python
CA_rmsd(wt_pose, mut_pose)
```

#### Interpretation

```text
Low RMSD
    Minimal structural perturbation

High RMSD
    Larger conformational change
```

---

### Ramachandran Analysis

Backbone conformational properties are visualized using a Ramachandran plot.

Phi (φ) and Psi (ψ) torsion angles are extracted for all protein residues.

Secondary structure is assigned using DSSP:

```python
Dssp(pose)
```

Residues are colored according to their secondary structure:

```text
Red   = α-Helix
Blue  = β-Sheet
Gray  = Loop/Coil
```

Reference regions corresponding to common secondary-structure conformations are displayed on the plot.

---

## Output Files

### Mutation Results

```text
mutation_results.csv
```

Contains:

- Mutation name
- Mutant energy (Mut_dG)
- RMSD
- ΔΔG

### Mutant Structures

All generated mutant structures are stored in:

```text
mutants/
```

Example:

```text
mutants/
├── A25V.pdb
├── A25L.pdb
├── A25W.pdb
└── ...
```

### Composition Statistics

```text
Count_aa_per_chain
Count_aa_type_in_chain
Amino_acid_count_whole_protein
```

---

## Workflow Summary

```text
Load PDB
    ↓
Sequence and Composition Analysis
    ↓
Wild-Type Relaxation
    ↓
Generate Single Mutants
    ↓
Relax Mutants
    ↓
Calculate Energy Scores
    ↓
Calculate ΔΔG
    ↓
Calculate RMSD
    ↓
Save Mutant Structures
    ↓
Visualize Results
```

---

## Limitations

- Only standard amino acids are supported.
- Only single-point mutations are evaluated.
- Residue 1 is intentionally excluded from mutation scanning.
- Energies are reported in Rosetta Energy Units (REU).
- Predictions are computational estimates and should be experimentally validated.
- Runtime increases substantially with protein size and mutation count.

---

## Dependencies

- PyRosetta
- BioPython
- Pandas
- PySide6
- Matplotlib
- Seaborn
- mplcursors

---

## Author

MutantFinderStruc was developed as a protein engineering and computational structural biology tool for mutation screening, stability assessment, and structural analysis.
