#!/usr/bin/env python3
from pyrosetta.io import pose_from_pdb
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
from pyrosetta import get_fa_scorefxn
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.rosetta.core.scoring import rama_prepro
from pyrosetta.rosetta.core.scoring.dssp import Dssp
import pandas as pd
import sys
from PySide6.QtWidgets import (
    QApplication,
    QMainWindow,
    QPushButton,
    QTextEdit,
    QFileDialog,
    QVBoxLayout,
    QWidget,
    QLabel,
    QSpinBox
)
from Bio.PDB import PDBParser
import os
from collections import Counter
from PySide6.QtWidgets import QMainWindow, QWidget, QVBoxLayout
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
import mplcursors
from pathlib import Path
import pyrosetta


def get_pyrosetta_database():
    if hasattr(sys, "_MEIPASS"):
        return Path(sys._MEIPASS) / "pyrosetta_database"
    else:
        return Path(pyrosetta.__file__).parent / "database"


db_path = get_pyrosetta_database()

pyrosetta.init(extra_options=f"-database {db_path}")

amino_acids = [
    "A","R","N","D","C","E","Q","G","H","I",
    "L","K","M","F","P","S","T","W","Y","V"
]

polar = ["S", "T", "N", "Q", "Y", "C"]
nonpolar = ["A", "V", "L", "I", "M", "F", "W", "P", "G"]
positive = ["K", "R", "H"]
negative = ["D", "E"]

three_to_one = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLU": "E",
    "GLN": "Q",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V"
}

one_to_three = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "E": "GLU",
    "Q": "GLN",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL"
}

output_dir = "mutants"
os.makedirs(output_dir, exist_ok=True)

class PreliminaryDataCollection:
    def __init__(self, pdb_path):
        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure("protein", pdb_path)
        self.pdb_path = pdb_path
    
    
    def summary(self):
        df3 = pd.DataFrame(columns=["Chain", "AA count"])
        chains = []
        sequences = {}
        for model in self.structure:
            for chain in model:
                seq = []
                for residue in chain:
                    if residue.id[0] != " ":
                        continue
                    resname = residue.get_resname()
                    if resname in three_to_one:
                        seq.append(three_to_one[resname])
                if not seq:
                    continue
                chain_id = chain.id
                chains.append(chain_id)
                sequences[chain_id] = seq
                df3.loc[len(df3)] = [chain_id, len(seq)]
            break
        df4 = pd.DataFrame(index=amino_acids, columns=chains)
        for chain_id, seq in sequences.items():
            counts = Counter(seq)
            for aa in amino_acids:
                df4.loc[aa, chain_id] = counts.get(aa, 0)
        return df3, df4

    def AminoAcidCount(self):
        pose = pyrosetta.pose_from_pdb(self.pdb_path)
        sequence = pose.sequence()
        amino_acid_counts = {aa: sequence.count(aa) for aa in amino_acids}
        df1 = pd.DataFrame(
            amino_acid_counts.items(),
            columns=["Amino Acid", "Count"]
        )
        return df1
    
    def charge(self):
        pose = pyrosetta.pose_from_pdb(self.pdb_path)
        sequence = pose.sequence()
        polar_list = []
        nonpolar_list = []
        positive_list = []
        negative_list = []
        unknown = []
        for aa in sequence:
            if aa in polar:
                polar_list.append(aa)
            elif aa in nonpolar:
                nonpolar_list.append(aa)
            elif aa in positive:
                positive_list.append(aa)
            elif aa in negative:
                negative_list.append(aa)
            else:
                unknown.append(aa)
        Polar = len(polar_list)
        Nonpolar = len(nonpolar_list)
        Positive = len(positive_list)
        Negative = len(negative_list)
        return  Polar, Nonpolar, Positive, Negative

class Mutations:
    def __init__(self, pdb_path, progress_callback):
        self.pdb_path = pdb_path
        self.df_overall = pd.DataFrame(columns = ["Mutation_name", "Mut_dG", "RMSD", "ddG"])
        self.progress_callback = progress_callback

    def make_mutated_pdbs(self):
        # Load WT structure
        wt_pose = pose_from_pdb(self.pdb_path)
        # Score function (once)
        scorefxn = get_fa_scorefxn()
        # Relax WT ONCE
        relax = FastRelax()
        relax.set_scorefxn(scorefxn)
        relax.apply(wt_pose)
        # WT reference energy
        WT_dG = scorefxn(wt_pose)
        sequence = wt_pose.sequence()
        for position, native_aa in enumerate(sequence[1:], start=2):
            for mut_aa in amino_acids:
                if mut_aa == native_aa:
                    continue
                # Clone from RELAXED WT
                mut_pose = wt_pose.clone()
                # Mutate residue
                def safe_mutate_residue(pose, resi, new_res):
                    mut = MutateResidue(resi, new_res)
                    mut.apply(pose)
                mut_aa_3 = one_to_three[mut_aa]

                safe_mutate_residue(mut_pose, position, mut_aa_3)

                #relax
                relax.apply(mut_pose)
                Mut_dG = scorefxn(mut_pose)
                ddG = Mut_dG - WT_dG
                rmsd = CA_rmsd(wt_pose, mut_pose)
                mutation_name = f"{native_aa}{position}{mut_aa}"
                self.df_overall.loc[len(self.df_overall)] = [
                    mutation_name,
                    Mut_dG,
                    rmsd,
                    ddG
                ]
                line = (
                    f"{mutation_name:8s} | "
                    f"Mut_dG={Mut_dG:7.2f} | "
                    f"RMSD={rmsd:5.2f} | "
                    f"ΔΔG={ddG:6.2f}"
                )
                if self.progress_callback:
                    self.progress_callback(line)
                    QApplication.processEvents()
                outfile = os.path.join(
                    output_dir,
                    f"{native_aa}{position}{mut_aa}.pdb"
                )
                mut_pose.dump_pdb(outfile)

class SubWidget1(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Graph 1")
        self.resize(600, 500)

        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)

        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(self.canvas)

        self.plot()

    def plot(self):
        df = pd.read_csv("mutation_results.csv")

        ax = self.figure.add_subplot(111)
        ax.clear()

        # ✅ Use matplotlib scatter (required for mplcursors)
        scatter = ax.scatter(
            df["RMSD"],
            df["ddG"],
            s=30,
            alpha=0.7
        )

        ax.set_title("RMSD vs ΔΔG")
        ax.set_xlabel("RMSD")
        ax.set_ylabel("ΔΔG")

        # ✅ Attach cursor to the scatter artist
        cursor = mplcursors.cursor(scatter, hover=True)

        @cursor.connect("add")
        def on_add(sel):
            idx = sel.index
            sel.annotation.set_text(
                f"{df.iloc[idx]['Mutation_name']}\n"
                f"RMSD: {df.iloc[idx]['RMSD']:.2f}\n"
                f"ΔΔG: {df.iloc[idx]['ddG']:.2f}"
            )
            sel.annotation.get_bbox_patch().set(alpha=0.85)

        self.canvas.draw()

    

class SubWidget2(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.setWindowTitle("Graph 2: RMSD Histogram")
        self.resize(600, 500)

        central = QWidget(self)
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)

        self.bin_label = QLabel("Number of bins:")
        self.bin_selector = QSpinBox()
        self.bin_selector.setRange(5, 200)
        self.bin_selector.setValue(30)
        self.bin_selector.valueChanged.connect(self.plot)

        layout.addWidget(self.bin_label)
        layout.addWidget(self.bin_selector)

        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(self.canvas)

        self.ax = self.figure.add_subplot(111)

        # Load data once
        self.df = pd.read_csv("mutation_results.csv")

        self.plot()

    def plot(self):
        bins = self.bin_selector.value()

        self.ax.clear()

        sns.histplot(
            data=self.df,
            x="RMSD",
            bins=bins,
            ax=self.ax
        )

        self.ax.set_title(f"RMSD Distribution (bins = {bins})")
        self.ax.set_xlabel("RMSD")
        self.ax.set_ylabel("Count")

        self.canvas.draw()

class SubWidget3(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Ramachandran Plot Viewer")
        self.resize(700, 700)

        central = QWidget(self)
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)

        # Load button
        self.load_btn = QPushButton("Load PDB and Plot Ramachandran")
        self.load_btn.clicked.connect(self.load_pdb)
        layout.addWidget(self.load_btn)

        # Matplotlib canvas
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(self.canvas)

        self.ax = self.figure.add_subplot(111)

    def load_pdb(self):
        pdb_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select PDB file",
            "",
            "PDB Files (*.pdb)"
        )
        if not pdb_path:
            return

        pose = pose_from_pdb(pdb_path)
        self.plot_ramachandran(pose)
    
    def draw_ramachandran_regions(self):
        regions = [
            (-90, -80, 60, 60, "red", "α-helix"),
            (-150, 90, 60, 90, "blue", "β-sheet"),
            (30, 0, 60, 90, "green", "left α")
        ]

        for x, y, w, h, color, _ in regions:
            self.ax.add_patch(
                Rectangle((x, y), w, h, color=color, alpha=0.12, zorder=0))

    def plot_ramachandran(self, pose):
        phis, psis, colors = [], [], []
        # Compute secondary structure
        dssp = Dssp(pose)
        dssp.insert_ss_into_pose(pose)
        ss = pose.secstruct()

        color_map = {
            'H': 'red',
            'E': 'blue',
            'L': 'gray'
        }

        for i in range(1, pose.total_residue() + 1):
            if pose.residue(i).is_protein():
                phis.append(pose.phi(i))
                psis.append(pose.psi(i))
                colors.append(color_map.get(ss[i - 1], 'black'))

        self.ax.clear()
        self.ax.scatter(phis, psis, c=colors, s=12, alpha=0.7)

        self.draw_ramachandran_regions()
        self.ax.set_xlim(-180, 180)
        self.ax.set_ylim(-180, 180)
        self.ax.set_xlabel("Phi (°)")
        self.ax.set_ylabel("Psi (°)")
        self.ax.set_title("Ramachandran Plot")
        self.ax.legend(handles=[
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=8, label='Helix'),
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=8, label='Sheet'),
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=8, label='Loop')
        ], loc='upper right')

        self.canvas.draw()        

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("MutantFinderStruc")
        self.resize(800, 800)
        self.setMinimumSize(900, 600)
        self.setWindowTitle("MutantFinderStruc")
        self.current_file = None
        # Widgets
        self.text_display = QTextEdit()
        self.text_display.setReadOnly(True)

        self.load_button = QPushButton("Load File")
        self.load_button.clicked.connect(self.load_file)

        self.load_button2 = QPushButton("Get Data")
        self.load_button2.clicked.connect(self.get_data)

        self.text_display2 = QTextEdit()
        self.text_display2.setReadOnly(True)

        self.load_button3 = QPushButton("Mutate")
        self.load_button3.clicked.connect(self.mutate)

        self.text_display3 = QTextEdit()
        self.text_display3.setReadOnly(True)

        self.load_button4 = QPushButton("Graphs")
        self.load_button4.clicked.connect(self.graphics)
        
        # Layout
        layout = QVBoxLayout()
        layout.addWidget(self.load_button)
        layout.addWidget(self.text_display)
        layout.addWidget(self.load_button2)
        layout.addWidget(self.text_display2)
        layout.addWidget(self.load_button3)
        layout.addWidget(self.text_display3)
        layout.addWidget(self.load_button4)
        # Central widget
        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)

    def load_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Open File",
            "",
            "PDB Files (*.pdb);;All Files (*)"
        )
        if not file_path:
            return
        self.current_file = file_path
        with open(file_path, "r") as f:
            self.text_display.setText(f.read())

    def get_data(self):
        if not self.current_file:
            self.text_display2.append("\n No PDB file loaded.")
            return
        data = PreliminaryDataCollection(self.current_file)
        df3, df4 = data.summary()
        self.text_display2.append(f"\nThe counts of amino acids per chain\n {df3} \nthe counts of amino acid types per chain\n {df4}" )
        Polar, Nonpolar, Positive, Negative = data.charge()
        self.text_display2.append(f"\n The overall counts for polar, nonpolar, positive, and negative amino acids \nPolar: {Polar} | Nonpolar {Nonpolar} | Positive {Positive} | Negative {Negative}")
        AminoAcidCount = data.AminoAcidCount()
        self.text_display2.append(f"\n Amino Acid Counts in whole protein: \n {AminoAcidCount}")
        df3.to_csv("Count_aa_per_chain")
        df4.to_csv("Count_aa_type_in_chain")
        AminoAcidCount.to_csv("Amino_acid_count_whole_protein")

    def mutate(self):
        if not self.current_file:
            self.text_display3.append("\n No PDB file loaded.")
            return
        self.text_display3.clear()
        self.text_display3.append("Mutation scan started...this may take a while so feel free to grab a coffee.\nBy the way, the calculations for one mutation might take anywhere between three seconds to 30 minutes. Bear with me\n")
        QApplication.processEvents()
        mutator = Mutations(
            self.current_file,
            progress_callback=self.text_display3.append
        )
        mutator.make_mutated_pdbs()
        self.text_display3.append("\nScan complete.")
        mutator.df_overall.to_csv("mutation_results.csv", index=False)

    def graphics(self):
        if not self.current_file:
            return
        if not hasattr(self, "graph1"):
            self.graph1 = SubWidget1(self)
            self.graph2 = SubWidget2(self)
            self.graph3 = SubWidget3()

        self.graph1.show()
        self.graph2.show()
        self.graph3.show()

        return None

def main():
    app = QApplication(sys.argv)
    window1 = MainWindow()
    window1.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()