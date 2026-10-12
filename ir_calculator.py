"""
IR Calculator PyMOL Plugin (PyQt5 Version)
==========================================
Calculate IR wavenumber and reduced mass for chemical bonds

Installation:
1. In PyMOL: Plugin → Plugin Manager → Install New Plugin
2. Select this file (ir_calculator.py)
3. Restart PyMOL
4. Access via: Plugin → IR Calculator

Authors: Ella Wydemon, Donna H. Perygin, Phillip B. Ogden
License: MIT (see LICENSE)
Repository: https://github.com/dperygin/IR_Calculator_PyMOL
Version: 1.0 (PyQt5)
"""

try:
    from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, 
                                 QPushButton, QLineEdit, QTextEdit, QTabWidget,
                                 QWidget, QListWidget, QGroupBox, QMessageBox, QComboBox)
    from PyQt5.QtCore import Qt
    from PyQt5.QtGui import QFont
except ImportError:
    try:
        from PyQt4.QtGui import (QDialog, QVBoxLayout, QHBoxLayout, QLabel,
                                QPushButton, QLineEdit, QTextEdit, QTabWidget,
                                QWidget, QListWidget, QGroupBox, QMessageBox, QFont, QComboBox)
        from PyQt4.QtCore import Qt
    except ImportError:
        print("Error: PyQt5 or PyQt4 is required for this plugin")
        raise

import math
from pymol import cmd

# Keep a persistent reference to the plugin dialog when launched from PyMOL
plugin_dialog = None

# Physical constants
SPEED_OF_LIGHT = 2.99792458e10  # cm/s
AMU_TO_KG = 1.66053906660e-27   # kg/amu

# Common atomic masses (amu)
ATOMIC_MASSES = {
    # Hydrogen isotopes
    'H':   1.008,    # Protium — Primary Isotope (99.99%) 
    'D':   2.014,    # Deuterium (²H) — NMR solvents, isotope labeling, KIE studies
    'T':   3.016,    # Tritium (³H) — radiotracer labeling
    # Carbon isotopes
    'C':   12.000,   # Carbon-12 — Primary Isotope (98.94%)
    '13C': 13.003,   # Carbon-13 — ¹³C NMR, mechanistic labeling
    '14C': 14.003,   # Carbon-14 — radiotracer, dating
    # Nitrogen isotopes
    'N':   14.003,   # Nitrogen-14 — Primary Isotope (99.64%)
    '15N': 15.000,   # Nitrogen-15 — ¹⁵N NMR, metabolic labeling
    # Oxygen isotopes
    'O':   15.995,   # Oxygen-16 — Primary Isotope (99.76%)
    '18O': 17.999,   # Oxygen-18 — mechanistic labeling (ester hydrolysis, etc.)
    # Fluorine (one stable isotope; ¹⁸F included for PET radiochemistry)
    'F':   18.998,   # Fluorine-19 — Primary Isotope (100.00%)
    '18F': 18.001,   # Fluorine-18 — PET radiotracers
    # Phosphorus (one stable isotope; ³²P for biochemistry)
    'P':   30.974,   # Phosphorus-31 — Primary Isotope (100.00%)
    '32P': 31.974,   # Phosphorus-32 — radiotracer (biochemistry / nucleotides)
    # Sulfur isotopes
    'S':   31.972,   # Sulfur-32 — Primary Isotope (94.99%)
    '34S': 33.967,   # Sulfur-34 — isotope labeling, SIMS
    # Chlorine isotopes
    'Cl':   34.968,  # Chlorine-35 — Most Common Isotope (75.76%)
    '37Cl': 36.966,  # Chlorine-37 — mass spec isotope pattern (24.24%)
    # Bromine isotopes
    'Br':   78.918,  # Bromine-79 — Majority Isotope (50.69%)
    '81Br': 80.916,  # Bromine-81 — mass spec isotope pattern (49.31%)
    # Iodine (one stable isotope; ¹²⁵I for biochemistry)
    'I':    126.904, # Iodine-127 — Primary Isotope (100.00%)
    '125I': 124.905, # Iodine-125 — radiotracer (protein labeling)
    # Silicon isotopes
    'Si':   27.977,  # Silicon-28 — Primary Isotope (92.22%)
    '29Si': 28.976,  # Silicon-29 — ²⁹Si NMR
    # Boron isotopes
    'B':   11.009,   # Boron-11 — ¹¹B NMR Primary Isotope (80.17%)
    '10B': 10.013,   # Boron-10 — BNCT (boron neutron capture therapy)
    # Selenium
    'Se':  79.917  # Selenium-80 — ¹¹B NMR Most Common Isotope (49.61%)
}

# Common bond force constants (N/m)
# NOTE: These are empirical values for idealized bonds. Actual values depend on chemical context.
# This calculator uses the harmonic oscillator approximation (single harmonic potential).
# Isotope-labeled bonds use the same force constant as their unlabeled counterparts;
# the frequency shift arises entirely from the reduced mass.
BOND_DATA = {
    # C-H / C-D by hybridization
    'C-H':   {'k': 477,  'atoms': ['C', 'H'],  'description': 'Carbon-Hydrogen sp3 (alkane -> ~2950 cm⁻¹)',   'multiplicity': 1},
    '=C-H':  {'k': 500,  'atoms': ['C', 'H'],  'description': 'Carbon-Hydrogen sp2 (vinyl / aromatic -> ~3020 cm⁻¹)', 'multiplicity': 1},
    '≡C-H':  {'k': 597,  'atoms': ['C', 'H'],  'description': 'Carbon-Hydrogen sp (terminal alkyne -> ~3300 cm⁻¹)',   'multiplicity': 1},
    'C-D':   {'k': 500,  'atoms': ['C', 'D'],  'description': 'Carbon-Deuterium sp3 (~2100 cm⁻¹)',         'multiplicity': 1},
    '=C-D':  {'k': 530,  'atoms': ['C', 'D'],  'description': 'Carbon-Deuterium sp2 (~2250 cm⁻¹)',         'multiplicity': 1},
    '≡C-D':  {'k': 590,  'atoms': ['C', 'D'],  'description': 'Carbon-Deuterium sp (~2575 cm⁻¹)',          'multiplicity': 1},
    # Carbon-Carbon
    'C-C':   {'k': 450,  'atoms': ['C', 'C'],  'description': 'Carbon-Carbon (single)',                    'multiplicity': 1},
    'C=C':   {'k': 960,  'atoms': ['C', 'C'],  'description': 'Carbon-Carbon (double)',                    'multiplicity': 2},
    'C≡C':   {'k': 1600, 'atoms': ['C', 'C'],  'description': 'Carbon-Carbon (triple)',                    'multiplicity': 3},
    'C:C':   {'k': 700,  'atoms': ['C', 'C'],  'description': 'Carbon-Carbon (aromatic)',                  'multiplicity': 1.5},
    # Carbon-Oxygen
    'C-O':   {'k': 530,  'atoms': ['C', 'O'],  'description': 'Carbon-Oxygen (single)',                    'multiplicity': 1},
    'C=O':   {'k': 1200, 'atoms': ['C', 'O'],  'description': 'Carbon-Oxygen (double, carbonyl)',          'multiplicity': 2},
    # Carbon-Nitrogen
    'C-N':   {'k': 480,  'atoms': ['C', 'N'],  'description': 'Carbon-Nitrogen (single)',                  'multiplicity': 1},
    'C=N':   {'k': 870,  'atoms': ['C', 'N'],  'description': 'Carbon-Nitrogen (double, imine)',           'multiplicity': 2},
    'C≡N':   {'k': 1930, 'atoms': ['C', 'N'],  'description': 'Carbon-Nitrogen (triple, nitrile)',         'multiplicity': 3},
    # O-H / O-D
    'O-H':   {'k': 745,  'atoms': ['O', 'H'],  'description': 'Oxygen-Hydrogen (alcohol / acid)',          'multiplicity': 1},
    'O-D':   {'k': 745,  'atoms': ['O', 'D'],  'description': 'Oxygen-Deuterium (~2660 cm⁻¹)',             'multiplicity': 1},
    # N-H / N-D
    'N-H':   {'k': 630,  'atoms': ['N', 'H'],  'description': 'Nitrogen-Hydrogen (amine / amide)',         'multiplicity': 1},
    'N-D':   {'k': 630,  'atoms': ['N', 'D'],  'description': 'Nitrogen-Deuterium (~2400 cm⁻¹)',           'multiplicity': 1},
    # Carbon-Halogen
    'C-F':   {'k': 570,  'atoms': ['C', 'F'],  'description': 'Carbon-Fluorine',                           'multiplicity': 1},
    'C-Cl':  {'k': 390,  'atoms': ['C', 'Cl'], 'description': 'Carbon-Chlorine',                           'multiplicity': 1},
    'C-Br':  {'k': 330,  'atoms': ['C', 'Br'], 'description': 'Carbon-Bromine',                            'multiplicity': 1},
    'C-I':   {'k': 270,  'atoms': ['C', 'I'],  'description': 'Carbon-Iodine',                             'multiplicity': 1},
    # Nitrogen-Nitrogen
    'N-N':   {'k': 350,  'atoms': ['N', 'N'],  'description': 'Nitrogen-Nitrogen (single, hydrazine)',     'multiplicity': 1},
    'N=N':   {'k': 980,  'atoms': ['N', 'N'],  'description': 'Nitrogen-Nitrogen (double, azo)',           'multiplicity': 2},
    'N≡N':   {'k': 2250, 'atoms': ['N', 'N'],  'description': 'Nitrogen-Nitrogen (triple, diazo / azide)','multiplicity': 3},
    # Nitrogen-Oxygen
    'N-O':   {'k': 310,  'atoms': ['N', 'O'],  'description': 'Nitrogen-Oxygen (single, hydroxylamine)',   'multiplicity': 1},
    'N=O':   {'k': 1550, 'atoms': ['N', 'O'],  'description': 'Nitrogen-Oxygen (double, nitroso / nitro)', 'multiplicity': 2},
    # Oxygen-Oxygen
    'O-O':   {'k': 210,  'atoms': ['O', 'O'],  'description': 'Oxygen-Oxygen (single, peroxide)',          'multiplicity': 1},
    'O=O':   {'k': 1140, 'atoms': ['O', 'O'],  'description': 'Oxygen-Oxygen (double)',                    'multiplicity': 2},
    # Sulfur
    'S-H':   {'k': 410,  'atoms': ['S', 'H'],  'description': 'Sulfur-Hydrogen (thiol)',                   'multiplicity': 1},
    'S-D':   {'k': 410,  'atoms': ['S', 'D'],  'description': 'Sulfur-Deuterium (~1800 cm⁻¹)',             'multiplicity': 1},
    'C-S':   {'k': 260,  'atoms': ['C', 'S'],  'description': 'Carbon-Sulfur (single, thioether)',         'multiplicity': 1},
    'C=S':   {'k': 800,  'atoms': ['C', 'S'],  'description': 'Carbon-Sulfur (double, thione)',            'multiplicity': 2},
    'S-S':   {'k': 430,  'atoms': ['S', 'S'],  'description': 'Sulfur-Sulfur (single, disulfide)',         'multiplicity': 1},
    'S-O':   {'k': 370,  'atoms': ['S', 'O'],  'description': 'Sulfur-Oxygen (single, sulfonate)',         'multiplicity': 1},
    'S=O':   {'k': 950,  'atoms': ['S', 'O'],  'description': 'Sulfur-Oxygen (double, sulfoxide)',         'multiplicity': 2},
    # Phosphorus
    'P-H':   {'k': 320,  'atoms': ['P', 'H'],  'description': 'Phosphorus-Hydrogen (phosphine)',           'multiplicity': 1},
    'P-D':   {'k': 320,  'atoms': ['P', 'D'],  'description': 'Phosphorus-Deuterium (~1630 cm⁻¹)',         'multiplicity': 1},
    'P-C':   {'k': 280,  'atoms': ['P', 'C'],  'description': 'Phosphorus-Carbon (single)',                'multiplicity': 1},
    'P-O':   {'k': 400,  'atoms': ['P', 'O'],  'description': 'Phosphorus-Oxygen (single, ester)',         'multiplicity': 1},
    'P=O':   {'k': 980,  'atoms': ['P', 'O'],  'description': 'Phosphorus-Oxygen (double)',                'multiplicity': 2},
    'P-N':   {'k': 360,  'atoms': ['P', 'N'],  'description': 'Phosphorus-Nitrogen (single)',              'multiplicity': 1},
    'P=N':   {'k': 880,  'atoms': ['P', 'N'],  'description': 'Phosphorus-Nitrogen (double, iminophosphorane)', 'multiplicity': 2},
    # Silicon
    'Si-H':  {'k': 320,  'atoms': ['Si', 'H'],  'description': 'Silicon-Hydrogen (silane)',                'multiplicity': 1},
    'Si-D':  {'k': 320,  'atoms': ['Si', 'D'],  'description': 'Silicon-Deuterium (~1560 cm⁻¹)',           'multiplicity': 1},
    'Si-C':  {'k': 360,  'atoms': ['Si', 'C'],  'description': 'Silicon-Carbon (organosilicon)',           'multiplicity': 1},
    'Si-O':  {'k': 600,  'atoms': ['Si', 'O'],  'description': 'Silicon-Oxygen (siloxane)',                'multiplicity': 1},
    'Si-F':  {'k': 830,  'atoms': ['Si', 'F'],  'description': 'Silicon-Fluorine (fluorosilane)',          'multiplicity': 1},
    'Si-Cl': {'k': 370,  'atoms': ['Si', 'Cl'], 'description': 'Silicon-Chlorine (chlorosilane)',          'multiplicity': 1},
    # Boron
    'B-H':   {'k': 310,  'atoms': ['B', 'H'],  'description': 'Boron-Hydrogen (borane)',                   'multiplicity': 1},
    'B-D':   {'k': 310,  'atoms': ['B', 'D'],  'description': 'Boron-Deuterium (~1800 cm⁻¹)',              'multiplicity': 1},
    'B-C':   {'k': 290,  'atoms': ['B', 'C'],  'description': 'Boron-Carbon (organoboron)',                'multiplicity': 1},
    'B-O':   {'k': 510,  'atoms': ['B', 'O'],  'description': 'Boron-Oxygen (boronic acid / ester)',       'multiplicity': 1},
    'B-N':   {'k': 440,  'atoms': ['B', 'N'],  'description': 'Boron-Nitrogen (borazine / amine-borane)',  'multiplicity': 1},
    # Selenium
    'Se-H':  {'k': 330,  'atoms': ['Se', 'H'], 'description': 'Selenium-Hydrogen (selenol)',               'multiplicity': 1},
    'Se-D':  {'k': 330,  'atoms': ['Se', 'D'], 'description': 'Selenium-Deuterium (~1590 cm⁻¹)',           'multiplicity': 1},
    'C-Se':  {'k': 230,  'atoms': ['C', 'Se'], 'description': 'Carbon-Selenium (selenide)',                'multiplicity': 1},
}

# Approximate bond lengths for order estimation (Å)
# Note: Isotope substitution (H → D) does not meaningfully change bond length;
# use the same ('X','H') entry for deuterated analogues.
BOND_LENGTH_TEMPLATES = {
    ('C', 'C'):  {'single': 1.54, 'double': 1.34, 'triple': 1.20, 'aromatic': 1.40},
    ('C', 'N'):  {'single': 1.47, 'double': 1.25, 'triple': 1.16},
    ('C', 'O'):  {'single': 1.43, 'double': 1.23},
    ('C', 'H'):  {'single': 1.09},
    ('C', 'S'):  {'single': 1.82, 'double': 1.60},
    ('C', 'F'):  {'single': 1.35},
    ('C', 'Cl'): {'single': 1.77},
    ('C', 'Br'): {'single': 1.94},
    ('C', 'I'):  {'single': 2.10},
    ('C', 'Se'): {'single': 1.95},
    ('N', 'H'):  {'single': 1.01},
    ('N', 'N'):  {'single': 1.45, 'double': 1.25, 'triple': 1.10},
    ('N', 'O'):  {'single': 1.44, 'double': 1.22},
    ('O', 'H'):  {'single': 0.96},
    ('O', 'O'):  {'single': 1.48},
    ('S', 'H'):  {'single': 1.34},
    ('S', 'S'):  {'single': 2.05},
    ('S', 'O'):  {'single': 1.65, 'double': 1.43},
    ('P', 'O'):  {'single': 1.62, 'double': 1.48},
    ('P', 'H'):  {'single': 1.42},
    ('P', 'C'):  {'single': 1.84},
    ('P', 'N'):  {'single': 1.69, 'double': 1.54},
    ('Si', 'H'): {'single': 1.48},
    ('Si', 'C'): {'single': 1.87},
    ('Si', 'O'): {'single': 1.63},
    ('Si', 'F'): {'single': 1.60},
    ('Si', 'Cl'):{'single': 2.02},
    ('B', 'H'):  {'single': 1.19},
    ('B', 'C'):  {'single': 1.56},
    ('B', 'O'):  {'single': 1.36},
    ('B', 'N'):  {'single': 1.44},
    ('Se', 'H'): {'single': 1.47},
}

# Common Calculator Functions: Reduced Mass, IR Wavenumber, Force Constant
def calculate_reduced_mass(mass1, mass2):
    """Calculate reduced mass for a two-body system."""
    if mass1 <= 0 or mass2 <= 0:
        raise ValueError("Masses must be positive values")
    return (mass1 * mass2) / (mass1 + mass2)


def calculate_wavenumber(force_constant, reduced_mass):
    """Calculate IR wavenumber from force constant and reduced mass."""
    if force_constant <= 0 or reduced_mass <= 0:
        raise ValueError("Force constant and reduced mass must be positive")
    
    reduced_mass_kg = reduced_mass * AMU_TO_KG
    frequency_hz = (1 / (2 * math.pi)) * math.sqrt(force_constant / reduced_mass_kg)
    wavenumber = frequency_hz / SPEED_OF_LIGHT
    return wavenumber


def calculate_force_constant(wavenumber, reduced_mass):
    """Calculate force constant from wavenumber and reduced mass."""
    if wavenumber <= 0 or reduced_mass <= 0:
        raise ValueError("Wavenumber and reduced mass must be positive")
    
    reduced_mass_kg = reduced_mass * AMU_TO_KG
    frequency_hz = wavenumber * SPEED_OF_LIGHT
    force_constant = reduced_mass_kg * (2 * math.pi * frequency_hz) ** 2
    return force_constant


class IRCalculatorDialog(QDialog):
    """Main dialog for IR Calculator"""
    
    def __init__(self, parent=None):
        super(IRCalculatorDialog, self).__init__(parent)
        self.setWindowTitle("IR Wavenumber Calculator")
        self.setMinimumSize(750, 650)
        
        # Bond picking state
        self.custom_functional_group_data = None
        
        self.setup_ui()
        
    def setup_ui(self):
        """Setup the user interface"""
        layout = QVBoxLayout()
        
        # Create tab widget
        self.tabs = QTabWidget()
        
        # Add tabs
        self.tabs.addTab(self.create_common_bond_tab(), "Common Bonds")
        self.tabs.addTab(self.create_custom_tab(), "Custom Calculation")
        self.tabs.addTab(self.create_reverse_tab(), "Find Force Constant")
        
        layout.addWidget(self.tabs)
        
        # Close button
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        layout.addWidget(close_btn)
        
        self.setLayout(layout)
    
    def create_common_bond_tab(self):
        """Create common bond calculation tab"""
        widget = QWidget()
        layout = QVBoxLayout()
        
        # Title
        title = QLabel("Calculate from Common Bonds")
        title_font = QFont()
        title_font.setPointSize(14)
        title_font.setBold(True)
        title.setFont(title_font)
        title.setAlignment(Qt.AlignCenter)
        layout.addWidget(title)
        
        # Approximation note
        note = QLabel("Note: Uses harmonic oscillator approximation. Values are estimates.")
        note.setStyleSheet("color: #666666; font-size: 9pt; font-style: italic;")
        note.setAlignment(Qt.AlignCenter)
        layout.addWidget(note)
        
        # Bond list group
        group = QGroupBox("Select Bond Type")
        group_layout = QVBoxLayout()
        
        self.bond_list = QListWidget()
        mono_font = QFont("Courier")
        mono_font.setPointSize(9)
        self.bond_list.setFont(mono_font)
        
        # Populate list of common bonds
        for bond_name, data in BOND_DATA.items():
            display_text = f"{bond_name:8s} - {data['description']:35s} (k={data['k']} N/m)"
            self.bond_list.addItem(display_text)
        
        group_layout.addWidget(self.bond_list)
        group.setLayout(group_layout)
        layout.addWidget(group)
        
        # Calculate button
        calc_btn = QPushButton("Calculate Wavenumber")
        calc_btn.setStyleSheet("background-color: #4CAF50; color: white; font-weight: bold; padding: 10px;")
        calc_btn.clicked.connect(self.calculate_common_bond)
        layout.addWidget(calc_btn)
        
        # Results of Wavenumber Calculation from Common Bond
        self.common_results = QTextEdit()
        self.common_results.setReadOnly(True)
        self.common_results.setFont(mono_font)
        self.common_results.setMaximumHeight(200)
        layout.addWidget(self.common_results)
        
        widget.setLayout(layout)
        return widget
    
    def create_custom_tab(self):
        """Create custom calculation tab"""
        widget = QWidget()
        layout = QVBoxLayout()
        
        # Title
        title = QLabel("Calculate Wavenumber from Custom Values")
        title_font = QFont()
        title_font.setPointSize(14)
        title_font.setBold(True)
        title.setFont(title_font)
        title.setAlignment(Qt.AlignCenter)
        layout.addWidget(title)
        
        # Approximation note
        note = QLabel("Harmonic oscillator approximation. Uses standard atomic masses (protium, ¹²C, etc).")
        note.setStyleSheet("color: #666666; font-size: 9pt; font-style: italic;")
        note.setAlignment(Qt.AlignCenter)
        layout.addWidget(note)

        # Two options for Custom Bonds:
            # Model: Pick 2 adjacent and bonded atoms from a model
            # Custom: Select two elements and the Force Constant for related bond

        # ==== Model Bonds =====
        # Use PyMol Picker to select two adjacent atoms.
        # The two atoms must be selected such that they are pk1 and pk2.  Deselect doesn't work.
        # When loaded, the atoms and force constant are displaced in a message and presented in the input parameters frame.
        
        # Bond Selection button
        pick_layout = QHBoxLayout()
        pick_bond_btn = QPushButton("📍 Select Bond in PyMOL")
        pick_bond_btn.setStyleSheet("background-color: #9C27B0; color: white; font-weight: bold; padding: 10px; font-size: 11px;")
        pick_bond_btn.clicked.connect(self.start_bond_picking)
        pick_layout.addWidget(pick_bond_btn)
        
        load_bond_btn = QPushButton("⬇️ Load Selected Bond")
        load_bond_btn.setStyleSheet("background-color: #E91E63; color: white; font-weight: bold; padding: 10px; font-size: 11px;")
        load_bond_btn.clicked.connect(self.process_picked_bond)
        pick_layout.addWidget(load_bond_btn)
        
        layout.addLayout(pick_layout)
        
        # Status label for picking mode
        self.pick_status_label = QLabel("")
        self.pick_status_label.setAlignment(Qt.AlignCenter)
        self.pick_status_label.setStyleSheet("color: #9C27B0; font-weight: bold;")
        layout.addWidget(self.pick_status_label)
        
        # Input group
        group = QGroupBox("Input Parameters")
        group_layout = QVBoxLayout()

        # Alternatively, a bond can be generated using a pair of elements and the intervening bond's Force Constant
        # When using the Model Picker method, the elements and force constant are presented in the input parameters
                
        # Element 1 via a drop-down
        h1a = QHBoxLayout()
        h1a.addWidget(QLabel("Element 1:"))
        self.custom_elem1 = QComboBox()
        self.custom_elem1.addItems(['C', 'H', 'D', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I'])
        self.custom_elem1.setMaximumWidth(80)
        self.custom_elem1.currentTextChanged.connect(self.update_custom_masses)
        h1a.addWidget(self.custom_elem1)
        h1a.addWidget(QLabel("   Mass (amu):"))
        self.custom_mass1 = QLineEdit()
        self.custom_mass1.setMaximumWidth(100)
        self.custom_mass1.setReadOnly(True)
        self.custom_mass1.setStyleSheet("background-color: #f0f0f0;")
        h1a.addWidget(self.custom_mass1)
        h1a.addStretch()
        group_layout.addLayout(h1a)
        
        # Element 2 via a drop-down
        h2a = QHBoxLayout()
        h2a.addWidget(QLabel("Element 2:"))
        self.custom_elem2 = QComboBox()
        self.custom_elem2.addItems(['C', 'H', 'D', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I'])
        self.custom_elem2.setMaximumWidth(80)
        self.custom_elem2.currentTextChanged.connect(self.update_custom_masses)
        h2a.addWidget(self.custom_elem2)
        h2a.addWidget(QLabel("   Mass (amu):"))
        self.custom_mass2 = QLineEdit()
        self.custom_mass2.setMaximumWidth(100)
        self.custom_mass2.setReadOnly(True)
        self.custom_mass2.setStyleSheet("background-color: #f0f0f0;")
        h2a.addWidget(self.custom_mass2)
        h2a.addStretch()
        group_layout.addLayout(h2a)

        # Element 1-2 via a drop-down
        h12 = QHBoxLayout()
        h12.addWidget(QLabel("Bond Order:"))
        self.custom_bond = QComboBox()
        self.custom_bond.addItems(['single', 'double', 'triple'])
        self.custom_bond.setMaximumWidth(80)
        h12.addWidget(self.custom_bond)
        h12.addWidget(QLabel("C-Adjacent Bond:"))
        self.custom_bond_adj = QComboBox()
        self.custom_bond_adj.addItems(['non-Carbon', 'single', 'double', 'triple'])
        self.custom_bond_adj.setMaximumWidth(80)
        h12.addWidget(self.custom_bond_adj)
        h12.addStretch()
        group_layout.addLayout(h12)
        
        # Enter the Force Constant manually or using the load method.
        h3 = QHBoxLayout()
        h3.addWidget(QLabel("Force Constant (N/m):"))
        self.custom_k = QLineEdit()
        self.custom_k.setMaximumWidth(80)
        h3.addWidget(self.custom_k)
        h3_bond_btn = QPushButton("Update Force Constant")
        h3_bond_btn.clicked.connect(self.update_force_constant)
        h3.addWidget(h3_bond_btn)
        h3.addStretch()
        group_layout.addLayout(h3)
        
        group.setLayout(group_layout)
        layout.addWidget(group)
        
        # Initialize masses
        self.update_custom_masses()
        
        # Calculate button
        # The input parameters of element masses and force constant to calculate the wavenumber.
        # The calculation and results are presented via messages.
        
        calc_btn = QPushButton("Calculate")
        calc_btn.setStyleSheet("background-color: #2196F3; color: white; font-weight: bold; padding: 10px;")
        calc_btn.clicked.connect(self.calculate_custom)
        layout.addWidget(calc_btn)
        
        # Results showing Wavenumber Calculation
        mono_font = QFont("Courier")
        mono_font.setPointSize(9)
        self.custom_results = QTextEdit()
        self.custom_results.setReadOnly(True)
        self.custom_results.setFont(mono_font)
        self.custom_results.setMaximumHeight(200)
        layout.addWidget(self.custom_results)
        
        layout.addStretch()
        widget.setLayout(layout)
        return widget
    
    def create_reverse_tab(self):
        """Create Tab for the Reverse Calculation of the Force Constant from the IR Wavenumber"""
        widget = QWidget()
        layout = QVBoxLayout()
        
        # Title
        title = QLabel("Calculate Force Constant from Wavenumber")
        title_font = QFont()
        title_font.setPointSize(14)
        title_font.setBold(True)
        title.setFont(title_font)
        title.setAlignment(Qt.AlignCenter)
        layout.addWidget(title)
        
        # Approximation note
        note = QLabel("Reverse calculation from observed wavenumber. Harmonic approximation may limit accuracy.")
        note.setStyleSheet("color: #666666; font-size: 9pt; font-style: italic;")
        note.setAlignment(Qt.AlignCenter)
        layout.addWidget(note)
        
        # Input group
        group = QGroupBox("Input Parameters")
        group_layout = QVBoxLayout()
        
        # Mass of Atom 1
        h1 = QHBoxLayout()
        h1.addWidget(QLabel("Mass of Atom 1 (amu):"))
        self.reverse_mass1 = QLineEdit()
        self.reverse_mass1.setMaximumWidth(150)
        h1.addWidget(self.reverse_mass1)
        h1.addStretch()
        group_layout.addLayout(h1)
        
        # Mass of Atom 2
        h2 = QHBoxLayout()
        h2.addWidget(QLabel("Mass of Atom 2 (amu):"))
        self.reverse_mass2 = QLineEdit()
        self.reverse_mass2.setMaximumWidth(150)
        h2.addWidget(self.reverse_mass2)
        h2.addStretch()
        group_layout.addLayout(h2)
        
        # Wavenumber manually entered from IR
        h3 = QHBoxLayout()
        h3.addWidget(QLabel("Wavenumber (cm⁻¹):"))
        self.reverse_wn = QLineEdit()
        self.reverse_wn.setMaximumWidth(150)
        h3.addWidget(self.reverse_wn)
        h3.addStretch()
        group_layout.addLayout(h3)
        
        group.setLayout(group_layout)
        layout.addWidget(group)
        
        # Calculate button to determine Force Constant
        calc_btn = QPushButton("Calculate")
        calc_btn.setStyleSheet("background-color: #FF9800; color: white; font-weight: bold; padding: 10px;")
        calc_btn.clicked.connect(self.calculate_reverse)
        layout.addWidget(calc_btn)
        
        # Results showing Force Constant Calculation
        mono_font = QFont("Courier")
        mono_font.setPointSize(9)
        self.reverse_results = QTextEdit()
        self.reverse_results.setReadOnly(True)
        self.reverse_results.setFont(mono_font)
        self.reverse_results.setMaximumHeight(200)
        layout.addWidget(self.reverse_results)
        
        layout.addStretch()
        widget.setLayout(layout)
        return widget
    # ====== Functions to Determine Force Constant based on Selected Bond from Model ======
    # The following are a series of support functions for determining the correct Force Constant based:
    #    Two bonded atoms and their elements
    #    Bond Order/Type
    #    Hybridization of the atoms
    #    Adjacent Atoms and Bonds
    # The calculation is designed for providing harmonic oscillator approximations based on the common bonds and environments.
    # The Force Constant is generated based on a list and is not rigorously determined nor are all adjacent bond types considered.
    
    def update_custom_masses(self):
        """Update mass fields based on selected elements in Custom Tab"""
        self.custom_functional_group_data = None
        elem1 = self.custom_elem1.currentText()
        elem2 = self.custom_elem2.currentText()
        
        if elem1 in ATOMIC_MASSES:
            self.custom_mass1.setText(f"{ATOMIC_MASSES[elem1]:.3f}")
        if elem2 in ATOMIC_MASSES:
            self.custom_mass2.setText(f"{ATOMIC_MASSES[elem2]:.3f}")

    def update_force_constant(self):
        elem1 = self.custom_elem1.currentText()
        elem2 = self.custom_elem2.currentText()
        bond_order = self.custom_bond.currentText()
        bond_order_adj = self.custom_bond_adj.currentText()

        if bond_order_adj == 'single':
            hybrid = 'sp3'
        elif bond_order_adj == 'double' or bond_order_adj == 'aromatic':
            hybrid = 'sp2'
        elif bond_order_adj == 'triple':
            hybrid = 'sp'
        else:
            hybrid = None
            
        suggested_k = self.get_force_constant_for_bond(elem1, elem2, bond_order, hybrid)
        self.custom_k.setText(str(suggested_k))
        
    def get_explicit_bond_order(self, state):
        """Read explicit bond order from PyMOL model data when available."""
        # Attempt to pull Bond Order directly from the Model using the Picked Atoms
        try:
            model = cmd.get_model("pk1 or pk2", state)
            try:
                order = getattr(model.bond[0], 'order', None)
            except:
                print("Bond order unable to be retrived from selection.")
                order = None
            if order is not None:
                return int(order)
            else:
                return None
        except Exception as e:
            print(f"Bond order detection failed: {e}")
        return None

    def get_explicit_aromatic_flag(self):
        """Check whether either picked atom is marked aromatic in PyMOL."""
        try:
            model = cmd.get_model("pk1 or pk2")
            for atom in model.atom:
                if getattr(atom, 'aromatic', False):
                    print("Explicit aromatic flag detected for pk1/pk2")
                    return True
        except Exception as e:
            print(f"Aromatic flag detection failed: {e}")

        # If the model does not expose aromatic flags, assume False.
        return False

    def get_neighbor_atoms(self, pick_sel, exclude_sel, radius=1.9):
        """Collect nearby atoms around a picked atom as a neighbor estimate."""
        # For each Picked Atom, collect all of the other bonded atoms. 
        neighbors = []
        try:
            cmd.iterate(
                f"{pick_sel} around {radius} and not ({pick_sel} or {exclude_sel})",
                "neighbors.append((index, elem, name, ID))",
                space={'neighbors': neighbors}
            )
        except Exception as e:
            print(f"Neighbor detection failed: {e}")
        return neighbors

    def get_neighbors_by_index(self, atom_index, radius=1.9, exclude_sel=None):
        """Collect nearby atoms around a specific atom index."""
        # Using the indicies for each selected Atom, collect all atoms bonded to the selected pair.
        # Deprecated due to index / ID mismatches
        neighbors = []
        try:
            selector = f"id {atom_index} around {radius} and not id {atom_index}"
            if exclude_sel:
                selector += f" and not ({exclude_sel})"
            cmd.iterate(
                selector,
                "neighbors.append((index, elem, name, ID))",
                space={'neighbors': neighbors}
            )
        except Exception as e:
            print(f"Neighbor detection failed for atom {atom_index}: {e}")
        return neighbors

    def get_neighbors_by_ID(self, atom_ID, radius=1.9, exclude_sel=None):
        """Collect nearby atoms around a specific atom ID."""
        # Using the ID for each selected Atom, collect all atoms bonded to the selected pair. 
        neighbors = []
        try:
            selector = f"id {atom_ID} around {radius} and not id {atom_ID}"
            if exclude_sel:
                selector += f" and not ({exclude_sel})"
            cmd.iterate(
                selector,
                "neighbors.append((index, elem, name, ID))",
                space={'neighbors': neighbors}
            )
        except Exception as e:
            print(f"Neighbor detection failed for atom {atom_ID}: {e}")
        return neighbors

    def estimate_bond_order_by_length(self, elem1, elem2, dist, aromatic=False):
        """Estimate bond order from atom pair and bond length if explicit order is absent."""
        # Determine the Bond Order based on
        #   Bond Length from Model
        #   Bond Length from Template List using Element Identity
        
        # Find Atom Pair among Bond Length Template,
        #   if element pair is not present, then assume the Bond is Single
        if elem1 == 'D':
            elem1 = 'H'

        if elem2 == 'D':
            elem2 = 'H'

        if aromatic and {elem1, elem2} == {'C', 'C'}:
            return 'aromatic'

        pair = tuple(sorted([elem1, elem2]))
        templates = BOND_LENGTH_TEMPLATES.get(pair)
        if not templates:
            print(f"No length template for {pair}; defaulting to single")
            return 'single'

        # Using Bond Distance and Bond Length Templates
        #   ID which Bond Order matches has the closest matching Bond Length to the entered Length
        best = None
        best_diff = None
        for order, template_dist in templates.items():
            if order == 'aromatic' and not aromatic:
                continue
            diff = abs(dist - template_dist)
            if best_diff is None or diff < best_diff:
                best_diff = diff
                best = order

        if best is None:
            print(f"No applicable order found for {pair}; defaulting to single")
            return 'single'

        print(f"Length-based bond order estimate for {elem1}-{elem2}: {best} (dist={dist:.2f} Å; aromatic={aromatic})")
        return best or 'single'

    def infer_hybridization(self, elem, neighbor_count, aromatic=False, bond_order=None):
        """Infer approximate hybridization from connectivity and bond order."""
        # Determine hybridization for submitted Atom using 
        #   Bond Order
        #   Neighbor Count
        #   Otherwise it is unknown
        # Note: This program is designed for organic IR spectra.
        #   Hybridization is determined intended for Carbon, Oxygen, terminal Halogens, and Nitrogen.
        #   Halogens should be single-bonded to a Carbon and thus return sp3.
        if elem in ('H', 'D'):
            return 'unhybridized'

        if aromatic:
            return 'sp2'

        if bond_order in (3, 'triple'):
            return 'sp'
        elif bond_order in (2, 'double'):
            return 'sp2'
        elif bond_order in ('aromatic',):
            return 'sp2'


        if elem == "C":
            if neighbor_count >= 4:
                return 'sp3'
            elif neighbor_count == 3:
                return 'sp2'
            elif neighbor_count == 2:
                return 'sp'

        else:
            if neighbor_count == 1:
                #Conjugated atoms already covered including carboxylate ion.
                #This instance is for single-bonded atoms with three lone pairs
                return 'unhybridized'
            elif neighbor_count >= 2:
                #Primarily covers the various single-bonded instances of oxygen and nitrogens
                return 'sp3'
            
        return 'unknown'

    def detect_functional_group_exception(self, elem1, elem2, bond_order, neighbors1, neighbors2,
                                        atom_id1=None, atom_id2=None, state=None):
        """
        Detect functional group exceptions and return override parameters.
        
        Returns dict with keys:
        - 'force_constant': override force constant (N/m) if applicable
        - 'functional_group': name of detected functional group
        - 'corrected_hybridization': corrected hybridization if applicable
        - 'note': explanatory note for user
        - 'frequency_range': expected IR range when available
        
        Returns empty dict if no exception applies.
        """
        # Due to the presence of neighboring atoms influences C-C and C-H signals,
        #   the following are series of exception detectors for special functional groups.
       
        exception = {}
        
        # ALDEHYDE C-H DETECTION
        # Aldehyde has: H bonded to sp2 C that is also bonded to O (C=O)
        # Pattern: C-H where C has O neighbor and is sp2 hybridized
        if (elem1 == 'C' and elem2 == 'H') or (elem1 == 'H' and elem2 == 'C'):
            # Identify which is the carbon and which is hydrogen
            if elem1 == 'C':
                c_neighbors = neighbors1
                h_neighbors = neighbors2
            else:
                c_neighbors = neighbors2
                h_neighbors = neighbors1
            
            # Check if carbon has oxygen neighbor
            neighbor_elements = [n[1] for n in c_neighbors] if c_neighbors else []
            has_oxygen = 'O' in neighbor_elements
            
            if has_oxygen and len(c_neighbors) == 2:
                # This is likely an aldehyde C-H
                # Aldehyde C-H appears at 2720-2820 cm⁻¹ (lower than alkyl C-H)
                # and involves coupling with C=O
                # Use sp2 hybridization constant for aldehyde
                exception = {
                    'force_constant': 540,  # Aldehyde C-H (sp2 carbon bonded to H)
                    'functional_group': 'aldehyde C-H',
                    'corrected_hybridization': 'C=sp2 (carbonyl carbon)',
                    'note': 'Aldehyde C-H detected (CHO group). Frequencies ~2720-2820 cm⁻¹.',
                    'frequency_range': '2720-2820 cm⁻¹'
                }
                print(f"✓ Functional group detected: ALDEHYDE C-H")
                print(f"  Override force constant: {exception['force_constant']} N/m")
                print(f"  Note: {exception['note']}")
            elif has_oxygen and len(c_neighbors) == 3:
                # This is likely an aldehyde C-H
                # Aldehyde C-H appears at 2720-2820 cm⁻¹ (lower than alkyl C-H)
                # and involves coupling with C=O
                # Use sp2 hybridization constant for aldehyde
                exception = {
                    'force_constant': 540,  # Aldehyde C-H (sp2 carbon bonded to H)
                    'functional_group': 'aldehyde C-H',
                    'corrected_hybridization': 'C=sp2 (carbonyl carbon)',
                    'note': 'Aldehyde C-H detected (CHO group). Frequencies ~2720-2820 cm⁻¹.',
                    'frequency_range': '2720-2820 cm⁻¹'
                }
                print(f"✓ Functional group detected: ALDEHYDE C-H")
                print(f"  Override force constant: {exception['force_constant']} N/m")
                print(f"  Note: {exception['note']}")
                            
        
        # CARBOXYLIC ACID O-H DETECTION
        # Carboxylic acid O-H: O bonded to C(=O)
        elif (elem1 == 'O' and elem2 == 'H') or (elem1 == 'H' and elem2 == 'O'):
            if elem1 == 'O':
                o_neighbors = neighbors1
            else:
                o_neighbors = neighbors2
            
            carbon_neighbors = [n for n in o_neighbors if n[1] == 'C']
            if carbon_neighbors and bond_order == 'single':
                for carbon in carbon_neighbors:
                    carbon_neighbors_of_carbon = self.get_neighbors_by_ID(carbon[3], radius=1.9)
                    oxygen_count = sum(1 for n in carbon_neighbors_of_carbon if n[1] == 'O')
                    if oxygen_count >= 2:
                        exception = {
                            'force_constant': 500,  # Carboxylic acid O-H (usually broad ~2500-3300 cm⁻¹)
                            'functional_group': 'carboxylic acid O-H',
                            'note': 'Carboxylic acid O-H detected (COOH group). Broad absorption ~2500-3300 cm⁻¹.',
                            'frequency_range': '2500-3300 cm⁻¹'
                        }
                        print(f"✓ Functional group detected: CARBOXYLIC ACID O-H")
                        print(f"  Override force constant: {exception['force_constant']} N/m")
                        print(f"  Note: {exception['note']}")
                        break

        # CONJUGATED CARBONYL C=O DETECTION
        elif (elem1 == 'C' and elem2 == 'O') or (elem1 == 'O' and elem2 == 'C'):
            if bond_order == 'double':
                # Double-check the actual C-O distance before treating as a carbonyl.
                # Some esters/anhydrides can have ambiguous metadata; use geometry
                # to avoid misclassifying C-O singles as C=O doubles.

                if elem1 == 'C':
                    c_neighbors = neighbors1
                    c_index = atom_id1
                else:
                    c_neighbors = neighbors2
                    c_index = atom_id2


                conjugated = False
                for neighbor in c_neighbors:
                    if neighbor[1] != 'C':
                        # Ignore non-carbon atoms adjacent to C in the Carbonyl
                        continue

                    neighbor_id = neighbor[3]
                    if elem1 == 'C':
                        carbonyl_id = atom_id1
                    else:
                        carbonyl_id = atom_id2

                    # Get the neighbors to the adjacent Carbon, excluding the Carbonyl Carbon
                    neighbor_neighbors = self.get_neighbors_by_ID(neighbor_id, radius=1.9, exclude_sel=f"id {carbonyl_id}")

                    for nn in neighbor_neighbors:
                        # Checking for a double Carbon-Carbon bond on adjacent carbon
                        if nn[1] != 'C':
                            continue
                        if nn[3] == neighbor_id:
                            continue

                        try:
                            dist_cc = cmd.get_distance(f"id {neighbor_id}", f"id {nn[3]}")
                        except Exception:
                            print(f"Unable to collect distance between Carbon Atoms ID: {neighbor_id} and {nn[3]}")
                            dist_cc = None

                        if dist_cc is not None and dist_cc < 1.47:
                            conjugated = True
                            break

                    if conjugated:
                        break

                if conjugated:
                    exception = {
                        'force_constant': 1050,  # Lower k for conjugated carbonyl
                        'functional_group': 'conjugated carbonyl C=O',
                        'note': 'Conjugated carbonyl detected. Conjugation lowers the C=O stretch frequency.',
                        'frequency_range': '1660-1715 cm⁻¹'
                    }
                    print(f"✓ Functional group detected: CONJUGATED CARBONYL C=O")
                    print(f"  Override force constant: {exception['force_constant']} N/m")
                    print(f"  Note: {exception['note']}")
        
        # HANDLE OXYGENS BONDED TO HYDROGEN ON C-O PAIRS
        # If oxygen has an H covalently bonded (distance <1.2 Å) and the C-O
        # distance is consistent with a single bond (>= 1.30 Å), treat as C-O single.
            else:
                if elem1 == 'O':
                    o_id = atom_id1
                    c_id = atom_id2
                    o_neighbors_local = neighbors1
                else:
                    o_id = atom_id2
                    c_id = atom_id1
                    o_neighbors_local = neighbors2

                h_neighbors = [n for n in o_neighbors_local if n[1] in ('H', 'D')]
                bonded_h = False
                for h in h_neighbors:
                    try:
                        d_oh = cmd.get_distance(f"id {o_id}", f"id {h[3]}")
                    except Exception:
                        print(f"Unable to collect distance between Hydrogen Atom ID: {h[3]} and Oxygen Atom ID: {o_id}")
                        d_oh = None
                    if d_oh is not None and d_oh < 1.20:
                        bonded_h = True
                        break

                if bonded_h:
                    try:
                        d_co = cmd.get_distance(f"id {c_id}", f"id {o_id}")
                    except Exception:
                        print(f"Unable to collect distance between Carbon Atom ID: {c_id} and Oxygen Atom ID: {o_id}") 
                        d_co = None

                    carbon_neighbors = self.get_neighbors_by_ID(c_id, radius=1.9, exclude_sel=f"id {o_id}")
                    other_oxygen_on_carbon = any(nn[1] == 'O' for nn in carbon_neighbors)

                    if d_co is None or d_co >= 1.30 or other_oxygen_on_carbon:
                        exception = {
                            'force_constant': BOND_DATA['C-O']['k'],
                            'functional_group': 'C-O single (O-H present)',
                            'note': 'Oxygen bonded to hydrogen; treating as C-O single bond (~1145 cm⁻¹).',
                            'frequency_range': '≈1145 cm⁻¹'
                        }
                        print(f"✓ Functional group detected: C-O single (oxygen bonded to H)")
                        print(f"  Override force constant: {exception['force_constant']} N/m")
                        print(f"  Note: {exception['note']}")

        return exception

    def get_force_constant_for_bond(self, elem1, elem2, bond_order, hybridization=None):
        """Map inferred bond order and hybridization to a bond constant from BOND_DATA."""
        # Force Constant searching based on bond order and hybridization.         
        if bond_order == 'aromatic':
            if {elem1, elem2} == {'C', 'C'}:
                return BOND_DATA['C-C']['k']
            bond_order = 'single'

        # C-H and C-D force constants depend on carbon hybridization.
        if bond_order == 'single' and set([elem1, elem2]) == {'C', 'H'}:
            hybrid = hybridization.strip().lower() if isinstance(hybridization, str) else hybridization
            if hybrid == 'sp3':
                return BOND_DATA['C-H']['k']
            if hybrid == 'sp2':
                return BOND_DATA['=C-H']['k']
            if hybrid == 'sp':
                return BOND_DATA['≡C-H']['k']
            return BOND_DATA['C-H']['k']

        if bond_order == 'single' and set([elem1, elem2]) == {'C', 'D'}:
            hybrid = hybridization.strip().lower() if isinstance(hybridization, str) else hybridization
            if hybrid == 'sp3':
                return BOND_DATA['C-D']['k']
            if hybrid == 'sp2':
                return BOND_DATA['=C-D']['k']
            if hybrid == 'sp':
                return BOND_DATA['≡C-D']['k']
            return BOND_DATA['C-D']['k']

        # All remaining element pairs are handled by element identity and bond order.
        key = None
        if bond_order == 'single':
            if f"{elem1}-{elem2}" in BOND_DATA:
                key = f"{elem1}-{elem2}"
            elif f"{elem2}-{elem1}" in BOND_DATA:
                key = f"{elem2}-{elem1}"
        elif bond_order == 'double':
            if f"{elem1}={elem2}" in BOND_DATA:
                key = f"{elem1}={elem2}"
            elif f"{elem2}={elem1}" in BOND_DATA:
                key = f"{elem2}={elem1}"
        elif bond_order == 'triple':
            if f"{elem1}≡{elem2}" in BOND_DATA:
                key = f"{elem1}≡{elem2}"
            elif f"{elem2}≡{elem1}" in BOND_DATA:
                key = f"{elem2}≡{elem1}"

        if key and key in BOND_DATA:
            return BOND_DATA[key]['k']

        for bond_name, data in BOND_DATA.items():
            if set(data['atoms']) == {elem1, elem2}:
                if bond_order == 'single' and data['multiplicity'] == 1:
                    return data['k']
                if bond_order == 'double' and data['multiplicity'] == 2:
                    return data['k']
                if bond_order == 'triple' and data['multiplicity'] == 3:
                    return data['k']
        return None

    def start_bond_picking(self):
        """Give instructions for bond picking"""
        self.pick_status_label.setText("👉 Click 2 atoms in PyMOL, then click 'Load Selected Bond'")
        
        cmd.deselect()
        
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        msg.setWindowTitle("Bond Selection Instructions")
        msg.setText("How to pick a bond:")
        msg.setInformativeText(
            "PyMOL uses atom PICKING (not selection) for this plugin.\n\n"
            "1. Simply CLICK on the first atom in the bond\n"
            "2. Then CLICK on the second atom\n"
            "   (Do NOT hold Shift - just plain clicks)\n"
            "3. The two atoms should turn green (picked indicator)\n"
            "4. Click 'Load Selected Bond' button\n\n"
            "The elements and masses will auto-fill.\n"
            "If bond multiplicity is ambiguous (e.g., C-C could\n"
            "be single, double, or triple), you will be asked to specify it."
        )
        msg.setWindowModality(Qt.WindowModal)
        msg.exec_()
    
    def stop_bond_picking(self):
        """Stop bond picking mode"""
        self.pick_status_label.setText("")
        cmd.deselect()
    
    def process_picked_bond(self):
        """Process the atoms selected in PyMOL"""
        try:
            # Get current state info
            current_state = cmd.get_state()
            if current_state > 1:
                print(f"Note: Reading from state {current_state} (NMR/trajectory)")
            
            # Get picked atoms using pk1 and pk2 (most recent picks)
            atom1 = []
            atom2 = []

            cmd.iterate("pk1", "atom1.append((index, elem, name, ID))", space={'atom1': atom1})
            cmd.iterate("pk2", "atom2.append((index, elem, name, ID))", space={'atom2': atom2})

            if (not atom1 or not atom2):
                QMessageBox.warning(self, "Picking Error", 
                                  "Please pick exactly 2 atoms by clicking them in PyMOL.\n\n"
                                  "PyMOL should be in picking mode.\n"
                                  "Simply click on one atom, then click on another.")
                return

            # Combine the two picked atoms
            atoms = atom1 + atom2
            
            # Get first two atoms and their internal PyMOL indices
            atom_idx1, elem1, name1, atom_id1 = atoms[0]
            atom_idx2, elem2, name2, atom_id2 = atoms[1]
            
            # Validate element detection
            if not elem1 or not elem2:
                QMessageBox.warning(self, "Element Detection Failed",
                                  f"Could not determine element types.\n"
                                  f"elem1={elem1}, elem2={elem2}\n"
                                  f"Check PDB file is properly formatted.")
                return
            
            # Check if atoms are bonded using distance
            try:
                dist = cmd.get_distance(f"id {atom_id1}", f"id {atom_id2}")
                # Rough covalent bond distance check (typical range 1.0-2.0 Å)
                if dist > 2.5:
                    response = QMessageBox.warning(self, "Distance Check",
                        f"Atoms are {dist:.2f}Å apart.\n"
                        f"Are they actually bonded?\n\n"
                        f"Click OK to continue anyway, Cancel to abort.",
                        QMessageBox.Ok | QMessageBox.Cancel)
                    if response != QMessageBox.Ok:
                        return
            except:
                # Distance calculation failed (maybe atoms in different models)
                print(f"Warning: Could not calculate distance between atoms {atom_id1} and {atom_id2}")
            
            print(f"\n{'='*60}")
            print(f"BOND PICKED: {elem1}({name1}) - {elem2}({name2})")
            print(f"Distance: {dist:.2f}Å" if 'dist' in locals() else "Distance: [could not calculate]")
            print(f"Atom IDs: {atom_id1}, {atom_id2}")
            print(f"{'='*60}\n")
            
            # Auto-fill the custom tab
            if elem1 in ATOMIC_MASSES and elem2 in ATOMIC_MASSES:
                # Switch to custom tab
                self.tabs.setCurrentIndex(1)
                
                # Set elements
                elem_idx1 = self.custom_elem1.findText(elem1)
                if elem_idx1 >= 0:
                    self.custom_elem1.setCurrentIndex(elem_idx1)
                
                elem_idx2 = self.custom_elem2.findText(elem2)
                if elem_idx2 >= 0:
                    self.custom_elem2.setCurrentIndex(elem_idx2)
                
                # Try to infer bond order and hybridization
                aromatic_flag = self.get_explicit_aromatic_flag()
                explicit_order = self.get_explicit_bond_order(current_state)
                if explicit_order is not None:
                    bond_order = {1: 'single', 2: 'double', 3: 'triple'}.get(int(explicit_order), 'single')
                    print(f"Explicit bond order found: {bond_order}")
                else:
                    bond_order = self.estimate_bond_order_by_length(elem1, elem2, dist if dist is not None else 0.0, aromatic_flag)
                    print(f"Estimated bond order: {bond_order}")
                if aromatic_flag:
                    print(f"Aromatic flag active for this bond")

                neighbors1 = self.get_neighbor_atoms('pk1', 'pk2')
                neighbors2 = self.get_neighbor_atoms('pk2', 'pk1')

                # NOTE: Do not force C-O double->single here; handle O-H bonded cases
                # inside detect_functional_group_exception so carbonyl detections
                # (including conjugation) are not preempted by a premature override.
                hybrid1 = self.infer_hybridization(elem1, len(neighbors1) + 1, aromatic_flag, bond_order)
                hybrid2 = self.infer_hybridization(elem2, len(neighbors2) + 1, aromatic_flag, bond_order)

                print(f"pk1 neighbors: {[n[1] for n in neighbors1]} (count={len(neighbors1)})")
                print(f"pk2 neighbors: {[n[1] for n in neighbors2]} (count={len(neighbors2)})")
                print(f"Inferred hybridization: {elem1}={hybrid1}, {elem2}={hybrid2}")

                # Change Bond_Orders
                if bond_order == 'single':
                    self.custom_bond.setCurrentIndex(0)
                elif bond_order == 'double' or bond_order == 'aromatic':
                    self.custom_bond.setCurrentIndex(1)
                elif bond_order == 'triple':
                    self.custom_bond.setCurrentIndex(2)

                if elem1 == 'C' and elem2 == 'C':
                    if hybrid1 == 'sp' or hybrid2 == 'sp':
                        self.custom_bond_adj.setCurrentIndex(3)
                    elif hybrid1 == 'sp2' or hybrid2 == 'sp2':
                        self.custom_bond_adj.setCurrentIndex(2)
                    elif hybrid1 == 'sp3' or hybrid2 == 'sp3':
                        self.custom_bond_adj.setCurrentIndex(1)

                elif elem1 == 'C':
                    if hybrid1 == 'sp3':
                        self.custom_bond_adj.setCurrentIndex(1)
                    elif hybrid1 == 'sp2':
                        self.custom_bond_adj.setCurrentIndex(2)
                    elif hybrid1 == 'sp':
                        self.custom_bond_adj.setCurrentIndex(3)

                elif elem2 == 'C':
                    if hybrid2 == 'sp3':
                        self.custom_bond_adj.setCurrentIndex(1)
                    elif hybrid2 == 'sp2':
                        self.custom_bond_adj.setCurrentIndex(2)
                    elif hybrid2 == 'sp':
                        self.custom_bond_adj.setCurrentIndex(3)
                else:
                    self.custom_bond_adj.setCurrentIndex(0)
            
                # Check for functional group exceptions (e.g., aldehyde C-H)
                functional_group_data = self.detect_functional_group_exception(
                    elem1, elem2, bond_order, neighbors1, neighbors2,
                    atom_id1, atom_id2, current_state
                )
                functional_group_note = ""
                
                if functional_group_data and 'force_constant' in functional_group_data:
                    suggested_k = functional_group_data['force_constant']
                    functional_group_name = functional_group_data.get('functional_group', 'special group')
                    functional_group_note = f"\n[{functional_group_name}] {functional_group_data.get('note', '')}"
                    if functional_group_data.get('corrected_hybridization'):
                        print(f"Corrected hybridization: {functional_group_data['corrected_hybridization']}")
                else:
                    suggested_k = self.get_force_constant_for_bond(elem1, elem2, bond_order, hybrid1 if elem1 == 'C' else hybrid2)
                
                self.custom_functional_group_data = functional_group_data or None
                if suggested_k is not None:
                    self.custom_k.setText(str(suggested_k))
                    print(f"Suggested force constant: {suggested_k} N/m for {bond_order} bond")
                else:
                    self.custom_functional_group_data = None
                    print(f"No default force constant available for {elem1}-{elem2} as a {bond_order} bond")
                    print("Please enter force constant manually")

                self.stop_bond_picking()
                QMessageBox.information(self, "Bond Loaded", 
                                        f"Bond {elem1}-{elem2} loaded!\n"
                                        f"Inferred bond order: {bond_order}\n"
                                        f"Inferred hybridization: {elem1}={hybrid1}, {elem2}={hybrid2}\n\n"
                                        f"Elements and masses auto-filled.\n"
                                        f"{'Force constant suggested.' if suggested_k else 'Please enter force constant.'}"
                                        f"{functional_group_note}\n\n"
                                        f"Note: Uses harmonic oscillator approximation.")
                return
            else:
                QMessageBox.warning(self, "Unsupported Element", 
                                    f"One or both elements not supported: {elem1}, {elem2}")
                self.stop_bond_picking()
                return
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error processing bond: {str(e)}")
            self.stop_bond_picking()
            return

    def calculate_common_bond(self):
        """Calculate from selected common bond"""
        try:
            current_row = self.bond_list.currentRow()
            if current_row < 0:
                QMessageBox.warning(self, "No Selection", "Please select a bond from the list")
                return
            
            bond_name = list(BOND_DATA.keys())[current_row]
            bond_data = BOND_DATA[bond_name]
            
            atom1, atom2 = bond_data['atoms']
            mass1 = ATOMIC_MASSES[atom1]
            mass2 = ATOMIC_MASSES[atom2]
            force_constant = bond_data['k']
            
            reduced_mass = calculate_reduced_mass(mass1, mass2)
            wavenumber = calculate_wavenumber(force_constant, reduced_mass)
            
            results = f"""
{'='*60}
Results for {bond_name} bond: {bond_data['description']}
{'='*60}
Atoms: {atom1} ({mass1:.3f} amu) - {atom2} ({mass2:.3f} amu)
Reduced mass (μ) = {reduced_mass:.6f} amu
Force constant (k) = {force_constant} N/m
Predicted wavenumber (ν̃) = {wavenumber:.2f} cm⁻¹
{'='*60}
"""
            self.common_results.setText(results)
            
        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))
    
    def calculate_custom(self):
        """Calculate from custom inputs"""
        try:
            elem1 = self.custom_elem1.currentText()
            elem2 = self.custom_elem2.currentText()
            mass1 = float(self.custom_mass1.text())
            mass2 = float(self.custom_mass2.text())
            force_constant = float(self.custom_k.text())
            
            reduced_mass = calculate_reduced_mass(mass1, mass2)
            functional_group_data = self.custom_functional_group_data
            
            if functional_group_data and functional_group_data.get('frequency_range'):
                results = f"""
{'='*60}
Results for {elem1}-{elem2} bond
{'='*60}
Element 1: {elem1} ({mass1:.3f} amu)
Element 2: {elem2} ({mass2:.3f} amu)
Reduced mass (μ) = {reduced_mass:.6f} amu
Force constant (k) = {force_constant:.2f} N/m

Functional group override: {functional_group_data.get('functional_group')}
Note: {functional_group_data.get('note')}
Predicted IR frequency range: {functional_group_data.get('frequency_range')}
{'='*60}
"""
            else:
                wavenumber = calculate_wavenumber(force_constant, reduced_mass)
                results = f"""
{'='*60}
Results for {elem1}-{elem2} bond
{'='*60}
Element 1: {elem1} ({mass1:.3f} amu)
Element 2: {elem2} ({mass2:.3f} amu)
Reduced mass (μ) = {reduced_mass:.6f} amu
Force constant (k) = {force_constant:.2f} N/m

*** Predicted IR Wavenumber (ν̃) = {wavenumber:.2f} cm⁻¹ ***
{'='*60}
"""
            self.custom_results.setText(results)
            
        except ValueError:
            QMessageBox.critical(self, "Input Error", "Please enter valid numerical values")
        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))
    
    def calculate_reverse(self):
        """Calculate force constant from wavenumber"""
        try:
            mass1 = float(self.reverse_mass1.text())
            mass2 = float(self.reverse_mass2.text())
            wavenumber = float(self.reverse_wn.text())
            
            reduced_mass = calculate_reduced_mass(mass1, mass2)
            force_constant = calculate_force_constant(wavenumber, reduced_mass)
            
            results = f"""
{'='*60}
Results
{'='*60}
Mass 1 = {mass1:.3f} amu
Mass 2 = {mass2:.3f} amu
Reduced mass (μ) = {reduced_mass:.6f} amu
Wavenumber (ν̃) = {wavenumber:.2f} cm⁻¹
Force constant (k) = {force_constant:.2f} N/m
{'='*60}
"""
            self.reverse_results.setText(results)
            
        except ValueError:
            QMessageBox.critical(self, "Input Error", "Please enter valid numerical values")
        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))


def __init_plugin__(app=None):
    """
    PyMOL plugin initialization function.
    """
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('IR Calculator', run_plugin_gui)


def run_plugin_gui():
    """Launch the plugin GUI"""
    global plugin_dialog
    plugin_dialog = IRCalculatorDialog()
    plugin_dialog.show()


# For standalone testing
if __name__ == "__main__":
    import sys
    from PyQt5.QtWidgets import QApplication
    
    app = QApplication(sys.argv)
    dialog = IRCalculatorDialog()
    dialog.show()
    sys.exit(app.exec_())
