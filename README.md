# IR Wavenumber Calculator Plugin for PyMOL

[![DOI](https://zenodo.org/badge/1415915394.svg)](https://doi.org/10.5281/zenodo.23324150)

An educational PyMOL plugin that calculates approximate IR stretching wavenumbers, reduced masses, and force constants using the harmonic oscillator model, with direct bond selection from 3D structures in PyMOL.

> This tool is intended for teaching. Predicted values are harmonic-oscillator estimates based on approximate force constants, not substitutes for experimental spectra or quantum-chemical frequency calculations.

## Features

The plugin opens a non-modal window with three tabs:

- **Common Bonds**: pick a bond type from a built-in table of approximate force constants and calculate its predicted wavenumber.
- **Custom Calculation**: click two atoms in PyMOL (or enter elements manually). The plugin reads element identities, estimates bond order and hybridization, applies a few functional-group adjustments (e.g. aldehyde C–H, conjugated C=O, carboxylic acid O–H), and suggests a force constant. Isotopic masses (D, ¹³C, ¹⁵N, ¹⁸O, etc.) are supported.
- **Find Force Constant**: enter an experimental wavenumber and two atomic masses to back-calculate an effective force constant.

## Requirements

- PyMOL 3.x (incentive or open-source) with its bundled PyQt5
- Tested on Ubuntu 24.04 and Windows 10/11. Other platforms and PyMOL versions have not been systematically tested.

## Installation

1. Download `ir_calculator.py` from this repository (or the latest [release](https://github.com/dperygin/IR_Calculator_PyMOL/releases)).
2. In PyMOL: **Plugin → Plugin Manager → Install New Plugin → Choose file…** and select `ir_calculator.py`.
3. Restart PyMOL.
4. Open it from **Plugin → IR Calculator**.


## Quick start: calculate a bond from a structure

1. Load or build a molecule in PyMOL (e.g. `fragment ethane`, or load a `.mol`/`.sdf` with explicit bond orders).
2. Open **Plugin → IR Calculator → Custom Calculation**.
3. Click **📍 Select Bond in PyMOL**, then click the two atoms of the bond in the viewer.
4. Click **⬇️ Load Selected Bond**. Elements, masses, bond order, and a suggested force constant are filled in.
5. Click **Calculate**.

Diagnostic messages (detected bond order, hybridization, functional-group overrides) are printed to the PyMOL console.

## Method

The predicted wavenumber is

ν̃ = (1 / 2πc) · √(k / μ),  μ = m₁m₂ / (m₁ + m₂)

with c = 2.99792458 × 10¹⁰ cm/s and 1 amu = 1.66053906660 × 10⁻²⁷ kg. Masses are monoisotopic (¹H, ¹²C, ¹⁶O, …). Force constants are approximate empirical values for idealized bonds; see the paper for the full table and a discussion of limitations (anharmonicity, vibrational coupling, hydrogen bonding, environment effects).

## Teaching materials

The `activities/` folder contains three student handouts described in the paper:

1. Exploring Bond Strength and Vibrational Frequency
2. The Deuterium Effect
3. Structure Determination from IR Data

Instructors who would like answer keys can contact the authors (see below).

## Citation

If you use this plugin in teaching or research, please cite:

> Wydemon, E.; Perygin, D. H.; Ogden, P. B. Development and Implementation of an Interactive IR Wavenumber Calculator Plugin for PyMOL: Enhancing Student Understanding of Vibrational Spectroscopy Through 3D Molecular Visualization. Manuscript in preparation.

To cite the software itself, please cite the archived release:

> Wydemon, E.; Perygin, D. H.; Ogden, P. B. *IR Wavenumber Calculator Plugin for PyMOL*, version 1.0.1; Zenodo, 2026. https://doi.org/10.5281/zenodo.23325492

The DOI [10.5281/zenodo.23324150](https://doi.org/10.5281/zenodo.23324150) always resolves to the latest version. Full citation metadata is in [`CITATION.cff`](CITATION.cff).

## License

Code: MIT License (see [`LICENSE`](LICENSE)). Teaching materials in `activities/`: [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/). You may share and adapt them for any purpose, including commercial use, with attribution to the authors.

## Contact

Phillip Ogden, pogden@jsu.edu. Bug reports and suggestions are welcome via [GitHub Issues](https://github.com/dperygin/IR_Calculator_PyMOL/issues).

## Acknowledgments

PyMOL is a trademark of Schrödinger, LLC. This plugin is an independent project and is not affiliated with or endorsed by Schrödinger.
