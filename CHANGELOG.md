# Changelog

All notable changes to this project are documented here. This project uses [Semantic Versioning](https://semver.org/).

## [1.0.1] - 2026-10-12

### Fixed
- C–D force constants now equal the corresponding C–H values (sp³ 477, sp² 500, sp 597 N/m; previously 500, 530, 590). Isotopic substitution does not change the force constant, so predicted C–H/C–D shifts now arise entirely from the change in reduced mass, consistent with the isotope-effect discussion in the paper and Activity 2.
- Updated the approximate wavenumbers shown in the C–D bond descriptions (~2170, ~2220, ~2420 cm⁻¹).

## [1.0.0] - 2026-10-07

Initial public release; the version described in the accompanying paper.

### Added
- Common Bonds, Custom Calculation, and Find Force Constant tabs
- Direct two-atom bond selection from PyMOL (pk1/pk2)
- Bond-order detection from explicit PyMOL bond orders, with a distance-based fallback
- Hybridization-aware C–H / C–D force constants
- Functional-group adjustments: aldehyde C–H, carboxylic acid O–H, conjugated C=O, C–O single bonds
- Isotope-specific masses for H, C, N, O, F, P, S, Cl, Br, I, Si, B
