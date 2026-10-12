# Changelog

All notable changes to this project are documented here. This project uses [Semantic Versioning](https://semver.org/).

## [1.0.0] - 2026-10-07

Initial public release; the version described in the accompanying paper.

### Added
- Common Bonds, Custom Calculation, and Find Force Constant tabs
- Direct two-atom bond selection from PyMOL (pk1/pk2)
- Bond-order detection from explicit PyMOL bond orders, with a distance-based fallback
- Hybridization-aware C–H / C–D force constants
- Functional-group adjustments: aldehyde C–H, carboxylic acid O–H, conjugated C=O, C–O single bonds
- Isotope-specific masses for H, C, N, O, F, P, S, Cl, Br, I, Si, B
