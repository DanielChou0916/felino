---
title: "Felino: Extension of an open-source phase-field framework to geomaterial fracture"
tags:
  - C++
  - Finite Element Analysis
  - Phase Field Fracture Mechanics
  - Geomaterials
authors:
  - name: Daniel T. Chou
    orcid: 0000-0001-7358-8819
    affiliation: 1
affiliations:
  - name: Georgia Institute of Technology, United States
    index: 1
date: 28 July 2025
bibliography: paper.bib
---


# Summary
This work presents an extended version of Felino [@chou2025felino], an open-source phase-field fracture framework.
Felino is implemented as an application built on top of the MOOSE finite-element framework [@permann2020moose], which uses libMesh as its underlying numerical library [@kirk2006libmesh].
The extension introduces constitutive models for geomaterials, enabling simulations where tensile and compressive strengths differ. 
- The official website of felino: [Felino official website](https://danielchou0916.github.io/felino.github.io/#).
- Benchmark example relevant to this extension :[Uniaxial Compression on Composite Material](https://danielchou0916.github.io/felino.github.io/tutorials/4_composite_uc2D/)

# Statement of need
Phase-field fracture models have been widely used in metallic fatigue simulations. However, geomaterials exhibit asymmetric mechanical behavior, especially under compressive-shear loading, which requires more advanced energy-splitting formulations. This updated version implements three constitutive models that capture this asymmetry:

1. **Representative Crack Element (RCE)** – interpolates between intact and cracked states using strain jump projections[@storm2020concept].
2. **Drucker–Prager Decomposition** – derives activated and inactivated energy parts from a pressure-dependent failure criterion[@navidtehrani2022general].
3. **Extra Driving Force Formulation** – introduces an additional compressive-shear resistance term to the phase-field equation[@liu2025emergence].

Details of each model: [Tension-Compression Asymmetry](https://danielchou0916.github.io/felino.github.io/technical_contents/decomposition/)
Details of programming objects: [(AD)LinearElasticPFFractureStress and (AD)ComputePFFStress](https://danielchou0916.github.io/felino.github.io/feature_objects/crack_stress/)
# Key features
- Extended constitutive models for geomaterials.
- Support for asymmetric tensile/compressive fracture behavior.
- Benchmark examples for uniaxial compression tests.
- Fully integrated with MOOSE automatic differentiation.


# References