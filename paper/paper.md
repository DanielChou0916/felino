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

- Installation instruction: [README](https://github.com/DanielChou0916/felino/blob/main/README.md) or [Felino official website](https://danielchou0916.github.io/felino.github.io/installation/)

- Official website of Felino: [Felino official website](https://danielchou0916.github.io/felino.github.io/#)

- Benchmark example (this extension): [Uniaxial Compression on Composite Material](https://danielchou0916.github.io/felino.github.io/tutorials/4_composite_uc2D/)


# Statement of need
Phase-field fracture models have been increasingly used to study fatigue and fracture in metals and other solids due to their ability to represent complex crack patterns without explicit tracking [@li2023review] [@kalina2023overview]. However, geomaterials often exhibit strong tension–compression asymmetry and pressure-dependent failure, especially under compressive–shear loading, which motivates specialized energy-splitting and crack-driving-force formulations. This release extends Felino with three constitutive options targeting such asymmetry:

1. **Representative Crack Element (RCE)** – interpolates between intact and cracked states using strain jump projections [@storm2020concept].
2. **Drucker–Prager Decomposition** – derives activated and inactivated energy parts from a pressure-dependent failure criterion [@navidtehrani2022general].
3. **Extra Driving Force Formulation** – introduces an additional compressive–shear resistance term in the phase-field evolution equation [@liu2025emergence].

Details of each model: [Tension-Compression Asymmetry](https://danielchou0916.github.io/felino.github.io/technical_contents/decomposition/)
Details of programming objects: [(AD)LinearElasticPFFractureStress and (AD)ComputePFFStress](https://danielchou0916.github.io/felino.github.io/feature_objects/crack_stress/)


# State of the field
Phase-field fracture simulations are commonly performed either by (i) implementing custom weak forms in open-source finite-element toolchains (e.g., FEniCS) [@BarattaEtal2023; @ScroggsEtal2022; @BasixJoss], or (ii) using commercial multiphysics/FEA platforms where phase-field formulations are realized via user-defined PDE interfaces and/or user subroutines and user elements (e.g., COMSOL and Abaqus-based implementations) [@COMSOL; @ABSstd2009; @molnar20172d; @molnar2020open]. In parallel, the MOOSE ecosystem provides reusable infrastructure for phase-field modeling and multiphysics coupling [@permann2020moose; @kirk2006libmesh].

Felino is positioned as a MOOSE-based, reproducible workflow for geomaterial fracture, focusing on tension–compression asymmetry and pressure-dependent constitutive behavior. Rather than proposing a new solver stack, Felino “builds” on the MOOSE application model to deliver (i) modular constitutive formulations implemented as reusable objects, (ii) automatic-differentiation-based residual/Jacobian consistency, and (iii) curated benchmark and validation examples that lower the barrier for both research prototyping and engineering-oriented studies.

# Software design
Felino is implemented as a MOOSE application, leveraging MOOSE’s modular kernels/material system and automatic differentiation to provide consistent residual and Jacobian evaluations for coupled mechanics–phase-field problems [@permann2020moose]. Geomaterial-specific constitutive options are exposed through input-file configuration, enabling users to switch among formulations without rewriting governing equations in code. This design prioritizes extensibility (adding new splits/driving forces as reusable objects) and reproducibility (versioned examples and documented workflows) over specialization to a narrow set of built-in material models.

# Limitations
At present, Felino focuses on (nonlinear) elasticity-based formulations and does not yet provide plasticity models. Large-scale performance benchmarking against commercial codes is also not included in this work.


# Key features

- Extended constitutive models for geomaterials.

- Support for asymmetric tensile/compressive fracture behavior.

- Benchmark examples for uniaxial compression tests.

- Fully integrated with MOOSE automatic differentiation.

# Research impact statement
Felino provides documented, executable benchmark cases and validation-style examples for geomaterial phase-field fracture, including uniaxial compression on composite/heterogeneous configurations. These examples are intended to serve as reproducible starting points for researchers developing new constitutive splits or crack-driving-force terms, and for engineers seeking transparent, scriptable workflows for parametric studies. By packaging the geomaterial formulations as reusable objects within the MOOSE ecosystem, Felino supports reuse in broader multiphysics settings and facilitates extension to additional constitutive families in future work.

# Benchmark examples
This section provides cookbook-style example cases that demonstrate how Felino can be applied to representative fracture and fatigue scenarios. Each example page includes a brief description of the purpose, required inputs, expected outputs, and a quick-start command to run the case.

1. [Single Edge-Notch Tension (SENT) fatigue analysis](https://danielchou0916.github.io/felino.github.io/tutorials/1_sent_cyclic_pulling/)
2. [Fatigue behavior under thermal shock](https://danielchou0916.github.io/felino.github.io/tutorials/2_thermal_shock_metal/)
3. [Asymmetric three-point bending fatigue test](https://danielchou0916.github.io/felino.github.io/tutorials/3_3D_anisotropic_bending/)
4. [Uniaxial compression on composite material](https://danielchou0916.github.io/felino.github.io/tutorials/4_composite_uc2D/)



# AI usage disclosure
Generative AI tools were used to assist with English language editing of the manuscript. All technical content, claims, and citations were written and verified by the authors.

# References