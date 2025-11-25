---
title: 'dKMC: Delocalised kinetic Monte Carlo for simulating fundamental transport processes involving partially delocalised carriers in disordered materials'
tags:
  - Julia
  - organic semiconductors
  - organic photovoltaics
  - charge transport
  - charge separation
  - materials science
  - disordered materials
authors:
  - name: Daniel Balzer
    orcid: 0000-0002-4640-1059
    affiliation: 1
  - name: Ivan Kassal
    orcid: 0000-0002-8376-0819
    affiliation: 1
    corresponding: true
affiliations:
 - name: School of Chemistry, University of Sydney, NSW 2006, Australia
   index: 1
date: 28 April 2025
bibliography: paper.bib
---

# Summary
The movement of charge and energy is a fundamental process in materials science, underpinning technologies such as solar cells, light-emitting diodes, batteries, and electronics. Transport is well understood in both highly ordered materials (band conduction) and highly disordered ones (hopping conduction). However, in moderately disordered materials—including many organic semiconductors—transport lies in the intermediate transport regime between these well-understood extremes. Accurately modelling intermediate-regime conduction is difficult because describing wavefunction delocalisation requires a fully quantum-mechanical treatment, which is challenging in disordered materials that lack periodicity. We describe delocalised kinetic Monte Carlo (`dKMC`), the first theoretical approach to treat, in three dimensions, all the processes crucial in organic semiconductors: disorder, delocalisation, and polaron formation. As a result, it can treat the intermediate transport regime between band and hopping conduction. `dKMC` reveals that the fundamental physics of transport in moderately disordered materials is that of charges and excitons hopping between partially delocalised electronic states. In this work, we release the `dKMC.jl` package, which contains modules for simulating the fundamental processes of charge and exciton transport as well as charge separation and generation.

![Figure 1: While transport mechanisms are well understood in the extremes, coherent band conduction through extended states and incoherent hopping through localised states, they remain poorly understood in the intermediate regime where many organic semiconductors lie. Figure adapted with permission from [@Balzer2021].](intermediate_regime.png)

# Statement of need
The need for dKMC is demonstrated by its ability to solve two important problems. 

First, `dKMC` is the first computational technique that is able to include all of the processes crucial in organic semiconductors while remaining computationally tractable enough to treat realistic, three-dimensional systems on mesoscopic time and length scales [@Balzer2021]. Before `dKMC`, the difficulty in modelling transport in the intermediate regime had prevented the development of such a theory. Instead, prior approaches either needed to exclude one of the key ingredients (disorder, delocalisation, and polaron formation), which could lead to inaccurate results, or included them all but restricted the application to smaller systems (dimension, time, and length scales). `dKMC` solves these problems, striking a balance between the level of approximation and the size of the system it can simulate. Therefore,`dKMC` provides a simulation tool that can both accurately and efficiently simulate fundamental transport processes involving partially delocalised carriers.

Second, `dKMC` explains the often confusing behaviour of organic electronics, including organic photovoltaics (OPVs). Most models of transport in disordered organic semiconductors assume hopping transport, where charge carriers or excitons are localised onto individual molecules and move via thermally assisted hops from one molecule to the next. However, hopping transport fails to explain how charges and excitons move as fast as they do, or how charges in OPVs overcome their strong Coulomb attraction and separate from CT states as efficiently as they do, or how charges are generated so efficiently even in OPV devices with little to no energetic offsets. Hopping transport fails because, in many organic semiconductors, the charges and excitons remain delocalised across multiple molecules. By going beyond hopping, `dKMC` reveals that delocalisation improves each of the four fundamental transport processes in OPVs: charge transport [@Balzer2021], exciton transport [@Balzer2023], charge separation [@Balzer2022], and charge generation [@Balzer2024]. Delocalisation improves all of these important transport processes in essentially the same way, by enabling carriers to hop further and faster, explaining the failure of classical theories. Therefore, `dKMC` is a tool for explaining the otherwise unpredictable behaviour observed in devices.

`dKMC` was developed for organic semiconductors, but it can be more generally applied to other materials that lie in the intermediate regime. In particular, the first application of `dKMC` was to OPVs and therefore the package contains a module for each of the four fundamental processes in OPVs: charge transport, exciton transport, charge separation, and charge generation. However, the transport modules can be more generally applied to organic semiconductor materials used in other devices, such as organic light-emitting diodes or organic field-effect transistors.

# Acknowledgements
We were supported by the Westpac Scholars Trust (Research Fellowship and Future Leaders Scholarship), by the Australian Research Council (DP220103584), by the Australian Government Research Training Program, and by the University of Sydney Nano Institute Grand Challenge Computational Materials Discovery. We were also supported by computational resources from the National Computational Infrastructure (Gadi) and by the University of Sydney Informatics Hub (Artemis). 

# References
