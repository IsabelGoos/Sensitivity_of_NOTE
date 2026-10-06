# Neutrino Oscillation Tomography of the Earth

This repository contains the code associated with the publication [Sensitivity of neutrino oscillations to the Earth’s interior properties](https://www.sciencedirect.com/science/article/pii/S0012821X26005340?ref=pdf_download&fr=RR-2&rr=a45b4c3c4c340477), with the goal of making the analysis and results fully reproducible.

## Overview

Understanding the Earth's deep interior is difficult because seismic observations alone can be ambiguous: different combinations of **temperature, composition, and mass density** can produce similar seismic signatures. In this project, I investigate whether **atmospheric neutrinos** can provide an independent constraint on the Earth's interior properties.

Atmospheric neutrinos are produced by the constant flux of cosmic rays entering our atmosphere and interacting with air molecules. Neutrinos come in 3 types, ν<sub>e</sub>, ν<sub>μ</sub>, and ν<sub>τ</sub>, together with their corresponding antiparticles. As they propagate, neutrinos have a probability to change from one type to another, a quantum-mechanical phenomenon called **neutrino oscillations**. As neutrinos travel through the Earth, their oscillation behaviour is affected by the **electron density** of the material they cross. Because electron density depends on both mass density and chemical composition, neutrino measurements provide a complementary probe of the Earth's interior to conventional seismic observations.

## Analysis pipeline

This framework provides quantitative forecasts of the performance of next-generation neutrino detectors in probing the electron density of the Earth's deep interior.

The framework:

* **Forecasts the sensitivity** of the next-generation neutrino detectors **KM3NeT/ORCA, Hyper-Kamiokande, and DUNE**.
* Compares realistic detector performance with that of an idealized detector, which represents the intrinsic **sensitivity limit** of neutrino oscillation tomography.
* **Computes Fréchet derivatives**, which are required for a future joint inversion pipeline combining neutrino and seismic data.
* **Tests the linearity of the detector response** to electron-density perturbations, providing insight into the validity of a linearized inversion framework and helping inform the choice of inversion methodology.

## Results

The simulations show that an idealized detector has its strongest sensitivity to perturbations in the Earth's core, while for realistic detector configurations the sensitivity shifts toward the mantle, including the mantle transition zone. In particular, variations in electron density associated with changes in composition, such as hydrogen enrichment, may produce potentially measurable signatures.

These results demonstrate the potential of neutrino oscillations as a complementary probe of the Earth's electron-density profile and provide a foundation for future joint constraints from **neutrino and seismic observations**


