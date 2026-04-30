# Electroporation-Action-Potential-Modeling
# Modelling the Effect of Electroporation on Cardiomyocyte Action Potentials  
## Project Overview

This repository contains MATLAB code and simulation results from my Bachelor's thesis:

**Modelling the Effect of Electroporation on the Triggering of Action Potentials in Cardiomyocytes**

The project investigates how electroporation affects the electrophysiological response of cardiomyocytes, with a focus on action potential generation and intracellular calcium dynamics.

## Background

Electroporation is a process in which high-voltage electric pulses increase the permeability of the cell membrane. In cardiac electrophysiology, irreversible electroporation is used as a non-thermal ablation method for treating cardiac arrhythmias.

However, around irreversibly electroporated tissue, there is also a region of reversibly electroporated tissue. In this region, cardiomyocytes survive the electric pulse exposure, but their ability to generate action potentials may be affected.

## Thesis Goal

The goal of this project was to model the effect of electroporation on isolated cardiomyocytes and analyze how increased membrane permeability influences:

- action potential generation
- repolarization
- intracellular calcium concentration
- electrical stability of cardiomyocytes

## What I Implemented

The project uses two mathematical models of cardiomyocyte electrophysiology:

1. **Luo-Rudy model**
   - Used to simulate ventricular cardiac action potentials
   - Extended by adding an electroporation current

2. **Pandit et al. model**
   - Used to simulate rat ventricular cardiomyocytes
   - Includes intracellular calcium handling
   - Extended by adding ion-specific electroporation currents

Both models were implemented and simulated in MATLAB using the `ode15s` solver.

## Electroporation Model

Electroporation was modelled as an increase in membrane permeability caused by the formation of pores in the cell membrane.

The number of pores was varied in simulations to analyze how different levels of electroporation affect cardiomyocyte behavior.

## Key Results

The simulations showed that even a small number of membrane pores can significantly affect cardiomyocyte electrophysiology.

Main observations:

- action potential duration increases with the number of pores
- repolarization becomes slower
- at higher pore numbers, cells fail to repolarize
- cardiomyocytes may remain depolarized
- intracellular calcium concentration increases significantly
- action potential generation can be completely suppressed

These results suggest that even weak electroporation may strongly influence the electrical activity of cardiomyocytes.

## Repository Structure

```text
Electroporation-Cardiomyocyte-Model/
│
├── project-code/
│   ├── luo-rudy/
│   │   └── MATLAB code for the Luo-Rudy model
│   │
│   └── pandit/
│       └── MATLAB code for the Pandit et al. model
│
├── results/
│   ├── luo-rudy/
│   │   └── simulation results and figures
│   │
│   └── pandit/
│       └── simulation results and figures
│
├── docs/
│   └── thesis PDF or additional documentation
│
├── README.md
├── LICENSE
└── .gitignore
## Repository Structure

```text
Electroporation-Cardiomyocyte-Model/
│
├── project-code/
│   ├── luo-rudy/
│   │   └── MATLAB code for the Luo-Rudy model
│   │
│   └── pandit/
│       └── MATLAB code for the Pandit et al. model
│
├── results/
│   ├── luo-rudy/
│   │   └── simulation results and figures
│   │
│   └── pandit/
│       └── simulation results and figures
│
├── docs/
│   └── thesis PDF or additional documentation
│
├── README.md
├── LICENSE
└── .gitignore

## Repository Structure

```text
Electroporation-Cardiomyocyte-Model/
│
├── project-code/
│   ├── luo-rudy/
│   │   └── MATLAB code for the Luo-Rudy model
│   │
│   └── pandit/
│       └── MATLAB code for the Pandit et al. model
│
├── results/
│   ├── luo-rudy/
│   │   └── simulation results and figures
│   │
│   └── pandit/
│       └── simulation results and figures
│
├── docs/
│   └── thesis PDF or additional documentation
│
├── README.md
├── LICENSE
└── .gitignore
