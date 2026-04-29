# TransCal

Finite-element study repository for incompressible 2D flow problems, organized as coursework deliverables with mesh generation, simulation scripts, and exported VTK results.

## Overview

This repository contains two main work packages:

- **mesh-generation**: educational mesh generation utilities for 1D/2D structured grids.
- **flow-simulations**: flow simulations in streamfunction-vorticity form over different geometries (channel, cavity, cylinder, and step), plus generated meshes and time-dependent VTK outputs.

## Repository Structure

- `mesh-generation/gerador_de_malhas.py`: structured mesh generator (1D/2D, triangular or quadrilateral layout options).
- `flow-simulations/shared/correnteVorticidade2D.py`: Eulerian streamfunction-vorticity solver setup.
- `flow-simulations/shared/correnteVorticidade2D-lagrangeano.py`: Lagrangian/Eulerian variant for the same formulation.
- `flow-simulations/*/*.geo`: geometry definitions (Gmsh).
- `flow-simulations/*/*.msh`: generated meshes.
- `flow-simulations/*/sol-*/*.vtk`: simulation snapshots for post-processing.

## Requirements

- Python 3.10+
- Gmsh (for `.geo` to `.msh` workflows)
- Python packages:
  - `numpy`
  - `matplotlib`
  - `meshio`

Suggested installation:

```bash
pip install numpy matplotlib meshio
```

## Quick Start

Run from the repository root.

1. Generate/inspect structured meshes:

```bash
python mesh-generation/gerador_de_malhas.py
```

2. Run the Eulerian flow solver:

```bash
python flow-simulations/shared/correnteVorticidade2D.py
```

3. Run the Lagrangian/Eulerian variant:

```bash
python flow-simulations/shared/correnteVorticidade2D-lagrangeano.py
```

## Outputs

- Intermediate and transient fields are exported as `.vtk` files (`solucao-*.vtk`) for visualization in tools such as ParaView.

## Notes

- Scripts are currently configured with hardcoded paths and case parameters for coursework scenarios.
- If you move folders, update mesh input/output paths directly in the solver scripts.
