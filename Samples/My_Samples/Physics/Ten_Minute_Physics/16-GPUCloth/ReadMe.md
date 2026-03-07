# GPU Cloth Simulation

![Sample Preview](sample.gif)

Implementation of lesson [16 - GPU Cloth](https://github.com/matthias-research/pages/blob/master/tenMinutePhysics/16-GPUCloth.py) from the [Ten Minute Physics](https://matthias-research.github.io/pages/tenMinutePhysics/) series by [Matthias Müller](https://github.com/matthias-research).

**Created by:** [PAMinerva](https://github.com/PAMinerva)<br>
**Credits to:** [Matthias Müller](https://github.com/matthias-research) for the original concept and implementation.<br>
**Powered by:** [Wicked Engine](https://github.com/turanszkij/WickedEngine)<br>

## Overview

Real-time GPU cloth simulation using position-based dynamics with two solver modes:
- **Coloring (hybrid Gauss-Seidel / Jacobi)**: Independent constraint passes use direct writes; dependent passes use CAS-loop float atomic accumulation.
- **Full Jacobi**: All constraints solved simultaneously using CAS-loop float atomic accumulation.

## Features

- **GPU Compute Simulation** — All physics runs on the GPU via compute shaders
- **CAS-Loop Float Atomics** — Concurrent accumulation using `InterlockedCompareExchange` with `asuint` and `asfloat`
- **Two Solver Modes** — Toggle between coloring (faster convergence) and full Jacobi at runtime
- **5-Pass Constraint Coloring** — 4 independent stretch passes + 1 Jacobi pass for shear/bending
- **Sphere & Ground Collision** — With friction
- **Interactive Grabbing** — GPU raycast (Möller-Trumbore) for picking and dragging cloth vertices

## Controls

- **WASD / Arrow Keys** — Move camera
- **Right Mouse / Middle Mouse** — Rotate camera
- **Left Mouse** — Pick and drag a vertex
- **GUI Controls**:
  - **Run/Stop** — Start or pause the simulation
  - **Restart** — Reset the simulation to initial state
  - **Coloring / Jacobi** — Toggle solver mode
