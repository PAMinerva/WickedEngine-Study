# Soft Body Skinning Sample with Wavefront OBJ Model Loading

Implementation of lesson [12 - Soft Body Skinning](https://github.com/matthias-research/pages/blob/master/tenMinutePhysics/12-softBodySkinning.html) from the [Ten Minute Physics](https://matthias-research.github.io/pages/tenMinutePhysics/) series by [Matthias Müller](https://github.com/matthias-research).

**Created by:** [PAMinerva](https://github.com/PAMinerva)<br><br>
**Credits to:** [Matthias Müller](https://github.com/matthias-research) for the original concept and implementation.<br><br>
**Powered by:** [Wicked Engine](https://github.com/turanszkij/WickedEngine)<br>
Special thanks to [Turánszki János](https://github.com/turanszkij) for creating Wicked Engine and making it available under the MIT license.

## Overview

This project demonstrates real-time soft body simulation and skinning using tetrahedral meshes and barycentric coordinates, implemented with the Wicked Engine C++ API.<br>
**OBJ Model Loading:** This implementation uses the ModelImporter_OBJ module from the Wicked Engine editor to load OBJ files.<br><br>
The simulation is fully CPU-based and does not rely on external physics libraries such as Bullet, PhysX, or Havok.

## Key Features

- **Tetrahedral Soft Body Physics** - Simulates deformable objects using position-based dynamics (PBD/XPBD) on a tetrahedral mesh
- **Mesh Skinning** - High-resolution visual mesh is skinned to the simulated tetrahedral mesh using barycentric coordinates
- **Interactive Grabbing** - Vertices can be grabbed and moved with the mouse in real time
- **Wicked Engine Integration** - Full rendering, GUI, and camera controls via Wicked Engine

## Performance Notes

⚠️ **This implementation is suitable for educational purposes and moderate mesh sizes.**
- All simulation and skinning is performed on the CPU.
- No GPU acceleration or parallelization is used.
- Performance may degrade with very high-resolution meshes or many soft bodies.

**Recommended use:** Educational purposes and small to medium-scale soft body simulations.

## How It Works

The simulation loop follows this pattern:

1. **Integrate** - Update positions and velocities of tetrahedral mesh vertices with gravity
2. **Constraint Projection** - Enforce edge and volume constraints using PBD/XPBD
3. **Skinning** - Update the high-resolution visual mesh by interpolating tetrahedral vertex positions using barycentric coordinates
4. **Visualization** - Update mesh data and upload to GPU for rendering

Each frame is processed in real time. The skinning info is computed once at initialization.

## Controls

- **WASD / Arrow Keys** - Move camera
- **Right Mouse / Middle Mouse** - Rotate camera
- **Left Mouse** - Pick and drag a vertex
- **GUI Buttons** - Run/Stop simulation, Restart, Squash, Add Soft Body, Show/Hide Tetrahedral Mesh, Adjust Compliance

