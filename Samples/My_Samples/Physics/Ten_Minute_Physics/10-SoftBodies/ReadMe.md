# Soft Body Physics Sample

Implementation of lesson [10 - Soft Bodies](https://github.com/matthias-research/pages/blob/master/tenMinutePhysics/10-softBodies.html) from the [Ten Minute Physics](https://matthias-research.github.io/pages/tenMinutePhysics/) series by [Matthias Müller](https://github.com/matthias-research).

**Created by:** [PAMinerva](https://github.com/PAMinerva)<br>
**Credits to:** [Matthias Müller](https://github.com/matthias-research) for the original concept and implementation.

## Overview

This project demonstrates how to create simple tetrahedral meshes and simulate soft body physics using the Wicked Engine C++ API.<br>
It uses Position-Based Dynamics (PBD), and the physics is computed entirely within the application without relying on external physics libraries such as Bullet, PhysX, or Havok.

## Key Features

- **No Physics Engine** - Pure constraint-based simulation
- **Tetrahedral Mesh** - Volumetric representation using tetrahedrons
- **PBD Approach** - Position-based dynamics with edge and volume constraints
- **Interactive** - Click and drag to grab and manipulate soft bodies

## Performance Warning

⚠️ **This implementation is O(n) per soft body and does NOT scale well!**

Performance will degrade significantly after **3-5 soft bodies** depending on the hardware used and mesh complexity. This is because:
- Each soft body is solved independently without spatial optimization
- No broad-phase collision detection
- Constraint solving is done serially
- No parallel processing or GPU acceleration

**Recommended use:** Educational purposes and small-scale simulations only.

## How It Works

The simulation loop follows the PBD pattern:

1. **PreSolve** - Apply gravity and explicit integration
2. **Solve** - Enforce edge length and volume constraints
3. **PostSolve** - Update velocities based on position changes

Each frame is subdivided into multiple substeps for stability.

## Controls

- **WASD** - Move camera
- **Right Mouse** - Rotate camera
- **Left Mouse** - Grab and drag soft bodies
- **GUI Buttons** - Run/Stop, Restart, Squash, Add Body
