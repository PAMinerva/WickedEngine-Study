# Spatial Hashing Physics Sample

Implementation of lesson [11 - Hashing](https://github.com/matthias-research/pages/blob/master/tenMinutePhysics/11-hashing.html) from the [Ten Minute Physics](https://matthias-research.github.io/pages/tenMinutePhysics/) series by [Matthias Müller](https://github.com/matthias-research).

**Created by:** [PAMinerva](https://github.com/PAMinerva)<br><br>
**Credits to:** [Matthias Müller](https://github.com/matthias-research) for the original concept and implementation.<br><br>
**Powered by:** [Wicked Engine](https://github.com/turanszkij/WickedEngine)<br>
Special thanks to [Turánszki János](https://github.com/turanszkij) for creating Wicked Engine and making it available under the MIT license.

## Overview

This project demonstrates how to implement a simple spatial hash grid for efficient broad-phase collision detection between many moving spheres ("balls") using the Wicked Engine C++ API.<br>
The simulation is fully CPU-based and does not rely on external physics libraries such as Bullet, PhysX, or Havok.

## Key Features

- **Spatial Hash Grid** - Efficient neighbor search for collision detection
- **Broad-Phase Collisions** - Fast detection of potential overlapping pairs
- **Naive Narrow-Phase** - Simple sphere-sphere collision response

## Performance Notes

⚠️ **This implementation is O(n) for integration and O(k) for neighbor search, where k is the average number of neighbors per ball.**

Performance is suitable for hundreds of balls on modern hardware, but will degrade if the number of balls or cell density is too high. This is because:
- The hash function may produce collisions (cell aliasing)
- No parallel processing or GPU acceleration
- No advanced broad-phase (e.g. sweep and prune, BVH)

**Recommended use:** Educational purposes and small to medium-scale simulations.

## How It Works

The simulation loop follows this pattern:

1. **Integrate** - Update positions and velocities
2. **Spatial Hashing** - Assign balls to grid cells using a hash function
3. **Collision Detection** - Query neighboring cells for potential overlaps
4. **Collision Response** - Resolve sphere-sphere and world boundary collisions
5. **Visualization** - Update sphere positions and colors in the scene

Each frame is processed in real time. The hash grid is rebuilt every frame.

## Controls

- **WASD** - Move camera
- **Right Mouse** - Rotate camera
- **GUI Buttons** - Toggle collisions, reset, stop, run.
