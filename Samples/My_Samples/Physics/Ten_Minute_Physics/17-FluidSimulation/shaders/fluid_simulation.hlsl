// ============================================================================
// fluid_simulation.hlsl
//
// All simulation compute kernels for the project's 3D Eulerian fluid solver.
//
// This file does not implement one single monolithic shader. Instead, it
// contains the individual compute passes that are dispatched in sequence from
// the C++ simulation update code. Together they implement the project's full
// fluid step on a staggered MAC grid:
//
// 1. classify solid cells and obstacle cells
// 2. inject/source smoke and temperature, and apply external forces
// 3. enforce boundary conditions on walls and moving obstacles
// 4. compute velocity divergence
// 5. solve pressure with Red-Black Gauss-Seidel + SOR
// 6. project velocity to make the field approximately divergence-free
// 7. optionally compute curl and apply vorticity confinement
// 8. advect velocity, smoke, and temperature
// 9. optionally run MacCormack correction passes to reduce advection diffusion
// 10. re-apply boundary conditions where needed after projection/advection
//
// Each pass is compiled through permutation defines, so the same HLSL source
// provides all simulation kernels:
//   SET_OBSTACLE, APPLY_FORCES, BOUNDARY, DIVERGENCE,
//   PRESSURE_RED, PRESSURE_BLACK, PROJECT,
//   COMPUTE_CURL, APPLY_VORTICITY,
//   ADVECT_VELOCITY, ADVECT_SMOKE, ADVECT_TEMPERATURE,
//   MACCORMACK_VEL, MACCORMACK_SMOKE, MACCORMACK_TEMP
//
// In short: this file is the GPU-side implementation of the simulation
// pipeline actually used by the sample, not a direct transcript of an
// external tutorial.
// ============================================================================

#include "fluid_common.hlsli"

// ============================================================================
// KERNEL: setObstacleCS
// Mark cells as solid/fluid. For the obstacle sphere, also inject velocity
// on the faces of solid cells.
// ============================================================================
#ifdef SET_OBSTACLE
[numthreads(FLUID_THREADS_3D_X, FLUID_THREADS_3D_Y, FLUID_THREADS_3D_Z)]
void setObstacleCS(uint3 dtid : SV_DispatchThreadID)
{
	// Get cell indices from thread ID
	// Used to index into textures and compute the cell-center position in world space.
    int i = (int)dtid.x;
    int j = (int)dtid.y;
    int k = (int)dtid.z;

	// Bounds check: skip threads outside the grid dimensions.
    if (i >= cb.gridX || j >= cb.gridY || k >= cb.gridZ)
        return;

	// Compute the cell-center position in world space.
    float h = cb.cellSize;
    float cx = (float)i * h + h * 0.5;
    float cy = (float)j * h + h * 0.5;
    float cz = (float)k * h + h * 0.5;

    // Convert the obstacle center from world space to domain-local world space.
    float obsLocalX = cb.obsCenterX - (-cb.gridX * h * 0.5);
    float obsLocalY = cb.obsCenterY;
    float obsLocalZ = cb.obsCenterZ - (-cb.gridZ * h * 0.5);

	// Compute squared distance from the cell center to the obstacle center in domain-local world space.
	// Also compute squared obstacle radius for comparison (avoid sqrt for efficiency).
    float dx = cx - obsLocalX;
    float dy = cy - obsLocalY;
    float dz = cz - obsLocalZ;
    float dist2 = dx * dx + dy * dy + dz * dz;
    float r2 = cb.obsRadius * cb.obsRadius;

    // All walls are solid (closed box).
    bool isBoundary = (i == 0 || i == cb.gridX - 1 ||
                       j == 0 || j == cb.gridY - 1 ||
                       k == 0 || k == cb.gridZ - 1);

	// A cell is an interior obstacle if its center is within the sphere radius.
    bool isObs = (dist2 < r2);

	// If this cell is solid (either boundary or interior obstacle)...
    if (isObs || isBoundary)
    {
		// Mark as solid in the solid texture (0 = solid, 1 = fluid).
        texSolid[int3(i, j, k)] = 0;

		// If the current cell is an interior obstacle...
        if (isObs)
        {
			// Set velocity on ALL 6 faces of this cell to the obstacle velocity.
            texU[int3(i, j, k)]     = cb.obsVelX;
            texU[int3(i + 1, j, k)] = cb.obsVelX;
            texV[int3(i, j, k)]     = cb.obsVelY;
            texV[int3(i, j + 1, k)] = cb.obsVelY;
            texW[int3(i, j, k)]     = cb.obsVelZ;
            texW[int3(i, j, k + 1)] = cb.obsVelZ;

			// Clear pressure and smoke density in obstacle cells to prevent unwanted forces and visual artifacts.
            texPressure[int3(i, j, k)] = 0.0;
            texSmoke[int3(i, j, k)] = 0.0;
        }
    }
    else // otherwise, this cell is fluid.
    {
        texSolid[int3(i, j, k)] = 1;
    }
}
#endif

// ============================================================================
// KERNEL: applyForcesCS
// Apply gravity, temperature-based buoyancy, density weight,
// and Gaussian smoke+temperature injection.
// ============================================================================
#ifdef APPLY_FORCES
[numthreads(FLUID_THREADS_3D_X, FLUID_THREADS_3D_Y, FLUID_THREADS_3D_Z)]
void applyForcesCS(uint3 dtid : SV_DispatchThreadID)
{
	// Get cell indices from thread ID
    int i = (int)dtid.x;
    int j = (int)dtid.y;
    int k = (int)dtid.z;

	// Bounds check: skip threads outside the grid dimensions.
    if (i >= cb.gridX || j >= cb.gridY || k >= cb.gridZ)
        return;

	// If this cell is solid, skip it. Forces only apply to fluid cells.
    if (texSolid[int3(i, j, k)] == 0)
        return;

	// Cell size (h) is needed for computing positions.
    float h = cb.cellSize;

    // --- Gravity ---
    // Only apply to faces where BOTH cells sharing the face are fluid.
    // Faces touching solid cells are managed by boundaryCS.
	//
	// If this cell is not the top row, and the neighbor above is fluid,
	// apply gravity to the top face (v[j+1]) of the current cell (i, j, k),
	// since that's the face between this cell and the one above it (i, j+1, k),
	// which is where gravity acts in the y direction between these two cells.
    if (abs(cb.gravY) > 0.0 && j < cb.gridY - 1 && texSolid[int3(i, j + 1, k)] != 0)
    {
        texV[int3(i, j + 1, k)] += cb.gravY * cb.dt;
    }
	// Similarly for x and z directions: only apply if neighbor in that direction is fluid.
    if (abs(cb.gravX) > 0.0 && i < cb.gridX - 1 && texSolid[int3(i + 1, j, k)] != 0)
    {
        texU[int3(i + 1, j, k)] += cb.gravX * cb.dt;
    }
    if (abs(cb.gravZ) > 0.0 && k < cb.gridZ - 1 && texSolid[int3(i, j, k + 1)] != 0)
    {
        texW[int3(i, j, k + 1)] += cb.gravZ * cb.dt;
    }

    // --- Temperature-based buoyancy + density weight ---
    // F_y = (T - T_ambient) * buoyancy - |smoke| * weight
    // Only on faces between two fluid cells (neighbor above must be fluid).
    if (j < cb.gridY - 1 && texSolid[int3(i, j + 1, k)] != 0)
    {
		// Temperature is stored at cell centers, so average the temperature of
		// the current cell and the one above it to get the temperature at the face between them.
        float tempHere  = texTemp[int3(i, j, k)];
        float tempAbove = texTemp[int3(i, j + 1, k)];
        float avgTemp   = 0.5 * (tempHere + tempAbove);

		// Similarly, average the smoke density of the current cell and the one above it to get
		// the smoke density at the face.
        float smokeHere  = texSmoke[int3(i, j, k)];
        float smokeAbove = texSmoke[int3(i, j + 1, k)];
        float avgSmoke   = 0.5 * (smokeHere + smokeAbove);

		// Buoyancy is the force exerted by a fluid opposing the weight of a partially or fully
		// immersed object (which may also be a parcel of fluid).
		// For smoke, buoyancy is typically modeled as proportional to the temperature difference from ambient,
		// since hotter fluid is less dense and rises. However, in this simulation we assume constant density for simplicity,
		// so we will approximate buoyancy as a direct force proportional to the temperature difference and smoke concentration,
		// without explicitly modeling density variations.
		// The first term (avgTemp - ambientTemp) * buoyancy causes hotter fluid to rise,
		// while the second term -avgSmoke * weight causes denser smoke to sink.
        float force = (avgTemp - cb.ambientTemp) * cb.buoyancy - avgSmoke * cb.fluidWeight;

		// Apply the buoyancy force to the vertical velocity (v) at the face between this cell and the one above it.
		// If force > 0 (positive buoyancy), it will increase the upward velocity, causing the fluid to rise:
		// Remember that cell faces store velocity components (scalar values), so if a v component (stored at horizontal faces)
		// is positive, it will point in the positive y direction (upwards).
        texV[int3(i, j + 1, k)] += force * cb.dt;
    }

    // --- Gaussian smoke + temperature injection ---
	// if smoke effect is enabled...
    if (cb.enableSmoke)
    {
		// Compute the cell-center position in world space.
        float cx = (float)i * h + h * 0.5;
        float cy = (float)j * h + h * 0.5;
        float cz = (float)k * h + h * 0.5;

		// Convert the smoke source center from world space to domain-local world space.
        float srcX = cb.smokeSrcX - (-cb.gridX * h * 0.5);
        float srcY = cb.smokeSrcY;
        float srcZ = cb.smokeSrcZ - (-cb.gridZ * h * 0.5);

		// Compute distance from the cell center to the smoke source center in domain-local world space.
        float dx = cx - srcX;
        float dy = cy - srcY;
        float dz = cz - srcZ;
        float dist = sqrt(dx * dx + dy * dy + dz * dz);

		// If the cell is within the smoke source (the emitter sphere) radius...
        float r = cb.smokeSrcRadius;
        if (dist < r)
        {
            // Compute normalized distance (0 at center, 1 at radius edge)
            float normDist = dist / r;

            // Gaussian falloff: smooth spatial injection profile
            // exp(-4.0 * normDist^2) gives ~0.018 (≈1.8%) at normDist=1 (sphere edge),
            // ensuring smooth transition from full injection at center to ~zero at boundary.
			// The term injection is used (instead of emission) to reflect that we are directly
			// setting the density and temperature values in the grid based on the source properties,
			// rather than modeling a continuous emission process over time.
            float gaussian = exp(-normDist * normDist * 4.0);

            // Inject smoke density: take max to preserve existing smoke while adding new
            texSmoke[int3(i, j, k)] = max(texSmoke[int3(i, j, k)], cb.smokeInflowDensity * gaussian);

            // Inject temperature: similarly, take max to avoid cooling hot smoke with cold inflow
            texTemp[int3(i, j, k)]  = max(texTemp[int3(i, j, k)],  cb.tempInjection * gaussian);
        }
    }
}
#endif

// ============================================================================
// KERNEL: boundaryCS
// Enforce boundary conditions on ALL solid cells (walls + obstacles).
//
// This kernel performs two tasks:
//
//   1) NO-PENETRATION (normal velocity):
//      Set the normal velocity on all 6 faces of every solid cell.
//        - Domain walls       → 0        (static walls)
//        - Interior obstacles → obsVel   (moving sphere)
//
//      In a MAC grid, velocity components live on cell faces. A face shared
//      between a solid cell and a fluid cell stores the normal velocity at
//      the boundary. By writing sVelX/Y/Z to ALL 6 faces of every solid
//      cell, we directly impose:  v_normal(shared face) = solid velocity.
//
//      projectCS (the pressure gradient subtraction) will NOT overwrite
//      these values, because it only updates a face when BOTH cells on
//      either side are fluid:
//        if (texSolid[neighbor] != 0)   // != 0 means fluid
//            texU[face] -= scale * (p_this - p_neighbor);
//      If one side is solid, projectCS skips that face entirely. Therefore
//      the shared face keeps the value we set here.
//
//      The fluid adapts via the pressure solver: divergenceCS reads the
//      fixed shared face, any imbalance creates non-zero divergence in
//      the neighboring fluid cell, and the pressure solver builds up
//      pressure to cancel it. projectCS then redistributes the correction
//      onto the fluid cell's OTHER (fluid-fluid) faces, redirecting flow
//      sideways — never through the solid.
//
//   2) FREE-SLIP TANGENTIAL EXTRAPOLATION (domain walls only):
//      Copy tangential velocity from the nearest interior fluid cell.
//      This ensures that fluid slides freely along walls without friction.
//      This step overwrites the tangential components set in (1), but NOT
//      the normal ones.
//
// WHY THIS KERNEL RUNS 3 TIMES PER FRAME (see fluid_sim.cpp):
//
//   Step 5  — Pre-projection:
//     Steps 1–4 (external forces, buoyancy, vorticity confinement) write
//     directly into the velocity textures and can corrupt shared-face
//     values. We restore them here so that divergenceCS reads correct
//     boundary velocities and the pressure solver computes accurate
//     pressures for fluid cells adjacent to solids.
//
//   Step 9  — Post-projection:
//     projectCS preserves the normal velocities on shared faces (see above),
//     but it DOES update velocities in fluid cells adjacent to walls. The
//     free-slip tangential extrapolation (task 2) copies tangential velocity
//     from those fluid cells, so it must re-run after projection to pick up
//     the newly corrected fluid velocities. Without this pass, tangential
//     wall values would be stale (from pre-projection).
//
//   Step 15 — Post-advection:
//     Semi-Lagrangian advection interpolates velocities from arbitrary
//     positions and can write incorrect values onto shared faces. We
//     re-enforce all boundary conditions after advection to clean up.
//
// ============================================================================
#ifdef BOUNDARY
[numthreads(FLUID_THREADS_3D_X, FLUID_THREADS_3D_Y, FLUID_THREADS_3D_Z)]
void boundaryCS(uint3 dtid : SV_DispatchThreadID)
{
	// Get cell indices from thread ID
    int i = (int)dtid.x;
    int j = (int)dtid.y;
    int k = (int)dtid.z;

	// Bounds check: skip threads outside the grid dimensions.
    if (i >= cb.gridX || j >= cb.gridY || k >= cb.gridZ)
        return;

    // Skip fluid cells: only enforce at solid boundaries
    if (texSolid[int3(i, j, k)] != 0)
        return;

    // Classify: domain wall vs interior obstacle (sphere)
    bool isDomainWall = (i == 0 || i == cb.gridX - 1 ||
                         j == 0 || j == cb.gridY - 1 ||
                         k == 0 || k == cb.gridZ - 1);

	// For domain walls, normal velocity = 0 (static wall).
	// For interior obstacles, normal velocity = obsVel (moving sphere velocity; see SetObstacle call in SampleRenderPath::Update in sample.cpp).
    float sVelX = isDomainWall ? 0.0 : cb.obsVelX;
    float sVelY = isDomainWall ? 0.0 : cb.obsVelY;
    float sVelZ = isDomainWall ? 0.0 : cb.obsVelZ;

    // --- No-penetration: set normal velocity on all 6 faces ---
    texU[int3(i, j, k)]     = sVelX;   // left face
    texU[int3(i + 1, j, k)] = sVelX;   // right face
    texV[int3(i, j, k)]     = sVelY;   // bottom face
    texV[int3(i, j + 1, k)] = sVelY;   // top face
    texW[int3(i, j, k)]     = sVelZ;   // back face
    texW[int3(i, j, k + 1)] = sVelZ;   // front face

    // --- Free-slip tangential extrapolation (domain walls only) ---
    // Copy tangential velocity from the nearest interior fluid cell.
    // This runs AFTER the 6-face assignment above, overwriting the
    // tangential components (but NOT the normal ones).
    if (!isDomainWall)
        return;

    // Left wall (i == 0)
	// if this cell is at the left boundary and the neighbor to the right is fluid,
	// copy tangential velocities (V and W) from that neighbor.
	// Remember that (i, j, k) is the current solid cell, so the neighbor to the right is (i+1, j, k).
	// Moreover, texV[i, j, k] is the vertical velocity at the bottom face of the current cell,
	// texV[i, j+1, k] is the vertical velocity at the top face of the current cell,
	// texW[i, j, k] is the depth velocity at the back face of the current cell,
	// texW[i, j, k+1] is the depth velocity at the front face of the current cell.
    if (i == 0 && i + 1 < cb.gridX && texSolid[int3(i + 1, j, k)] != 0)
    {
        texV[int3(i, j, k)]     = texV[int3(i + 1, j, k)];
        texV[int3(i, j + 1, k)] = texV[int3(i + 1, j + 1, k)];
        texW[int3(i, j, k)]     = texW[int3(i + 1, j, k)];
        texW[int3(i, j, k + 1)] = texW[int3(i + 1, j, k + 1)];
    }
    // Right wall (i == gridX-1): neighbour (i-1, j, k) is to the left, tangential = V, W
    if (i == cb.gridX - 1 && i - 1 >= 0 && texSolid[int3(i - 1, j, k)] != 0)
    {
        texV[int3(i, j, k)]     = texV[int3(i - 1, j, k)];
        texV[int3(i, j + 1, k)] = texV[int3(i - 1, j + 1, k)];
        texW[int3(i, j, k)]     = texW[int3(i - 1, j, k)];
        texW[int3(i, j, k + 1)] = texW[int3(i - 1, j, k + 1)];
    }
// Bottom wall (j == 0): neighbor (i, j+1, k) is above, tangential = U, W
	// tangential = U, W
    if (j == 0 && j + 1 < cb.gridY && texSolid[int3(i, j + 1, k)] != 0)
    {
        texU[int3(i, j, k)]     = texU[int3(i, j + 1, k)];
        texU[int3(i + 1, j, k)] = texU[int3(i + 1, j + 1, k)];
        texW[int3(i, j, k)]     = texW[int3(i, j + 1, k)];
        texW[int3(i, j, k + 1)] = texW[int3(i, j + 1, k + 1)];
    }
    // Top wall (j == gridY-1): neighbor (i, j-1, k) is below, tangential = U, W
    if (j == cb.gridY - 1 && j - 1 >= 0 && texSolid[int3(i, j - 1, k)] != 0)
    {
        texU[int3(i, j, k)]     = texU[int3(i, j - 1, k)];
        texU[int3(i + 1, j, k)] = texU[int3(i + 1, j - 1, k)];
        texW[int3(i, j, k)]     = texW[int3(i, j - 1, k)];
        texW[int3(i, j, k + 1)] = texW[int3(i, j - 1, k + 1)];
    }
    // Back wall (k == 0): neighbor (i, j, k+1) is in front, tangential = U, V
    if (k == 0 && k + 1 < cb.gridZ && texSolid[int3(i, j, k + 1)] != 0)
    {
        texU[int3(i, j, k)]     = texU[int3(i, j, k + 1)];
        texU[int3(i + 1, j, k)] = texU[int3(i + 1, j, k + 1)];
        texV[int3(i, j, k)]     = texV[int3(i, j, k + 1)];
        texV[int3(i, j + 1, k)] = texV[int3(i, j + 1, k + 1)];
    }
    // Front wall (k == gridZ-1): neighbor (i, j, k-1) is in back, tangential = U, V
    if (k == cb.gridZ - 1 && k - 1 >= 0 && texSolid[int3(i, j, k - 1)] != 0)
    {
        texU[int3(i, j, k)]     = texU[int3(i, j, k - 1)];
        texU[int3(i + 1, j, k)] = texU[int3(i + 1, j, k - 1)];
        texV[int3(i, j, k)]     = texV[int3(i, j, k - 1)];
        texV[int3(i, j + 1, k)] = texV[int3(i, j + 1, k - 1)];
    }
}
#endif

// ============================================================================
// KERNEL: divergenceCS
// Compute the divergence of the velocity field for each fluid cell.
//
// ── Definition ─────────────────────────────────────────────────────────
// The divergence of a 3D velocity field v = (u, v, w) is a scalar:
//
//   ∇·v = ∂u/∂x + ∂v/∂y + ∂w/∂z
//
// It measures the net rate at which fluid flows out of (positive) or
// into (negative) an infinitesimal volume element.  For an
// incompressible fluid, ∇·v = 0 everywhere — what flows in must
// flow out.  Any non-zero divergence signals local compression or
// expansion, which the pressure projection step will correct.
//
// ── Forward differences on a staggered (MAC) grid ─────────────────────
// On a collocated grid the derivative ∂u/∂x would typically use
// *central* differences:  (u(i+1) − u(i−1)) / (2h).
//
// On a MAC grid, velocity components already live at cell *faces*:
//   u is stored at the left/right faces of each cell,
//   v at the bottom/top faces, w at the back/front faces.
//
// For cell (i,j,k) the two u-faces are exactly one cell apart:
//   u_right = texU[i+1, j, k]   (right face)
//   u_left  = texU[i,   j, k]   (left  face)
//
// Their midpoint falls precisely at the cell center, so
//
//   ∂u/∂x ≈ (u_right − u_left) / h
//
// is already a *centered* second-order approximation without needing
// the (i+1)/(i−1) stencil — the staggering provides the centering
// for free.  The same applies to ∂v/∂y and ∂w/∂z.
//
// ── Result ─────────────────────────────────────────────────────────────
// The full discrete divergence for cell (i,j,k) is:
//
//   d = (u_{i+1} - u_i + v_{j+1} - v_j + w_{k+1} - w_k) / h
//
// If d ≠ 0, the cell violates incompressibility. The pressure solver
// (pressureRedCS / pressureBlackCS) will compute pressures to cancel
// this imbalance, and projectCS will subtract the pressure gradient to
// produce a divergence-free velocity field.
//
// This kernel also clears the pressure texture to zero, so the
// iterative pressure solver starts from a clean slate each frame.
// ============================================================================
#ifdef DIVERGENCE
[numthreads(FLUID_THREADS_3D_X, FLUID_THREADS_3D_Y, FLUID_THREADS_3D_Z)]
void divergenceCS(uint3 dtid : SV_DispatchThreadID)
{
	// Get cell indices from thread ID
    int i = (int)dtid.x;
    int j = (int)dtid.y;
    int k = (int)dtid.z;

	// Bounds check: skip threads outside the grid dimensions.
    if (i >= cb.gridX || j >= cb.gridY || k >= cb.gridZ)
        return;

	// Skip solid cells: divergence is only defined for fluid cells.
    if (texSolid[int3(i, j, k)] == 0)
    {
        texDivergence[int3(i, j, k)] = 0.0;
        return;
    }

    float h = cb.cellSize;

    // Sample velocities on the 6 faces of cell (i,j,k)
    float uRight = texU[int3(i + 1, j, k)];
    float uLeft  = texU[int3(i, j, k)];
    float vTop   = texV[int3(i, j + 1, k)];
    float vBot   = texV[int3(i, j, k)];
    float wFront = texW[int3(i, j, k + 1)];
    float wBack  = texW[int3(i, j, k)];

    // Divergence: d = (u_right - u_left + v_top - v_bot + w_front - w_back) / h
    float div = (uRight - uLeft + vTop - vBot + wFront - wBack) / h;
    texDivergence[int3(i, j, k)] = div;

    // Clear pressure before solve (first iteration starts fresh)
    texPressure[int3(i, j, k)] = 0.0;
}
#endif

// ============================================================================
// KERNEL: pressureRedCS / pressureBlackCS
//
// ┌─────────────────────────────────────────────────────────────────────┐
// │  CONTEXT: WHERE ARE WE IN THE SIMULATION PIPELINE?                 │
// └─────────────────────────────────────────────────────────────────────┘
//
// At this point in the frame we have an intermediate velocity field v*
// (stored in texU/texV/texW) that was produced by applying external
// forces (gravity, buoyancy, etc.) to the previous frame's velocity.
//
// Problem: v* is NOT divergence-free.  Some cells have more fluid
// flowing in than out (compression) or vice versa (expansion).  For an
// incompressible fluid this is physically wrong — fluid cannot be
// created or destroyed.  We need to find a pressure field p that, when
// we subtract its gradient from v*, gives a new velocity v^{n+1} that
// IS divergence-free.  This kernel computes that pressure.
//
// ┌─────────────────────────────────────────────────────────────────────┐
// │  STEP 1: DERIVING THE PRESSURE POISSON EQUATION (PPE)              │
// └─────────────────────────────────────────────────────────────────────┘
//
// The projection step (see projectCS) corrects v* into a divergence-free
// field v^{n+1} by subtracting the pressure gradient:
//
//   v^{n+1} = v*  −  (Δt / ρ) · ∇p                                 (A)
//
// where Δt = time step, ρ = fluid density, ∇p = pressure gradient.
//
// We want v^{n+1} to be divergence-free, i.e.:
//
//   ∇ · v^{n+1} = 0                                                 (B)
//
// Take the divergence of both sides of (A):
//
//   ∇ · v^{n+1} = ∇ · v*  −  (Δt / ρ) · ∇ · (∇p)
//               = ∇ · v*  −  (Δt / ρ) · ∇²p
//
// Substitute constraint (B) (left side = 0):
//
//   0 = ∇ · v*  −  (Δt / ρ) · ∇²p
//
// Rearrange to isolate the Laplacian of pressure:
//
//   ∇²p  =  (ρ / Δt) · (∇ · v*)                                    (1)
//
// This is the Pressure Poisson Equation (PPE).  It says: the Laplacian
// of the pressure field equals the divergence of the intermediate
// velocity, scaled by ρ/Δt.  If we solve (1) for p and then apply
// equation (A) in projectCS, the resulting velocity is guaranteed to be
// divergence-free.
//
// ┌─────────────────────────────────────────────────────────────────────┐
// │  STEP 2: DISCRETIZING THE LAPLACIAN (LEFT SIDE OF THE PPE)         │
// └─────────────────────────────────────────────────────────────────────┘
//
// We need to approximate ∇²p on our discrete grid.  Start from the 1D
// definition of the second derivative:
//
//   d²p/dx²  ≈  (p(x+h) − 2·p(x) + p(x−h)) / h²
//
// This is the standard central finite-difference formula: we evaluate
// the function at three equally-spaced points and divide by h².
//
// In 3D the Laplacian is the sum of the second derivatives along each
// axis:
//
//   ∇²p = ∂²p/∂x² + ∂²p/∂y² + ∂²p/∂z²
//
// Applying the 1D formula to each axis for cell (i,j,k):
//
//   ∂²p/∂x² ≈ (p_{i+1,j,k} + p_{i−1,j,k} − 2·p_{i,j,k}) / h²
//   ∂²p/∂y² ≈ (p_{i,j+1,k} + p_{i,j−1,k} − 2·p_{i,j,k}) / h²
//   ∂²p/∂z² ≈ (p_{i,j,k+1} + p_{i,j,k−1} − 2·p_{i,j,k}) / h²
//
// Summing all three:
//
//   ∇²p ≈ (p_{i+1} + p_{i−1} + p_{j+1} + p_{j−1} + p_{k+1} + p_{k−1}
//           − 6·p_{i,j,k}) / h²
//
//       = (Σ_neighbors p  −  6 · p_{i,j,k}) / h²
//
// This is the 7-point Laplacian stencil: 6 neighbors plus the center.
//
// ── Handling solid cells (boundary conditions) ──────────────────────
//
// When a neighbor is solid (wall, obstacle), we use a Neumann boundary
// condition: ∂p/∂n = 0 at the solid face (no pressure gradient into
// the wall).  In discrete terms this means p_solid = p_{i,j,k}, so
// the contribution of that neighbor to the stencil becomes:
//
//   p_solid − p_{i,j,k} = 0
//
// That term simply vanishes.  Instead of always having 6 neighbors, we
// count only the fluid (non-solid) ones in s_tot.  The stencil becomes:
//
//   ∇²p ≈ (Σ_{fluid neighbors} p  −  s_tot · p_{i,j,k}) / h²      (2)
//
// ┌─────────────────────────────────────────────────────────────────────┐
// │  STEP 3: DISCRETIZING THE DIVERGENCE (RIGHT SIDE OF THE PPE)       │
// └─────────────────────────────────────────────────────────────────────┘
//
// The right-hand side of (1) is (ρ/Δt) · ∇·v*.  The divergenceCS kernel
// already computed the discrete divergence and stored it in texDivergence:
//
//   div = (uR − uL + vT − vB + wF − wB) / h
//
// where uR, uL, etc. are the face velocities on the 6 faces of the cell.
// This is the standard finite-difference approximation of ∇·v* on a MAC
// (staggered) grid — see the comments in divergenceCS for the full
// derivation.  So:
//
//   ∇·v*  ≈  div     (the value stored in texDivergence)
//
// ┌─────────────────────────────────────────────────────────────────────┐
// │  STEP 4: COMBINING BOTH SIDES → THE UPDATE FORMULA                 │
// └─────────────────────────────────────────────────────────────────────┘
//
// Substitute (2) and div into the PPE (1):
//
//   (Σ p_neighbors − s_tot · p_ijk) / h²  =  (ρ / Δt) · div
//
// Multiply both sides by h²:
//
//   Σ p_neighbors − s_tot · p_ijk  =  ρ · h² · div / Δt
//
// Isolate p_ijk:
//
//   p_ijk  =  (Σ p_neighbors  −  ρ · h² · div / Δt)  /  s_tot     (3)
//
// This is the formula implemented below:
//   p_new = (p_sum - density * h * h * div / dt) / s_tot
//
// ┌─────────────────────────────────────────────────────────────────────┐
// │  STEP 5: ITERATIVE SOLVING WITH GAUSS-SEIDEL                       │
// └─────────────────────────────────────────────────────────────────────┘
//
// Equation (3) expresses p_ijk in terms of its neighbors' pressures.
// But those neighbors' pressures are unknowns too!  Every cell's
// pressure depends on every other cell's pressure.  This is a huge
// system of coupled linear equations (one per cell).
//
// We solve it iteratively with Gauss-Seidel: start with p = 0
// everywhere (done in divergenceCS), then repeatedly sweep through all
// cells applying formula (3).  Each time we apply (3) to a cell, we
// read its neighbors' *latest* values — including values that were
// already updated earlier in the same sweep.  This "use the freshest
// data" property is what distinguishes Gauss-Seidel from Jacobi
// iteration (which only uses values from the previous sweep) and makes
// it converge roughly 2× faster.
//
// ┌─────────────────────────────────────────────────────────────────────┐
// │  STEP 6: RED-BLACK COLORING (GPU-PARALLEL GAUSS-SEIDEL)            │
// └─────────────────────────────────────────────────────────────────────┘
//
// Classic Gauss-Seidel is inherently serial: updating cell A changes
// the value that cell B reads.  On a GPU we need parallelism.
//
// Red-Black coloring solves this.  Color every cell like a 3D
// checkerboard:
//
//   Red  cells: (i + j + k) % 2 == 0
//   Black cells: (i + j + k) % 2 == 1
//
// Key property: every Red cell's 6 axis-aligned neighbors are ALL
// Black, and vice versa.  (Think of a chessboard — a white square is
// surrounded only by black squares.)
//
// This means:
//   1. We can update ALL Red cells in parallel — they only read Black
//      values, and no Red cell reads another Red cell.  No race
//      conditions.
//   2. After a GPU barrier (ensuring all Red writes are visible), we
//      update ALL Black cells in parallel — they only read the
//      just-updated Red values.
//
// One Red pass + one Black pass = one complete Gauss-Seidel iteration.
// We repeat for cb.numPressureIters iterations per time step.
//
// This is why there are TWO shader permutations (PRESSURE_RED and
// PRESSURE_BLACK) dispatched alternately from the CPU side.
//
// ┌─────────────────────────────────────────────────────────────────────┐
// │  STEP 7: SUCCESSIVE OVER-RELAXATION (SOR)                          │
// └─────────────────────────────────────────────────────────────────────┘
//
// Even with Gauss-Seidel, convergence can be slow for large grids.
// SOR accelerates it by "overshooting" each update:
//
//   p_sor = (1 − ω) · p_old  +  ω · p_new
//
// where ω (omega, cb.overRelaxation) controls how much we overshoot:
//
//   ω = 1.0  →  Standard Gauss-Seidel, no over-relaxation.
//               p_sor = p_new exactly.
//
//   1 < ω < 2  →  Over-relaxation.  We push the update PAST the
//                  Gauss-Seidel value, betting that the true solution
//                  lies further in that direction.  This gamble pays off
//                  for smooth error modes (low-frequency errors) which
//                  dominate in Poisson problems.
//
//   ω = 1.9  →  Near-optimal for typical 3D Poisson problems on uniform
//               grids.  Can reduce the number of iterations needed by
//               5-10× compared to ω = 1.  This is our default.
//
//   ω ≥ 2.0  →  The iteration DIVERGES.  Must be avoided.
//
// Intuitively: plain Gauss-Seidel corrects each cell "just enough" to
// satisfy equation (3) given current neighbors.  But since neighbors
// haven't all been updated yet, this correction is conservative.  SOR
// says "I know I'm undershooting, so let me go further."  The art is
// choosing ω large enough to speed things up without overshooting so
// much that it becomes unstable.
//
// ============================================================================
#if defined(PRESSURE_RED) || defined(PRESSURE_BLACK)

#ifdef PRESSURE_RED
[numthreads(FLUID_THREADS_3D_X, FLUID_THREADS_3D_Y, FLUID_THREADS_3D_Z)]
void pressureRedCS(uint3 dtid : SV_DispatchThreadID)
#else
[numthreads(FLUID_THREADS_3D_X, FLUID_THREADS_3D_Y, FLUID_THREADS_3D_Z)]
void pressureBlackCS(uint3 dtid : SV_DispatchThreadID)
#endif
{
	// Get cell indices from thread ID
    int i = (int)dtid.x;
    int j = (int)dtid.y;
    int k = (int)dtid.z;

	// Bounds check: skip threads outside the grid dimensions.
    if (i >= cb.gridX || j >= cb.gridY || k >= cb.gridZ)
        return;

    // ── Red-Black selection ─────────────────────────────────────────
    // Determine this cell's color from the checkerboard pattern.
    // Red = even parity, Black = odd parity.  Each permutation
    // (PRESSURE_RED / PRESSURE_BLACK) processes only its own color
    // and early-exits on the other.
    int parity = (i + j + k) % 2;
#ifdef PRESSURE_RED
    if (parity != 0) return;  // this is a Black cell → skip in Red pass
#else
    if (parity != 1) return;  // this is a Red cell → skip in Black pass
#endif

    // Solid cells don't participate in the pressure solve — they are
    // obstacles or walls, not fluid.
    if (texSolid[int3(i, j, k)] == 0)
        return;

	// Cell size
    float h = cb.cellSize;

    // ── Build the Laplacian stencil (see STEP 2 above) ──────────────
    // Walk the 6 axis-aligned neighbors.  For each fluid (non-solid)
    // neighbor, add its pressure to p_sum and increment s_tot.
    // Solid neighbors are skipped (Neumann BC: ∂p/∂n = 0).
    float s_tot = 0.0;
    float p_sum = 0.0;

    // −X neighbor (left)
    if (i > 0 && texSolid[int3(i - 1, j, k)] != 0)
    {
        s_tot += 1.0;
        p_sum += texPressure[int3(i - 1, j, k)];
    }
    // +X neighbor (right)
    if (i < cb.gridX - 1 && texSolid[int3(i + 1, j, k)] != 0)
    {
        s_tot += 1.0;
        p_sum += texPressure[int3(i + 1, j, k)];
    }
    // −Y neighbor (bottom)
    if (j > 0 && texSolid[int3(i, j - 1, k)] != 0)
    {
        s_tot += 1.0;
        p_sum += texPressure[int3(i, j - 1, k)];
    }
    // +Y neighbor (top)
    if (j < cb.gridY - 1 && texSolid[int3(i, j + 1, k)] != 0)
    {
        s_tot += 1.0;
        p_sum += texPressure[int3(i, j + 1, k)];
    }
    // −Z neighbor (back)
    if (k > 0 && texSolid[int3(i, j, k - 1)] != 0)
    {
        s_tot += 1.0;
        p_sum += texPressure[int3(i, j, k - 1)];
    }
    // +Z neighbor (front)
    if (k < cb.gridZ - 1 && texSolid[int3(i, j, k + 1)] != 0)
    {
        s_tot += 1.0;
        p_sum += texPressure[int3(i, j, k + 1)];
    }

    // If all neighbors are solid, there's nothing to solve for.
    if (s_tot < 1e-6)
        return;

    float div   = texDivergence[int3(i, j, k)];
    float p_old = texPressure[int3(i, j, k)];

    // ── Gauss-Seidel update — equation (3) from STEP 4 ─────────────
    //   p_new = (Σ p_neighbors  −  ρ · h² · div / Δt)  /  s_tot
    //
    // p_sum holds Σ p_neighbors, div is read from texDivergence, and
    // the ρ·h²/Δt factor converts the divergence into the pressure
    // units expected by the Laplacian stencil (see STEP 4 derivation).
    float p_new = (p_sum - cb.density * h * h * div / cb.dt) / s_tot;

    // ── SOR blending (STEP 7) ───────────────────────────────────────
    // Blend between the old pressure and the Gauss-Seidel update,
    // overshooting by factor ω to accelerate convergence.
    //   ω = 1.0 → no over-relaxation (pure Gauss-Seidel)
    //   ω = 1.9 → near-optimal for this type of problem (our default)
    float p_sor = (1.0 - cb.overRelaxation) * p_old + cb.overRelaxation * p_new;

    texPressure[int3(i, j, k)] = p_sor;
}
#endif

// ============================================================================
// KERNEL: projectCS
//
// ── What this kernel does ───────────────────────────────────────────────
// Now that the pressure field p has been computed by the pressure solver
// (pressureRedCS / pressureBlackCS), we apply equation (A) to correct
// the intermediate velocity v* into a divergence-free velocity v^{n+1}:
//
//   v^{n+1} = v*  −  (Δt / ρ) · ∇p                                 (A)
//
// ┌─────────────────────────────────────────────────────────────────────┐
// │  WHERE DOES EQUATION (A) COME FROM?                                 │
// └─────────────────────────────────────────────────────────────────────┘
//
// ── Background: the Navier-Stokes momentum equation ─────────────────────
// The motion of an incompressible fluid is governed by:
//
//   ∂v/∂t  =  −(v · ∇)v  −  (1/ρ) · ∇p  +  f                      (NS)
//
//   subject to:  ∇ · v = 0   (incompressibility)
//
// Note: the full Navier-Stokes equation also includes a viscosity term
// ν∇²v (diffusion of velocity).  Our simulator ignores viscosity
// entirely — the numerical diffusion introduced by semi-Lagrangian
// advection already provides enough damping in practice, and skipping
// the viscosity solve saves significant computation.
//
// Each term on the right of (NS) is an acceleration (m/s²) — a rate of
// change of velocity per unit time:
//
//   −(v · ∇)v     Advection acceleration (explained below).
//
//   −(1/ρ) · ∇p   Pressure acceleration (explained below).
//
//   f              External forces per unit mass (gravity, buoyancy,
//                  vorticity confinement, etc.).
//
// ── Why is −(v · ∇)v an acceleration? ───────────────────────────────────
// In Navier-Stokes, ∂v/∂t measures how velocity changes at a FIXED
// point in space (our grid cell).  At that fixed point, velocity can
// change for two reasons:
//
//   1. Forces push the fluid (pressure, gravity).
//   2. The flow carries fluid with a DIFFERENT velocity into that point.
//
// The −(v · ∇)v term captures reason 2.  To see how, consider a 1D
// example.  Suppose the flow moves to the right (u > 0) and velocity
// increases to the right (∂u/∂x > 0):
//
//       slower fluid          our point          faster fluid
//       ─────────────────────────●─────────────────────────→ x
//            u = 1 m/s                            u = 3 m/s
//                             flow direction →
//
// Question: at our fixed point, will velocity increase or decrease
// over time?
//
// The flow moves to the right, which means the fluid arriving at our
// point comes from the LEFT (upstream).  The fluid on the left is
// SLOWER (1 m/s).  So as this slower fluid flows in, the velocity at
// our fixed point DECREASES.  This is purely kinematic — no force
// is involved, just the flow replacing the fluid at our measurement
// point with different fluid from upstream.
//
// Now let's check the math.  The term is:
//
//   −(v · ∇)v   →   in 1D:  −u · (∂u/∂x)
//
// With u > 0 and ∂u/∂x > 0:
//
//   −u · (∂u/∂x)  <  0      →  ∂v/∂t < 0  →  velocity decreases ✓
//
// This matches our physical reasoning.  The formula works because:
//
//   (v · ∇)v  is the directional derivative of v in the direction of
//   the flow.  It measures how velocity changes as you move DOWNSTREAM
//   (in the direction the fluid is going).
//
//   The MINUS sign flips this: what matters for ∂v/∂t at a fixed point
//   is not what's downstream, but what's arriving from UPSTREAM.  If
//   velocity increases downstream (positive directional derivative),
//   then what's coming from upstream is slower, so velocity at our
//   point decreases — hence the minus.
//
// A second example to confirm:  u > 0, ∂u/∂x < 0 (velocity DECREASES
// to the right → faster fluid is upstream):
//
//   −u · (∂u/∂x)  >  0      →  ∂v/∂t > 0  →  velocity increases ✓
//
// Correct: the faster upstream fluid flows in and raises the velocity
// at our fixed point.
//
// Units check: u [m/s] × ∂u/∂x [s⁻¹] = [m/s²].  This is indeed an
// acceleration — a rate of change of velocity per unit time.
//
// In our simulator, this term is handled by semi-Lagrangian advection
// (advectVelCS): instead of computing the derivative explicitly, we
// trace each grid point backward along the velocity field and copy the
// velocity from the departure point — an equivalent but more stable
// approach (see the advectVelCS comment for details).
//
// ── Why is −(1/ρ) · ∇p an acceleration? ────────────────────────────────
// Consider a small cube of fluid with volume V and density ρ (mass = ρV).
// If the pressure on the left face is higher than on the right face,
// there is a net force pushing the cube to the right.  In general, the
// net pressure force on the cube is:
//
//   F_pressure = −∇p · V
//
// The minus sign is because force points from high pressure to low
// pressure — opposite to the gradient (which points toward increasing
// pressure).  Applying Newton's second law (F = m · a):
//
//   a = F / m = (−∇p · V) / (ρ · V) = −∇p / ρ = −(1/ρ) · ∇p
//
// The volume cancels out, giving an acceleration that depends only on
// the pressure gradient and the fluid density — not on the size of the
// fluid element.  This is the pressure acceleration term in (NS).
//
// ── What is operator splitting? ─────────────────────────────────────────
// Equation (NS) says the velocity changes due to ALL three effects
// simultaneously.  Solving them all coupled together in one step would
// be extremely complex.  Operator splitting is a strategy that breaks
// the problem into simpler sub-problems by handling one term at a time.
//
// Start from the discrete form of (NS).  If we discretize the left
// side with forward Euler, the full equation for one time step is:
//
//   (v^{n+1} − v^n) / Δt  =  A(v^n) + B(v^n) + C(v^n)
//
// where A, B, C are the three acceleration terms (advection, pressure,
// forces).  This says: the velocity change over Δt equals the total
// acceleration times Δt.  Rearranging:
//
//   v^{n+1} = v^n  +  Δt · [A(v^n) + B(v^n) + C(v^n)]
//
// Solving this directly would require evaluating all terms at once and
// dealing with their coupling (especially pressure, which involves a
// global constraint).  Operator splitting approximates this by solving
// each term sequentially, using the output of each sub-step as input
// to the next:
//
//   v^(1)   = v^n   + Δt · C(v^n)          [forces: gravity, etc.]
//   v^(2)   = v^(1) + Δt · B(v^(1))        [pressure correction]
//   v^{n+1} = v^(2) + Δt · A(v^(2))        [advection]
//
// Why does this work?  Each sub-step has the structure:
//
//   v_out = v_in + Δt · acceleration(v_in)
//
// This is just forward Euler applied to ∂v/∂t = acceleration(v).
// Each sub-step solves a simpler PDE where only ONE physical effect
// acts on the velocity.  The key insight is that if you ADD UP the
// velocity changes from all sub-steps, you get approximately the same
// total change as the original equation:
//
//   v^(1) − v^n   = Δt · C(v^n)            [change from forces]
//   v^(2) − v^(1) = Δt · B(v^(1))          [change from pressure]
//   v^{n+1}−v^(2) = Δt · A(v^(2))          [change from advection]
//   ─────────────────────────────────────
//   v^{n+1}− v^n  = Δt·C(v^n) + Δt·B(v^(1)) + Δt·A(v^(2))
//
// Compare with the unsplit version:
//
//   v^{n+1}− v^n  = Δt·C(v^n) + Δt·B(v^n) + Δt·A(v^n)
//
// The difference is that in the split version, B sees v^(1) instead of
// v^n, and A sees v^(2) instead of v^n.  Since v^(1) and v^(2) differ
// from v^n by O(Δt) (a small amount), these substitutions introduce an
// error of O(Δt²) per step — the same order as the Euler discretization
// itself.  So splitting doesn't make the approximation meaningfully
// worse, while making each sub-step dramatically simpler to solve.
//
// ── How our pipeline uses operator splitting ────────────────────────────
// Our simulation (see fluid_sim.cpp, Steps 1-15) splits (NS) into four
// sub-steps per frame, in this order:
//
//   STEP 2: applyForcesCS     v^(1) = v^n + f · Δt
//           Apply external forces (gravity, buoyancy, smoke injection).
//           This is just the f term from (NS), integrated with forward
//           Euler.
//
//   STEP 3-4: vorticity       v^(2) = v^(1) + f_vort · Δt
//           Apply vorticity confinement (an artificial force that
//           compensates for numerical dissipation of turbulence).
//           This is another external force, handled as a separate
//           sub-step.
//
//   STEP 6-8: pressure        v^(3) = v^(2) − (Δt/ρ) · ∇p
//           Divergence → pressure solve → project.
//           This handles the −(1/ρ)·∇p term from (NS).  The pressure is
//           chosen so that v^(3) is divergence-free (∇·v^(3) = 0).
//           This sub-step is equation (A) — what THIS kernel computes.
//           Derivation of (A) is explained in the comments below.
//
//   STEP 10: advectVelCS      v^{n+1} = advect(v^(3))
//           Semi-Lagrangian advection.  This handles the −(v·∇)v term
//           from (NS) by tracing particles backward through the
//           (now divergence-free) velocity field and copying the
//           velocity from the departure point.
//
// Note: projection happens BEFORE advection in our pipeline.  This
// means advection transports an already-divergence-free field, which
// is the ordering used by Stam ("Stable Fluids", 1999).  Advection
// can introduce small divergence due to numerical errors, so a final
// boundary pass (Step 15) cleans up afterward.
//
// ── Deriving equation (A) from the pressure sub-step ────────────────────
// As explained above, the pressure sub-step solves:
//
//   ∂v/∂t = −(1/ρ) · ∇p
//
// Discretizing with forward Euler over one time step Δt:
//
//   (v^{n+1} − v*) / Δt  =  −(1/ρ) · ∇p
//
// where v* is the velocity entering this sub-step (after forces and
// vorticity).  Rearranging:
//
//   v^{n+1} = v*  −  (Δt / ρ) · ∇p                                 (A)
//
// This is equation (A).  The pressure p is not arbitrary: the pressure
// solver (see pressureRedCS/pressureBlackCS, STEP 1-7) computed exactly
// the p that makes v^{n+1} divergence-free (∇ · v^{n+1} = 0),
// ensuring incompressibility.
//
// ── Discretizing the pressure gradient ──────────────────────────────────
// On a MAC grid, each velocity component lives on a cell face.  We need
// the pressure gradient at that face.  Here's why it's a simple
// difference divided by h:
//
// The pressure gradient in the x-direction is the derivative:
//
//   (∇p)_x  =  ∂p/∂x
//
// The standard finite-difference approximation of a first derivative is
// the central difference formula:
//
//   ∂p/∂x  ≈  (p(x + h/2) − p(x − h/2)) / h
//
// We evaluate the function at two points symmetric around x, separated
// by h, and divide by that distance.
//
// On a MAC grid, the u-velocity at face (i, j, k) sits exactly at the
// boundary between cell (i−1) and cell (i).  The cell centers — where
// pressure is stored — are at positions (i−1)·h + h/2 and i·h + h/2,
// which are exactly h/2 to the left and h/2 to the right of the face.
// So the central difference formula becomes:
//
//   (∇p)_x at face i  ≈  (p_{i,j,k} − p_{i−1,j,k}) / h
//
//              p_{i-1}           p_{i}
//                 ·                ·          ← cell centers (pressure)
//         |       :       |       :       |
//         |   cell i-1    |    cell i     |
//         |       :       |       :       |
//                         ↑
//                    u-face [i]
//                 (velocity lives here)
//
//              |←  h/2  →|←  h/2  →|
//              |←──────  h  ──────→|
//
// The two pressure samples are exactly h apart, centered on the face.
// This is second-order accurate — no approximation beyond the usual
// finite-difference truncation.  The same logic applies to v-faces (Y)
// and w-faces (Z).
//
// Substituting into (A) and precomputing  scale = Δt / (ρ · h):
//
//   u[i,j,k]^{n+1} = u[i,j,k]*  −  scale · (p_{i,j,k} − p_{i−1,j,k})
//
// which is what the code computes as:
//   texU[i,j,k] -= scale * (p − p_left)
//
// ── Owned-face approach (avoiding race conditions) ──────────────────────
// Each face on the MAC grid is shared between two adjacent cells.  If
// both cells tried to update the same face, we'd have a race condition
// (two threads doing a non-atomic read-modify-write to the same address).
//
// Solution: each cell "owns" only its 3 LOWER-index faces:
//   - u[i]  (left face in X)
//   - v[j]  (bottom face in Y)
//   - w[k]  (back face in Z)
//
// The 3 UPPER-index faces (u[i+1], v[j+1], w[k+1]) are owned by the
// neighboring cell in that direction — it will update them as its own
// lower-index face.  This guarantees every face is written by exactly
// one thread.
//
// Faces at the domain boundary (adjacent to solid) are handled by
// boundaryCS instead.
// ============================================================================
#ifdef PROJECT
[numthreads(FLUID_THREADS_3D_X, FLUID_THREADS_3D_Y, FLUID_THREADS_3D_Z)]
void projectCS(uint3 dtid : SV_DispatchThreadID)
{
	// Get cell indices from thread ID
    int i = (int)dtid.x;
    int j = (int)dtid.y;
    int k = (int)dtid.z;

	// Bounds check: skip threads outside the grid dimensions.
    if (i >= cb.gridX || j >= cb.gridY || k >= cb.gridZ)
        return;

    if (texSolid[int3(i, j, k)] == 0)
        return;  // skip solid

    float h = cb.cellSize;
    // scale = Δt / (ρ · h), the coefficient from equation (A) with the
    // 1/h from the discrete gradient folded in.
    float scale = cb.dt / (cb.density * h);
    float p = texPressure[int3(i, j, k)];

    // Left face (u[i,j,k]) — owned by this cell.
    // Only correct if left neighbor is fluid; solid faces are set by boundaryCS.
    if (i > 0 && texSolid[int3(i - 1, j, k)] != 0)
    {
        texU[int3(i, j, k)] -= scale * (p - texPressure[int3(i - 1, j, k)]);
    }

    // Bottom face (v[i,j,k]) — owned by this cell.
    if (j > 0 && texSolid[int3(i, j - 1, k)] != 0)
    {
        texV[int3(i, j, k)] -= scale * (p - texPressure[int3(i, j - 1, k)]);
    }

    // Back face (w[i,j,k]) — owned by this cell.
    if (k > 0 && texSolid[int3(i, j, k - 1)] != 0)
    {
        texW[int3(i, j, k)] -= scale * (p - texPressure[int3(i, j, k - 1)]);
    }

    // The right (u[i+1]), top (v[j+1]), and front (w[k+1]) faces are
    // owned by the neighboring cell in that direction — it will update
    // them as its own lower-index face.
}
#endif

// ============================================================================
// KERNEL: advectVelCS
// Semi-Lagrangian advection for velocity.
// For each velocity component, backtrace from its grid location using the
// averaged velocity field, then sample the old field at the departure point.
//
// ── What this kernel does ───────────────────────────────────────────────
// This kernel applies the advection part of the Navier-Stokes update to
// the velocity field itself.  The continuous meaning of the advection term
// −(v · ∇)v is explained in detail in projectCS; here we show how that term
// is implemented numerically with a Semi-Lagrangian scheme.
//
// People sometimes hear "self-advection" and assume it is something exotic.
// It is not.  The idea is the same as for a particle carrying any other
// property while it moves: a parcel of fluid moves with velocity v, and it
// carries that same velocity with it.  So to update the velocity stored at
// a grid location x at time t + Δt, we ask:
//
//   "Which parcel of fluid arrived here?"
//   "Where was it one timestep ago?"
//   "What velocity was it carrying there?"
//
// The advection-only equation is:
//
//   ∂v/∂t + (v · ∇)v = 0
//
// or equivalently:
//
//   Dv/Dt = 0
//
// meaning that, along a fluid parcel trajectory, the parcel keeps its
// velocity while it is transported by the flow.  This leads directly to the
// Semi-Lagrangian update:
//
//   x_depart = x - v(x) · Δt
//   v_new(x) = v_old(x_depart)
//
// So the whole kernel is just:
//   1. start from the grid-index-space location of a stored velocity sample,
//   2. trace backward through the flow,
//   3. read the old velocity at that departure point,
//   4. write it as the new velocity here.
//
// In practice x_depart almost never lands exactly on a stored sample, so
// we interpolate the old field at that position (sampleVelocity).
//
// ── Why backtrace instead of using finite differences directly? ─────────
// Directly discretizing (v · ∇)v with an explicit finite-difference scheme
// is possible, but it is much more sensitive to the timestep and can become
// unstable easily.  The Semi-Lagrangian approach is far more robust: it
// asks where the arriving fluid came from and copies the old value from
// there, which remains stable even for relatively large Δt.
//
// The price is numerical diffusion.  Because the departure point usually
// lies between grid samples, we must interpolate, and interpolation is a
// weighted average of neighboring values.  Repeating that process frame
// after frame smooths sharp velocity gradients and gradually weakens
// vortices.  See applyVorticityCS for a more detailed explanation of why
// this smoothing attenuates rotational motion and why vorticity confinement
// is added afterward.
//
// ── Why velocity advection is trickier than smoke advection ─────────────
// Smoke and temperature are simpler because they are cell-centered scalar
// fields: for each cell, we backtrace from the cell center in grid-index
// space and sample one
// scalar value.
//
// Velocity is more involved for two reasons:
//
// 1. The backtrace does NOT start from the cell center in grid-index space.
//    On a MAC grid, each velocity component is stored at a face center:
//      - u lives on X faces: (i,     j+0.5, k+0.5)
//      - v lives on Y faces: (i+0.5, j,     k+0.5)
//      - w lives on Z faces: (i+0.5, j+0.5, k)
//    So u, v, and w must each be advected from their own staggered location.
//
// 2. Backtracing requires the FULL 3D velocity at that location.
//    When advecting u, the u component is already stored exactly at the
//    current u-face, but v and w are not.  They must be reconstructed there
//    by averaging nearby staggered samples.  The same logic applies when
//    advecting v and w.
//
// Key insight for staggered grids:
//   - u lives at (i, j+0.5, k+0.5) → to advect u, we need v and w interpolated at u's location
//   - v lives at (i+0.5, j, k+0.5) → needs u and w interpolated at v's location
//   - w lives at (i+0.5, j+0.5, k) → needs u and v interpolated at w's location
//
// So, compared to smoke/temperature, velocity advection is harder not
// because the physical idea changes, but because the MAC layout stores the
// three scalar components at different spatial locations.
//
// ── Solid handling and boundary cleanup ────────────────────────────────
// A face is advected only if at least one of the two cells sharing that
// face is fluid.  Faces surrounded by solid are left unchanged here.
// After this pass, boundaryCS runs again because Semi-Lagrangian sampling
// near solids can produce face values that no longer satisfy the imposed
// no-penetration / free-slip boundary conditions exactly.
//
// ── Implementation notes ───────────────────────────────────────────────
// 1. Read the current velocity field from texU/texV/texW.
// 2. For each staggered face sample, reconstruct the full local velocity.
// 3. Backtrace that face position by Δt.
// 4. Sample the OLD velocity field at the departure point.
// 5. Write the result into texU_temp / texV_temp / texW_temp.
// 6. Apply a small multiplicative dissipation (cb.velDissipation).
//
// This kernel runs on a staggered MAC grid:
//  u lives on X faces and has (gridX + 1, gridY,     gridZ    ) samples,
//  v lives on Y faces and has (gridX,     gridY + 1, gridZ    ) samples,
//  w lives on Z faces and has (gridX,     gridY,     gridZ + 1) samples.
// This is why the three IF bounds differ: only the index along the
// component's own axis is allowed to reach the corresponding grid dimension.
//
// To backtrace any one component correctly, we need the full velocity vector
// at that component's own staggered face location. Example: when advecting u,
// the current sample position is
//  p = (i, j+0.5, k+0.5).
// The u component is stored exactly there, so
//  vel_u = texU[i,j,k]
// is a direct read. But v and w are stored on different staggered faces, so
// p does not contain direct v or w samples. We therefore cannot use
// texV[i,j,k] or texW[i,j,k] directly, because those indices refer to valid
// v- and w-sample locations in their own staggered textures, not to the
// u-face position p itself.
//
// Instead, we reconstruct the missing cross-components from nearby staggered
// sample positions around p. For vel_v at a u-face, we evaluate the v field
// at the 4 nearby Y-face positions surrounding that u-face:
//  (i-0.5, j,   k+0.5), (i+0.5, j,   k+0.5),
//  (i-0.5, j+1, k+0.5), (i+0.5, j+1, k+0.5)
// and average those 4 values. These are positions in grid-index space derived
// from p = (i, j+0.5, k+0.5). They are valid v-sample positions because v
// lives on Y-faces: its y coordinate is integer, while x and z are
// half-integer. For vel_w, we do the analogous thing with the 4 nearby
// Z-face positions around the same u-face.
//
// However, those 4 positions are still continuous positions in grid-index space,
// not integer texel indices into texV or texW. We therefore cannot read the
// textures directly at those coordinates, nor can we just pick an arbitrary
// nearest texel without introducing a crude nearest-neighbor approximation.
// Instead, each sampleVelocity(pos, component) call performs a local trilinear
// interpolation of that single scalar component at the requested position: it
// finds the 8 neighboring texels of the requested staggered field around that
// position and blends them to return the component value there. The outer
// average of 4 calls then combines those 4 reconstructed neighboring values
// to estimate the missing cross-component at p.
//
// Once the full local velocity is known, we backtrace the current face
// position by Δt, sample the old component field at that departure point, and
// write the result into the matching temp texture. Before writing, we multiply
// by velDissipation, which applies a small per-step damping to the velocity
// field. This slightly reduces kinetic energy over time, helping prevent tiny
// numerical oscillations from persisting indefinitely and providing simple
// control over how quickly motion dies out. A value of 1.0 means no extra
// damping.
//
// The temp textures are used for ping-ponging so that all threads read the
// old velocity field consistently while writing the advected one.
// ============================================================================
#ifdef ADVECT_VELOCITY
[numthreads(FLUID_THREADS_3D_X, FLUID_THREADS_3D_Y, FLUID_THREADS_3D_Z)]
void advectVelCS(uint3 dtid : SV_DispatchThreadID)
{
	// Get cell size and timestep from constant buffer
    float h = cb.cellSize;
    float dt = cb.dt;

	// Get cell indices from thread ID
    int i = (int)dtid.x;
    int j = (int)dtid.y;
    int k = (int)dtid.z;

    // --- Advect U component ---
    // u lives at face (i, j+0.5, k+0.5) in grid-index space
	//
	// If this face is adjacent to fluid (not solid) on at least one side, we backtrace and update it.
    if (i <= cb.gridX && j < cb.gridY && k < cb.gridZ)
    {
        // Skip faces touching solid cells on both sides
        bool leftFluid  = (i > 0 && i <= cb.gridX && isInBounds(i - 1, j, k)) ? (texSolid[int3(i - 1, j, k)] != 0) : false;
        bool rightFluid = (i < cb.gridX && isInBounds(i, j, k)) ? (texSolid[int3(i, j, k)] != 0) : false;

        if (leftFluid || rightFluid)
        {
            // Position of this u sample in grid-index space
            float px = (float)i;
            float py = (float)j + 0.5;
            float pz = (float)k + 0.5;

            // Get velocity at this location
            float vel_u = texU[int3(i, j, k)];
            // Average v at u's location: average of 4 neighboring v values
            float vel_v = 0.25 * (
                sampleVelocity(float3(px - 0.5, py - 0.5, pz), 1) +
                sampleVelocity(float3(px + 0.5, py - 0.5, pz), 1) +
                sampleVelocity(float3(px - 0.5, py + 0.5, pz), 1) +
                sampleVelocity(float3(px + 0.5, py + 0.5, pz), 1));
            // Average w at u's location (4 nearest w-face values, offset in x and z)
            float vel_w = 0.25 * (
                sampleVelocity(float3(px - 0.5, py, pz - 0.5), 2) +
                sampleVelocity(float3(px + 0.5, py, pz - 0.5), 2) +
                sampleVelocity(float3(px - 0.5, py, pz + 0.5), 2) +
                sampleVelocity(float3(px + 0.5, py, pz + 0.5), 2));

            // Backtrace
            float bx = px - vel_u * dt / h;
            float by = py - vel_v * dt / h;
            float bz = pz - vel_w * dt / h;

            // Clamp to valid range
            bx = clamp(bx, 0.0, (float)cb.gridX);
            by = clamp(by, 0.5, (float)cb.gridY - 0.5);
            bz = clamp(bz, 0.5, (float)cb.gridZ - 0.5);

			// Sample u at the backtraced position, apply dissipation,
			// and write the result to the temp texture at the current face.
			// This is the Semi-Lagrangian update: new u at the current face
			// is the old u sampled at the departure point.
			// This is the core of the advection step for velocity.
            texU_temp[int3(i, j, k)] = sampleVelocity(float3(bx, by, bz), 0) * cb.velDissipation;
        }
        else // otherwise, if this face is surrounded by solid, we only apply dissipation to the existing value without backtracing.
        {
            texU_temp[int3(i, j, k)] = texU[int3(i, j, k)] * cb.velDissipation;
        }
    }

    // --- Advect V component ---
    // v lives at face (i+0.5, j, k+0.5) in grid-index space
    if (i < cb.gridX && j <= cb.gridY && k < cb.gridZ)
    {
        bool belowFluid = (j > 0 && isInBounds(i, j - 1, k)) ? (texSolid[int3(i, j - 1, k)] != 0) : false;
        bool aboveFluid = (j < cb.gridY && isInBounds(i, j, k)) ? (texSolid[int3(i, j, k)] != 0) : false;

        if (belowFluid || aboveFluid)
        {
            float px = (float)i + 0.5;
            float py = (float)j;
            float pz = (float)k + 0.5;

            float vel_v = texV[int3(i, j, k)];
            // Average u at v's location
            float vel_u = 0.25 * (
                sampleVelocity(float3(px - 0.5, py - 0.5, pz), 0) +
                sampleVelocity(float3(px + 0.5, py - 0.5, pz), 0) +
                sampleVelocity(float3(px - 0.5, py + 0.5, pz), 0) +
                sampleVelocity(float3(px + 0.5, py + 0.5, pz), 0));
            float vel_w = 0.25 * (
                sampleVelocity(float3(px, py - 0.5, pz - 0.5), 2) +
                sampleVelocity(float3(px, py + 0.5, pz - 0.5), 2) +
                sampleVelocity(float3(px, py - 0.5, pz + 0.5), 2) +
                sampleVelocity(float3(px, py + 0.5, pz + 0.5), 2));

            float bx = px - vel_u * dt / h;
            float by = py - vel_v * dt / h;
            float bz = pz - vel_w * dt / h;

            bx = clamp(bx, 0.5, (float)cb.gridX - 0.5);
            by = clamp(by, 0.0, (float)cb.gridY);
            bz = clamp(bz, 0.5, (float)cb.gridZ - 0.5);

            texV_temp[int3(i, j, k)] = sampleVelocity(float3(bx, by, bz), 1) * cb.velDissipation;
        }
        else
        {
            texV_temp[int3(i, j, k)] = texV[int3(i, j, k)] * cb.velDissipation;
        }
    }

    // --- Advect W component ---
    // w lives at face (i+0.5, j+0.5, k) in grid-index space
    if (i < cb.gridX && j < cb.gridY && k <= cb.gridZ)
    {
        bool backFluid  = (k > 0 && isInBounds(i, j, k - 1)) ? (texSolid[int3(i, j, k - 1)] != 0) : false;
        bool frontFluid = (k < cb.gridZ && isInBounds(i, j, k)) ? (texSolid[int3(i, j, k)] != 0) : false;

        if (backFluid || frontFluid)
        {
            float px = (float)i + 0.5;
            float py = (float)j + 0.5;
            float pz = (float)k;

            float vel_w = texW[int3(i, j, k)];
            float vel_u = 0.25 * (
                sampleVelocity(float3(px - 0.5, py, pz - 0.5), 0) +
                sampleVelocity(float3(px + 0.5, py, pz - 0.5), 0) +
                sampleVelocity(float3(px - 0.5, py, pz + 0.5), 0) +
                sampleVelocity(float3(px + 0.5, py, pz + 0.5), 0));
            float vel_v = 0.25 * (
                sampleVelocity(float3(px, py - 0.5, pz - 0.5), 1) +
                sampleVelocity(float3(px, py + 0.5, pz - 0.5), 1) +
                sampleVelocity(float3(px, py - 0.5, pz + 0.5), 1) +
                sampleVelocity(float3(px, py + 0.5, pz + 0.5), 1));

            float bx = px - vel_u * dt / h;
            float by = py - vel_v * dt / h;
            float bz = pz - vel_w * dt / h;

            bx = clamp(bx, 0.5, (float)cb.gridX - 0.5);
            by = clamp(by, 0.5, (float)cb.gridY - 0.5);
            bz = clamp(bz, 0.0, (float)cb.gridZ);

            texW_temp[int3(i, j, k)] = sampleVelocity(float3(bx, by, bz), 2) * cb.velDissipation;
        }
        else
        {
            texW_temp[int3(i, j, k)] = texW[int3(i, j, k)] * cb.velDissipation;
        }
    }
}
#endif

// ============================================================================
// KERNEL: advectSmokeCS
// Semi-Lagrangian advection for smoke density.
//
// Smoke is a cell-centered scalar field, so this kernel is the simpler scalar
// counterpart of advectVelCS.  The physical idea is the same:
//
//   "Which parcel of fluid arrived at this cell center?"
//   "Where was that parcel one timestep ago?"
//   "How much smoke did it carry there?"
//
// The advection-only equation for smoke density s is:
//
//   ∂s/∂t + (v · ∇)s = 0
//
// or equivalently:
//
//   Ds/Dt = 0
//
// meaning that a fluid parcel keeps its smoke value while it is transported
// by the flow.  This gives the Semi-Lagrangian update:
//
//   x_depart = x - v(x) · Δt
//   s_new(x) = s_old(x_depart)
//
// For smoke, x is the cell center, so for cell (i,j,k) we start from:
//
//   x = (i+0.5, j+0.5, k+0.5)
//
// Unlike velocity advection, we do not have to advect three staggered
// component fields separately.  We only need:
// 1. the smoke value stored at cell centers, and
// 2. the velocity at the same cell center.
//
// The velocity at the cell center is reconstructed by averaging opposite MAC
// faces (avgVelAtCenter):
//   u = 0.5 * (u_left  + u_right)
//   v = 0.5 * (v_down  + v_up)
//   w = 0.5 * (w_back  + w_front)
//
// Once that center velocity is known, we backtrace the center by Δt.
//
// ── The important coordinate-convention detail ───────────────────────────
// advectSmokeCS computes the departure point in cell-center coordinates:
//
//   px = i + 0.5
//   py = j + 0.5
//   pz = k + 0.5
//
// so the backtraced position (bx,by,bz) is also in cell-center coordinates.
//
// But sampleSmoke() uses a different convention: it treats texSmoke[i,j,k] as
// if it were located at coordinate (i,j,k), not at (i+0.5,j+0.5,k+0.5).
// Therefore, before calling sampleSmoke(), we must convert from
// cell-center coordinates to the coordinates expected by that helper:
//
//   samplePos = (bx-0.5, by-0.5, bz-0.5)
//
// Numerical example in 1D:
//   suppose we are updating cell i = 3, whose center is at px = 3.5
//   and backtracing gives bx = 2.8.
//
// Physically, x = 2.8 lies between cell-center 2.5 (cell 2) and
// cell-center 3.5 (cell 3), 30% of the way from cell 2 toward cell 3.
//
// To sample correctly we call:
//   sampleSmoke(2.8 - 0.5) = sampleSmoke(2.3)
//
// Then sampleSmoke computes:
//   i0 = floor(2.3) = 2
//   i1 = 3
//   f  = 0.3
//
// and returns:
//   0.7 * smoke[2] + 0.3 * smoke[3]
//
// which is exactly the desired interpolation at physical position 2.8.
// Without the -0.5 conversion, we would instead sample at 2.8 in
// texture-index space, which corresponds to the wrong physical location and
// introduces a systematic half-cell shift.
//
// As in all Semi-Lagrangian advection, the departure point rarely lands
// exactly on a stored sample, so interpolation is necessary.  Repeated
// interpolation makes the method stable, but also numerically diffusive:
// sharp smoke features gradually get smoother over time.
//
// Implementation summary:
// 1. Skip cells outside the domain.
// 2. Force smoke to zero in solid cells.
// 3. Compute the current cell-center position.
// 4. Reconstruct the center velocity with avgVelAtCenter().
// 5. Backtrace by Δt.
// 6. Convert from cell-center coordinates to sampleSmoke coordinates.
// 7. Sample old smoke at the departure point.
// 8. Apply smoke dissipation and write to texSmoke_temp.
// ============================================================================
#ifdef ADVECT_SMOKE
[numthreads(FLUID_THREADS_3D_X, FLUID_THREADS_3D_Y, FLUID_THREADS_3D_Z)]
void advectSmokeCS(uint3 dtid : SV_DispatchThreadID)
{
	// Get cell indices from thread ID
    int i = (int)dtid.x;
    int j = (int)dtid.y;
    int k = (int)dtid.z;

	// Bounds check: skip threads outside the grid dimensions.
    if (i >= cb.gridX || j >= cb.gridY || k >= cb.gridZ)
        return;

	// If this cell is solid...
    if (texSolid[int3(i, j, k)] == 0)
    {
		// ...we don't backtrace or sample; we just set smoke to zero here.
		// This ensures that smoke cannot exist in solid cells.
        texSmoke_temp[int3(i, j, k)] = 0.0;
        return;
    }

	// Get cell size and timestep from constant buffer
    float h = cb.cellSize;
    float dt = cb.dt;

    // Smoke is at the cell center in grid-index space: (i+0.5, j+0.5, k+0.5)
    float px = (float)i + 0.5;
    float py = (float)j + 0.5;
    float pz = (float)k + 0.5;

    // Average velocity at cell center
    float3 vel = avgVelAtCenter(i, j, k);

    // Backtrace
    float bx = px - vel.x * dt / h;
    float by = py - vel.y * dt / h;
    float bz = pz - vel.z * dt / h;

    // Clamp in grid-index space (sampleSmoke also performs its own clamping)
    bx = clamp(bx, 0.5, (float)cb.gridX - 0.5);
    by = clamp(by, 0.5, (float)cb.gridY - 0.5);
    bz = clamp(bz, 0.5, (float)cb.gridZ - 0.5);

    // sampleSmoke() interprets texSmoke[i,j,k] as living at coordinate
    // (i,j,k), while bx/by/bz are still in cell-center coordinates
    // (i+0.5,j+0.5,k+0.5).  Shift by -0.5 to express the departure point in
    // the coordinate system expected by sampleSmoke().
	// In other words, we convert the backtraced position from cell-center
	// coordinates (i+0.5,j+0.5,k+0.5) back to the texture-index coordinates
	// used by sampleSmoke().
    float advected = sampleSmoke(float3(bx - 0.5, by - 0.5, bz - 0.5));

    // Dissipation: multiplicative decay (0.998 = keep 99.8% per step)
    advected *= cb.smokeDissipation;

    // Write the advected smoke value to the temp texture at the current cell.
    // This is the Semi-Lagrangian update: new smoke at the current cell =
    // old smoke sampled at the departure point.
    // This is the core of the advection step for smoke.
    texSmoke_temp[int3(i, j, k)] = advected;
}
#endif

// ============================================================================
// KERNEL: computeCurlCS
// Compute the curl (vorticity) ∇ × v  of the velocity field at each cell center.
// Uses cell-center averaged velocities and central differences.
// ============================================================================
#ifdef COMPUTE_CURL
[numthreads(FLUID_THREADS_3D_X, FLUID_THREADS_3D_Y, FLUID_THREADS_3D_Z)]
void computeCurlCS(uint3 dtid : SV_DispatchThreadID)
{
	// Get cell index from thread ID
    int i = (int)dtid.x;
    int j = (int)dtid.y;
    int k = (int)dtid.z;

	// Bounds check
    if (i >= cb.gridX || j >= cb.gridY || k >= cb.gridZ)
        return;

	// Skip solid cells (curl is only defined in fluid)
    if (texSolid[int3(i, j, k)] == 0)
    {
        texCurl[int3(i, j, k)] = float4(0, 0, 0, 0);
        return;
    }

    // ── Definition of curl (vorticity) ──────────────────────────────────
    // The curl of a 3D velocity field v = (u, v, w) is defined as:
    //
    //   ω = ∇ × v = ( ∂w/∂y − ∂v/∂z,
    //                 ∂u/∂z − ∂w/∂x,
    //                 ∂v/∂x − ∂u/∂y )
    //
    // Each component of ω measures local rotation around its namesake axis:
    //   ω.x  depends on v,w varying along y,z → rotation around X (YZ plane)
    //   ω.y  depends on u,w varying along z,x → rotation around Y (XZ plane)
    //   ω.z  depends on u,v varying along x,y → rotation around Z (XY plane)
    // Intuitively: look down the +X axis — ω.x tells you how much the fluid
    // swirls in the plane you see (YZ). Positive = counter-clockwise by the
    // right-hand rule.
    //
    // ── Central differences ─────────────────────────────────────────────
    // To approximate these partial derivatives on a uniform grid we use
    // *central differences*: for a generic field f along axis x,
    //
    //   ∂f/∂x ≈ ( f(i+1) − f(i−1) ) / (2·h)
    //
    // Central differences are second-order accurate (O(h²)), meaning the
    // error shrinks quadratically as the grid is refined — a significant
    // improvement over one-sided (forward/backward) differences which are
    // only first-order (O(h)).  They achieve this because the first-order
    // error terms cancel symmetrically around the evaluation point.
    //
    // Computing ∂f/∂x, ∂f/∂y, ∂f/∂z requires sampling f at the two
    // neighbors along each axis → 6 neighbors total:
    //   x-axis: (i−1,j,k) and (i+1,j,k)   → ∂/∂x
    //   y-axis: (i,j−1,k) and (i,j+1,k)   → ∂/∂y
    //   z-axis: (i,j,k−1) and (i,j,k+1)   → ∂/∂z
    //
    // ── Index-space vs. world-space (why h = 1, not cellSize) ───────────
    // The physical cell size is cb.cellSize = domainSize / gridX (e.g.
    // 4.0 / 128 = 0.03125 m).  The *true* world-space derivative would be:
    //
    //   ∂w/∂y_physical = (w(j+1) − w(j−1)) / (2 · cellSize)
    //
    // However, this kernel works in *grid-index coordinates*, where the
    // spacing between adjacent cells is 1 by definition.  So we compute:
    //
    //   ∂w/∂y_code = (w(j+1) − w(j−1)) / (2 · 1)
    //
    // The resulting curl is therefore scaled by cellSize relative to the
    // true physical curl:
    //
    //   ω_code = ω_physical · cellSize
    //
    // Numeric example: if (vT.z − vB.z) = 10.0 and cellSize = 0.03125:
    //   ω_physical = 10 / (2 · 0.03125) = 160.0   (true vorticity, 1/s)
    //   ω_code     = 10 / (2 · 1)       = 5.0     (what we compute here)
    //   check: 5.0 = 160.0 · 0.03125              ✓
    //
    // This is NOT a bug.  It is intentional because the only consumer of
    // this curl is applyVorticityCS, which uses the Fedkiw vorticity
    // confinement formula.  That formula has an explicit grid-spacing
    // factor h (= cellSize):
    //
    //   f = ε · cellSize · (η̂ × ω_physical)        [Fedkiw et al.]
    //
    // Since ω_code already contains ω_physical · cellSize, the code in
    // applyVorticityCS can simply write:
    //
    //   f = ε · (η̂ × ω_code)                       [what the code does]
    //     = ε · (η̂ × ω_physical · cellSize)
    //     = ε · cellSize · (η̂ × ω_physical)         [= Fedkiw formula ✓]
    //
    // The "missing" cellSize in the curl denominator here is exactly the h
    // factor that Fedkiw's formula needs.  The two cancel out perfectly.
    // See applyVorticityCS below for the full step-by-step derivation.
    // ─────────────────────────────────────────────────────────────────────

    // Sample cell-center velocity at the 6 axis-aligned neighbors.
    // Indices are clamped to grid bounds so boundary cells gracefully
    // degenerate to one-sided differences (the stencil shrinks but the
    // code path stays branchless).
    int il = max(i - 1, 0), ir = min(i + 1, cb.gridX - 1);
    int jb = max(j - 1, 0), jt = min(j + 1, cb.gridY - 1);
    int kb = max(k - 1, 0), kf = min(k + 1, cb.gridZ - 1);

    float3 vL  = avgVelAtCenter(il, j, k);  // x−1  (left)
    float3 vR  = avgVelAtCenter(ir, j, k);  // x+1  (right)
    float3 vB  = avgVelAtCenter(i, jb, k);  // y−1  (bottom)
    float3 vT  = avgVelAtCenter(i, jt, k);  // y+1  (top)
    float3 vBk = avgVelAtCenter(i, j, kb);  // z−1  (back)
    float3 vF  = avgVelAtCenter(i, j, kf);  // z+1  (front)

    // Apply the curl formula with central differences:
    //   ω.x = ∂w/∂y − ∂v/∂z  ≈ (vT.z − vB.z) − (vF.y − vBk.y)
    //   ω.y = ∂u/∂z − ∂w/∂x  ≈ (vF.x − vBk.x) − (vR.z − vL.z)
    //   ω.z = ∂v/∂x − ∂u/∂y  ≈ (vR.y − vL.y) − (vT.x − vB.x)
	//
	// Here curl is ω_code = ω_physical · cellSize, so the central
	// difference denominator is 2·1 = 2 (see note above).
    float3 curl;
    curl.x = (vT.z - vB.z) - (vF.y - vBk.y);
    curl.y = (vF.x - vBk.x) - (vR.z - vL.z);
    curl.z = (vR.y - vL.y) - (vT.x - vB.x);
    curl *= 0.5;  // 1 / (2·h) with h = 1 in index space (see note above)

    texCurl[int3(i, j, k)] = float4(curl, 0.0);
}
#endif

// ============================================================================
// KERNEL: applyVorticityCS
// Vorticity confinement — re-inject the rotational energy that numerical
// dissipation silently destroys every frame.
//
// Important clarification:
//   this pass acts only on the VELOCITY field; it does not affect smoke.
//   It is also not the same thing as MacCormack correction.
//
//   - applyVorticityCS adds an explicit confinement force that strengthens
//     existing swirling motion (curl / vortices).
//   - macCormackVelCS, later in the pipeline, corrects the numerical
//     advection error of the velocity field itself.
//
// So:
//   vorticity confinement = add rotational energy back where vortices exist
//   MacCormack velocity   = reduce Semi-Lagrangian diffusion of velocity
//
// ── Why is this needed? ─────────────────────────────────────────────────
// Real fluids conserve vorticity almost perfectly: smoke rings, turbulent
// eddies, and swirling plumes persist for a long time.  Our simulator,
// however, uses Semi-Lagrangian (SL) advection, which introduces
// *numerical diffusion* that quietly kills vortices (see advectVelCS above).
// Here is why, step by step:
//
// 1. SL advection starts at the grid position of the quantity being
//    advected. On our staggered grid, each velocity component lives
//    on its own cell face — not at the cell center like smoke and
//    temperature.
//
// 2. From that starting position, we *backtrace* by −v·dt to find
//    where the fluid at this grid point "came from" one timestep ago.
//    The velocity used for backtracing is interpolated at the starting
//    position (since on a staggered grid the other components are
//    stored at different locations).
//
// 3. The backtraced position almost never lands exactly on a grid node.
//    It falls *between* nodes, so we must interpolate the *old* field
//    (the same field being advected) at that point — trilinear in 3D.
//
// 4. That interpolated value is assigned as the new value at the
//    starting position.  In code (advectVelCS):
//      texU_temp[i,j,k] = sampleVelocity(backtraced_pos, 0);
//
// 5. The problem: interpolation is mathematically a weighted average
//    of neighboring values.  Averaging is a smoothing operation — it
//    pulls extreme values toward the local mean.
//
// 6. Concrete 1D example — a sharp velocity spike (the simplest analog
//    of a vortex's steep gradient):
//
//      Before advection:
//        Cell:  0   1   2   3   4   5   6
//        Vel:   0   0   0  10   0   0   0
//
//      The backtrace from cells 2 and 4 lands between cells that
//      straddle the spike → linear interpolation blends the peak with
//      its zero neighbors:
//
//      After one SL step (assuming backtrace offset ~0.3 cells):
//        Vel:   0   0   3   7   3   0   0
//                       ↑   ↑   ↑
//             leaked  peak  leaked
//             from 3  fell  from 3
//
//      The peak dropped from 10 → 7 and energy leaked to neighbors.
//      This is *numerical diffusion*: the interpolation behaves like
//      a physical viscosity that was never intended.
//
// 7. Each advection step repeats this smoothing.  High-frequency features
//    (sharp gradients = vortex boundaries) are attenuated the most,
//    because their variation is on the scale of a single cell — exactly
//    where interpolation error is largest.
//
// 8. After many frames, the velocity field becomes increasingly smooth:
//    vortices widen, weaken, and eventually vanish, leaving the flow
//    looking laminar and "blobby" — completely lacking the turbulent
//    detail real smoke exhibits.
//
// Vorticity confinement (below) counteracts this by detecting where
// rotation still survives and nudging the velocity field to reinforce it.
//
// ── What does vorticity confinement do? ─────────────────────────────────
// Fedkiw et al. (2001) proposed adding a small corrective acceleration
// that reinforces surviving vortices to compensate for numerical loss.
// The final velocity update (via operator splitting of Navier-Stokes) is:
//
//   Δv = ε · h · (η̂ × ω) · dt                                [Fedkiw et al.]
//
// where h = cellSize.  The next section derives this formula from scratch;
// the glossary after that explains every symbol in geometric detail.
//
// ── Deriving the Fedkiw formula ────────────────────────────────────────
//
// STEP 1: THE STARTING POINT — NAVIER-STOKES WITH AN EXTRA FORCE
//
//   Vorticity confinement is modeled as an additional body force f_conf
//   in the Navier-Stokes momentum equation (see projectCS, eq. (NS)):
//
//     ∂v/∂t = ... + f_conf
//
//   Our pipeline handles it via operator splitting (see projectCS,
//   "How our pipeline uses operator splitting", STEP 3-4):
//
//     v^(2) = v^(1) + f_conf · Δt
//
//   So the goal is to define f_conf — an acceleration (m/s²) that
//   reinforces existing rotation without inventing new vortices.
//
// STEP 2: WHAT DIRECTION SHOULD THE FORCE POINT?
//
//   We need the force to push fluid in the direction it is already
//   spinning.  Two pieces of information identify that direction:
//
//   (a) ω = ∇ × v  (the vorticity vector, computed by computeCurlCS).
//       Points along the vortex AXIS (right-hand rule: curl your
//       fingers in the spinning direction → your thumb = ω).
//
//   (b) ∇|ω|  (gradient of the vorticity magnitude).
//       |ω| peaks at vortex cores and decays outward, so ∇|ω|
//       points RADIALLY INWARD toward the nearest vortex center.
//
//   Define the unit radial vector:  η̂ = normalize(∇|ω|)
//
//   Now take the cross product η̂ × ω.  This combines a radial vector
//   (η̂, toward the core) with an axial vector (ω, along the vortex
//   axis).  The result is perpendicular to both — that is, TANGENTIAL
//   to the rotation.  This is exactly the direction the fluid is
//   already spinning:
//
//        ω (axis)
//        ↑
//        |    ╭───╮
//        |   ╱     ╲  ← η̂ × ω is tangential to this spin.
//     ───●──→ η̂    │    It reinforces the existing rotation
//        |   ╲     ╱    without pushing fluid toward/away from
//        |    ╰───╯     the core or along the axis.
//        (radial, toward core)
//
//   So the shape of the force must be:  f_conf ∝ η̂ × ω
//
// STEP 3: WHAT MAGNITUDE SHOULD IT HAVE?
//
//   The cross product η̂ × ω already has magnitude |ω| (since |η̂| = 1
//   and η̂ ⊥ ω for an ideal vortex), which means the correction is
//   automatically stronger where vorticity is high and vanishes where
//   there is none — exactly the self-limiting behavior we want.
//
//   Two scaling factors complete the formula:
//
//   · ε (epsilon) — a dimensionless user parameter that controls how
//     aggressively we fight numerical dissipation.
//     In the code: cb.vorticityStrength.
//     Too low → mushy, dissipated smoke.
//     Too high → spurious swirling artifacts.
//
//   · h (grid spacing = cellSize) — makes the correction resolution-
//     independent.  On a finer grid the vortices are already better
//     resolved (less dissipation per cell), so the boost per cell
//     should be smaller.  Including h ensures the total injected
//     energy doesn't blow up as you refine the grid.
//
//   Putting it all together, the confinement force (acceleration) is:
//
//     f_conf = ε · h · (η̂ × ω)          (m/s²)       [Fedkiw et al.]
//
// STEP 4: FROM ACCELERATION TO VELOCITY INCREMENT
//
//   Operator splitting integrates f_conf over one timestep Δt
//   (forward Euler, the same scheme used for all sub-steps — see
//   projectCS, "What is operator splitting?"):
//
//     Δv = f_conf · Δt
//        = ε · h · (η̂ × ω) · Δt                      [Fedkiw et al.]
//
//   This is the formula presented at the top of this section and
//   implemented in this kernel.
//
// ── Geometric meaning of each symbol ───────────────────────────────────
// Quick-reference glossary for the formula Δv = ε · h · (η̂ × ω) · dt:
//
//   ω  = curl of velocity (from computeCurlCS).  Points along the
//        *axis* of local rotation (right-hand rule: curl fingers in the
//        spinning direction → thumb = ω).
//
//   |ω| = magnitude of ω — how strongly the fluid spins at this cell.
//        Highest at the center of a vortex, fading outward.
//
//   ∇|ω| = gradient of the scalar field |ω|.  A gradient always points
//        in the direction where the field increases most steeply.
//        Since |ω| peaks at vortex cores and decays outward, ∇|ω|
//        points *radially inward* toward the nearest vortex center.
//
//   η̂  = normalize(∇|ω|).  Same direction as ∇|ω|, unit length.
//        Points radially inward toward the vortex core.
//
//   η̂ × ω = cross product of a radial vector (η̂, toward the core)
//        with an axial vector (ω, along the vortex axis).  The result
//        is perpendicular to both, aligned with the direction the fluid
//        is already spinning. That is, the correction nudges the velocity
//        in the same direction as the existing rotation, amplifying it
//        rather than creating new vortices from scratch.
//        This is a *linear* acceleration (units: m/s², not angular)
//        applied at each cell; its tangential direction means the
//        collective effect across all cells around a vortex is to
//        spin up the existing rotation.
//
//   ε  = user-tunable strength (cb.vorticityStrength).  Too low → mushy;
//        too high → spurious swirling artifacts.
//
//   h  = cellSize = grid-spacing factor so the correction scales correctly
//        with resolution (finer grids need less boost per cell).
//
//   dt = timestep.  The expression ε·h·(η̂ × ω) is an acceleration
//        (force per unit mass); multiplying by dt gives the velocity
//        increment Δv which is added directly to the velocity field.
//        (This is standard operator splitting: each Navier-Stokes term
//        is integrated separately as v_new = v_old + acceleration · dt.)
//
// The result: small eddies and turbulent detail are preserved even on
// coarse grids, giving the smoke a realistic, lively appearance.
//
// ── Index-space shortcut ────────────────────────────────────────────────
// computeCurlCS stores ω_code = ω_physical · cellSize (because it uses
// index-space h = 1 instead of world-space h = cellSize).  This kernel
// exploits that: it reads ω_code directly and skips the explicit cellSize
// multiply, arriving at the same result.  Step-by-step derivation below.
// ============================================================================
#ifdef APPLY_VORTICITY
[numthreads(FLUID_THREADS_3D_X, FLUID_THREADS_3D_Y, FLUID_THREADS_3D_Z)]
void applyVorticityCS(uint3 dtid : SV_DispatchThreadID)
{
    int i = (int)dtid.x;
    int j = (int)dtid.y;
    int k = (int)dtid.z;
    if (i >= cb.gridX || j >= cb.gridY || k >= cb.gridZ)
        return;

    if (texSolid[int3(i, j, k)] == 0)
        return;

    // ── Sample |ω_code| at the 6 axis-aligned neighbors ────────
    // We need ∇|ω| — the *gradient* of the vorticity magnitude.
    //
    // The gradient operator ∇ applied to a scalar field s(x,y,z) produces
    // a vector of its partial derivatives:
    //
    //   ∇s = ( ∂s/∂x, ∂s/∂y, ∂s/∂z )
    //
    // This vector points in the direction of steepest increase of s, and
    // its length is the rate of increase per unit distance. Here s = |ω|,
    // so ∇|ω| points from low-spin regions toward nearby vortex cores.
    //
    // Each partial derivative is again approximated with central
    // differences (same stencil used in computeCurlCS):
    //
    //   ∂|ω|/∂x ≈ ( |ω|(i+1,j,k) − |ω|(i−1,j,k) ) / (2·h)
    //
    // → we need |ω| at the two neighbors along each axis, 6 samples total.
    int il = max(i - 1, 0), ir = min(i + 1, cb.gridX - 1);
    int jb = max(j - 1, 0), jt = min(j + 1, cb.gridY - 1);
    int kb = max(k - 1, 0), kf = min(k + 1, cb.gridZ - 1);

    float magL  = length(texCurl[int3(il, j, k)].xyz);
    float magR  = length(texCurl[int3(ir, j, k)].xyz);
    float magB  = length(texCurl[int3(i, jb, k)].xyz);
    float magT  = length(texCurl[int3(i, jt, k)].xyz);
    float magBk = length(texCurl[int3(i, j, kb)].xyz);
    float magF  = length(texCurl[int3(i, j, kf)].xyz);

    // ── Compute η̂ = normalize(∇|ω_code|) ────────────────────
    // Central differences in index space (h = 1) → the 0.5 is 1/(2·h).
    // The magnitudes above come from ω_code, so ∇|ω_code| is scaled by cellSize
    // relative to the true ∇|ω_physical|. But normalization discards
    // magnitude — only the *direction* survives. Whether we compute
    // ∇|ω_code| or ∇|ω_physical|, the normalized vector η̂ is identical.
    float3 eta = 0.5 * float3(magR - magL, magT - magB, magF - magBk);
    eta = normalize(eta + float3(1e-5, 1e-5, 1e-5));

    // ── Read ω_code at this cell ─────────────────────────────
    // texCurl was written by computeCurlCS using index-space differences,
    // so it contains ω_code = ω_physical · cellSize (see derivation there).
    float3 omega = texCurl[int3(i, j, k)].xyz;

    // ── Compute velocity increment (Δv) ────────────────────────
    // Fedkiw's world-space formula (h = cellSize):
    //   Δv = ε · h · (η̂ × ω_physical) · dt
    //
    // What the code computes:
    //   Δv = ε · (η̂ × omega) · dt
    //      = ε · (η̂ × ω_physical · h) · dt            [substitute ω_code]
    //      = ε · h · (η̂ × ω_physical) · dt            [factor out h]
    //      = Fedkiw's formula                          ✓
    //
    // The h (= cellSize) that was "missing" in computeCurlCS's denominator
    // is now baked into omega, and it acts as the grid-spacing factor h
    // that Fedkiw's formula requires.  No explicit h multiply needed.
    //
    // Note: the variable is named "force" but it is really Δv — the
    // acceleration ε·h·(η̂ × ω) has already been multiplied by dt,
    // producing a velocity increment that is added directly below.
    float3 force = cb.vorticityStrength * float3(
        eta.y * omega.z - eta.z * omega.y,
        eta.z * omega.x - eta.x * omega.z,
        eta.x * omega.y - eta.y * omega.x
    ) * cb.dt;

    // ── Apply Δv to the staggered velocity field ─────────────
    // Each cell writes only to its "owned" faces (lower-index: u[i],
    // v[j], w[k]) to avoid race conditions — adjacent cells share faces,
    // and non-atomic += from two threads would be undefined.
    texU[int3(i, j, k)] += force.x;
    texV[int3(i, j, k)] += force.y;
    texW[int3(i, j, k)] += force.z;
}
#endif

// ============================================================================
// KERNEL: advectTemperatureCS
// Semi-Lagrangian advection of temperature field (no dissipation).
// Temperature is perfectly conserved — only transported by velocity.
//
// This follows the same logic and coordinate conventions as advectSmokeCS:
// backtrace from the current cell center using avgVelAtCenter(), then sample
// the old temperature at that departure point. As with smoke, the backtraced
// position is computed in cell-center coordinates, so we shift by -0.5 before
// calling sampleTemperature().
// ============================================================================
#ifdef ADVECT_TEMPERATURE
[numthreads(FLUID_THREADS_3D_X, FLUID_THREADS_3D_Y, FLUID_THREADS_3D_Z)]
void advectTemperatureCS(uint3 dtid : SV_DispatchThreadID)
{
    int i = (int)dtid.x;
    int j = (int)dtid.y;
    int k = (int)dtid.z;
    if (i >= cb.gridX || j >= cb.gridY || k >= cb.gridZ)
        return;

    if (texSolid[int3(i, j, k)] == 0)
    {
        texTemp_temp[int3(i, j, k)] = 0.0;
        return;
    }

    float h = cb.cellSize;
    float dt = cb.dt;

    // Temperature is at the cell center in grid-index space: (i+0.5, j+0.5, k+0.5)
    float px = (float)i + 0.5;
    float py = (float)j + 0.5;
    float pz = (float)k + 0.5;

    float3 vel = avgVelAtCenter(i, j, k);

    // Backtrace
    float bx = px - vel.x * dt / h;
    float by = py - vel.y * dt / h;
    float bz = pz - vel.z * dt / h;

    bx = clamp(bx, 0.5, (float)cb.gridX - 0.5);
    by = clamp(by, 0.5, (float)cb.gridY - 0.5);
    bz = clamp(bz, 0.5, (float)cb.gridZ - 0.5);

    // Convert from cell-center coordinates back to the texture-index
    // coordinates expected by sampleTemperature(). No dissipation is applied:
    // temperature is only transported, not decayed, in this pass.
    texTemp_temp[int3(i, j, k)] = sampleTemperature(float3(bx - 0.5, by - 0.5, bz - 0.5));
}
#endif


// ============================================================================
// KERNEL: macCormackVelCS
// MacCormack correction for velocity advection (Selle et al. 2007).
//
// This is the velocity counterpart of macCormackSmokeCS:
//   advectVelCS first writes the stable forward Semi-Lagrangian result into
//   texU_temp / texV_temp / texW_temp,
//   then this pass estimates the round-trip error and writes the corrected
//   result back into the original velocity textures texU / texV / texW.
//
// The core MacCormack idea is the same as for smoke:
//   phi_orig = original value at the current sample location
//   phi_fwd  = forward SL result at the current sample location
//   phi_bwd  = reverse-step estimate obtained from the already-advected field
//   phi_mc   = phi_fwd + 0.5 * (phi_orig - phi_bwd)
//
// The important extra detail for velocity is that this field is staggered.
// We do not advect one cell-centered scalar, but three different component
// fields living on different face locations:
//
//   u lives on X faces: (i,     j+0.5, k+0.5)
//   v lives on Y faces: (i+0.5, j,     k+0.5)
//   w lives on Z faces: (i+0.5, j+0.5, k)
//
// That means each component must be corrected separately at its own staggered
// sample position, with its own bounds and coordinate convention.
//
// Why traces use the _temp velocity field
// ---------------------------------------
// Unlike macCormackSmokeCS, here the quantity being corrected is itself the
// velocity field.  To compute the backtrace and forward-trace directions for
// each staggered component, we need the full local velocity vector at that
// staggered sample location.
//
// This pass reconstructs that tracing velocity from texU_temp / texV_temp /
// texW_temp, not from the original texU / texV / texW.
//
// That choice is important:
//   we are applying the reverse-step estimate on the already-advected forward
//   SL velocity field, so the tracing velocity must come from that same field.
//
// In other words:
//   advectVelCS built a forward-SL velocity field in _temp,
//   macCormackVelCS uses that forward-SL field both
//   1. to reconstruct the local velocity used for back/forward tracing, and
//   2. to sample phi_bwd from the forward-traced position.
//
// This is the velocity analog of "go one position forward to apply the
// inverse step to the already-advected field", except that on a staggered
// grid we must first reconstruct the full velocity vector at the face sample.
//
// Neighborhood clamp
// ------------------
// After building phi_mc, we apply the same Selle "revert to SL" limiter used
// for smoke and temperature:
//
//   getMinMaxVel(b, component, mn, mx)
//
// This works like getMinMaxSmoke(), but on the staggered component field
// being corrected:
//   component 0 reads texU_temp with U dimensions
//   component 1 reads texV_temp with V dimensions
//   component 2 reads texW_temp with W dimensions
//
// Important clarification:
// -----------------------
//   this pass does NOT replace applyVorticityCS and does not inject new
//   rotational energy by itself.
//
//   - applyVorticityCS adds a force that reinforces surviving vortices.
//   - macCormackVelCS reduces the numerical diffusion introduced by forward
//     Semi-Lagrangian velocity advection.
//
// In other words:
//   applyVorticityCS fights the loss of visible swirl,
//   macCormackVelCS fights the smoothing error of the advection step.
//
// Both operate on velocity, but they solve different problems and are useful
// together rather than being redundant.
// ============================================================================
#ifdef MACCORMACK_VEL
[numthreads(FLUID_THREADS_3D_X, FLUID_THREADS_3D_Y, FLUID_THREADS_3D_Z)]
void macCormackVelCS(uint3 dtid : SV_DispatchThreadID)
{
    float h = cb.cellSize;
    float dt = cb.dt;
    int i = (int)dtid.x;
    int j = (int)dtid.y;
    int k = (int)dtid.z;

    // --- Correct U component ---
    if (i <= cb.gridX && j < cb.gridY && k < cb.gridZ)
    {
        bool leftFluid  = (i > 0 && isInBounds(i-1,j,k)) ? (texSolid[int3(i-1,j,k)] != 0) : false;
        bool rightFluid = (i < cb.gridX && isInBounds(i,j,k)) ? (texSolid[int3(i,j,k)] != 0) : false;

        if (leftFluid || rightFluid)
        {
            float px = (float)i;
            float py = (float)j + 0.5;
            float pz = (float)k + 0.5;

            // Reconstruct the tracing velocity from the already-advected
            // forward-SL field in _temp, not from the original velocity
            // textures. Since this correction is evaluated on the advected
            // velocity field itself, the trace directions must come from that
            // same field.
            float vel_u = texU_temp[int3(i,j,k)];
            float vel_v = 0.25 * (
                sampleVelocityTemp(float3(px-0.5,py-0.5,pz), 1) +
                sampleVelocityTemp(float3(px+0.5,py-0.5,pz), 1) +
                sampleVelocityTemp(float3(px-0.5,py+0.5,pz), 1) +
                sampleVelocityTemp(float3(px+0.5,py+0.5,pz), 1));
            float vel_w = 0.25 * (
                sampleVelocityTemp(float3(px-0.5,py,pz-0.5), 2) +
                sampleVelocityTemp(float3(px+0.5,py,pz-0.5), 2) +
                sampleVelocityTemp(float3(px-0.5,py,pz+0.5), 2) +
                sampleVelocityTemp(float3(px+0.5,py,pz+0.5), 2));

            // Backtrace (same direction as forward SL)
            float bx = px - vel_u * dt / h;
            float by = py - vel_v * dt / h;
            float bz = pz - vel_w * dt / h;
            bx = clamp(bx, 0.0, (float)cb.gridX);
            by = clamp(by, 0.5, (float)cb.gridY - 0.5);
            bz = clamp(bz, 0.5, (float)cb.gridZ - 0.5);

            // Forward trace (reverse direction)
            float fx = px + vel_u * dt / h;
            float fy = py + vel_v * dt / h;
            float fz = pz + vel_w * dt / h;
            fx = clamp(fx, 0.0, (float)cb.gridX);
            fy = clamp(fy, 0.5, (float)cb.gridY - 0.5);
            fz = clamp(fz, 0.5, (float)cb.gridZ - 0.5);

            float phi_orig = texU[int3(i,j,k)];
            float phi_fwd  = texU_temp[int3(i,j,k)];
            float phi_bwd  = sampleVelocityTemp(float3(fx, fy, fz), 0);

            float phi_mc = phi_fwd + 0.5 * (phi_orig - phi_bwd);

            // If the corrected U value leaves the local min/max range around
            // the departure point (measured on the stable forward-SL U field
            // in texU_temp), discard the correction and fall back to phi_fwd.
            float mn, mx;
            getMinMaxVel(float3(bx, by, bz), 0, mn, mx);
            if (phi_mc < mn || phi_mc > mx)
                phi_mc = phi_fwd;

            texU[int3(i,j,k)] = phi_mc;
        }
    }

    // --- Correct V component ---
    if (i < cb.gridX && j <= cb.gridY && k < cb.gridZ)
    {
        bool belowFluid = (j > 0 && isInBounds(i,j-1,k)) ? (texSolid[int3(i,j-1,k)] != 0) : false;
        bool aboveFluid = (j < cb.gridY && isInBounds(i,j,k)) ? (texSolid[int3(i,j,k)] != 0) : false;

        if (belowFluid || aboveFluid)
        {
            float px = (float)i + 0.5;
            float py = (float)j;
            float pz = (float)k + 0.5;

            float vel_v = texV_temp[int3(i,j,k)];
            float vel_u = 0.25 * (
                sampleVelocityTemp(float3(px-0.5,py-0.5,pz), 0) +
                sampleVelocityTemp(float3(px+0.5,py-0.5,pz), 0) +
                sampleVelocityTemp(float3(px-0.5,py+0.5,pz), 0) +
                sampleVelocityTemp(float3(px+0.5,py+0.5,pz), 0));
            float vel_w = 0.25 * (
                sampleVelocityTemp(float3(px,py-0.5,pz-0.5), 2) +
                sampleVelocityTemp(float3(px,py+0.5,pz-0.5), 2) +
                sampleVelocityTemp(float3(px,py-0.5,pz+0.5), 2) +
                sampleVelocityTemp(float3(px,py+0.5,pz+0.5), 2));

            float bx = px - vel_u * dt / h;
            float by = py - vel_v * dt / h;
            float bz = pz - vel_w * dt / h;
            bx = clamp(bx, 0.5, (float)cb.gridX - 0.5);
            by = clamp(by, 0.0, (float)cb.gridY);
            bz = clamp(bz, 0.5, (float)cb.gridZ - 0.5);

            float fx = px + vel_u * dt / h;
            float fy = py + vel_v * dt / h;
            float fz = pz + vel_w * dt / h;
            fx = clamp(fx, 0.5, (float)cb.gridX - 0.5);
            fy = clamp(fy, 0.0, (float)cb.gridY);
            fz = clamp(fz, 0.5, (float)cb.gridZ - 0.5);

            float phi_orig = texV[int3(i,j,k)];
            float phi_fwd  = texV_temp[int3(i,j,k)];
            float phi_bwd  = sampleVelocityTemp(float3(fx, fy, fz), 1);

            float phi_mc = phi_fwd + 0.5 * (phi_orig - phi_bwd);

            float mn, mx;
            getMinMaxVel(float3(bx, by, bz), 1, mn, mx);
            if (phi_mc < mn || phi_mc > mx)
                phi_mc = phi_fwd;

            texV[int3(i,j,k)] = phi_mc;
        }
    }

    // --- Correct W component ---
    if (i < cb.gridX && j < cb.gridY && k <= cb.gridZ)
    {
        bool backFluid  = (k > 0 && isInBounds(i,j,k-1)) ? (texSolid[int3(i,j,k-1)] != 0) : false;
        bool frontFluid = (k < cb.gridZ && isInBounds(i,j,k)) ? (texSolid[int3(i,j,k)] != 0) : false;

        if (backFluid || frontFluid)
        {
            float px = (float)i + 0.5;
            float py = (float)j + 0.5;
            float pz = (float)k;

            float vel_w = texW_temp[int3(i,j,k)];
            float vel_u = 0.25 * (
                sampleVelocityTemp(float3(px-0.5,py,pz-0.5), 0) +
                sampleVelocityTemp(float3(px+0.5,py,pz-0.5), 0) +
                sampleVelocityTemp(float3(px-0.5,py,pz+0.5), 0) +
                sampleVelocityTemp(float3(px+0.5,py,pz+0.5), 0));
            float vel_v = 0.25 * (
                sampleVelocityTemp(float3(px,py-0.5,pz-0.5), 1) +
                sampleVelocityTemp(float3(px,py+0.5,pz-0.5), 1) +
                sampleVelocityTemp(float3(px,py-0.5,pz+0.5), 1) +
                sampleVelocityTemp(float3(px,py+0.5,pz+0.5), 1));

            float bx = px - vel_u * dt / h;
            float by = py - vel_v * dt / h;
            float bz = pz - vel_w * dt / h;
            bx = clamp(bx, 0.5, (float)cb.gridX - 0.5);
            by = clamp(by, 0.5, (float)cb.gridY - 0.5);
            bz = clamp(bz, 0.0, (float)cb.gridZ);

            float fx = px + vel_u * dt / h;
            float fy = py + vel_v * dt / h;
            float fz = pz + vel_w * dt / h;
            fx = clamp(fx, 0.5, (float)cb.gridX - 0.5);
            fy = clamp(fy, 0.5, (float)cb.gridY - 0.5);
            fz = clamp(fz, 0.0, (float)cb.gridZ);

            float phi_orig = texW[int3(i,j,k)];
            float phi_fwd  = texW_temp[int3(i,j,k)];
            float phi_bwd  = sampleVelocityTemp(float3(fx, fy, fz), 2);

            float phi_mc = phi_fwd + 0.5 * (phi_orig - phi_bwd);

            float mn, mx;
            getMinMaxVel(float3(bx, by, bz), 2, mn, mx);
            if (phi_mc < mn || phi_mc > mx)
                phi_mc = phi_fwd;

            texW[int3(i,j,k)] = phi_mc;
        }
    }
}
#endif

// ============================================================================
// KERNEL: macCormackSmokeCS
// MacCormack correction for smoke density.
//
// Why this pass exists
// --------------------
// The previous pass, advectSmokeCS, already transported smoke with a forward
// Semi-Lagrangian (SL) update and stored the result in texSmoke_temp:
//
//   phi_fwd = A(phi_orig)
//
// where:
//   phi_orig = original smoke at the current cell center
//   A(.)     = one forward SL advection step
//
// That method is robust and stable, but it is numerically diffusive:
// repeated interpolation smooths out sharp smoke features over time.
//
// MacCormack is a second pass that estimates how much detail the forward SL
// step lost, then adds part of that detail back.
//
// Intuition in one sentence:
//   advectSmokeCS gives us a stable answer;
//   macCormackSmokeCS asks "if I try to undo that transport, how far from the
//   original value do I end up?"
//
// Coordinate setup
// ----------------
// Smoke is cell-centered, so the current sample location for cell (i,j,k) is:
//
//   p = (px,py,pz) = (i+0.5, j+0.5, k+0.5)
//
// Using the center velocity vel = avgVelAtCenter(i,j,k), we form:
//
//   b = p - vel * dt / h    // backtraced position
//   f = p + vel * dt / h    // forward-traced position
//
// advectSmokeCS used b to answer:
//
//   "What smoke value from the past arrives at p now?"
//
// and wrote that answer into texSmoke_temp.
//
// The important distinction is:
//   advectSmokeCS samples the past of the value that arrives at p
//   macCormackSmokeCS samples the future of the already-advected field to
//   check whether, when going back, we reconstruct the original well
//
// Why phi_bwd is sampled at the FORWARD position
// ----------------------------------------------
// In the line:
//
//   phi_bwd = sampleSmokeTemp(f - 0.5)
//
// the name phi_bwd refers to the meaning of the value, not to the direction
// of the sample location.
//
// We go one position forward to apply the inverse step to the already-advected
// field.  In other words, texSmoke_temp is the field after the forward SL
// pass, and sampling it at f estimates:
//
//   phi_bwd ≈ A^R(phi_fwd)
//
// where A^R is the reverse / "go back" step.
//
// So:
//   the position f is forward,
//   but the value phi_bwd means "what do I get after trying to go back
//   through the already-advected field?"
//
// The value obtained should ideally reconstruct the original value at the
// current cell center, not simply repeat the sample performed in
// advectSmokeCS.
//
// 1D numerical example
// --------------------
// Suppose the current cell center is:
//
//   p = 1.5
//
// and the local displacement is:
//
//   vel * dt / h = 0.2
//
// Then:
//
//   b = 1.5 - 0.2 = 1.3
//   f = 1.5 + 0.2 = 1.7
//
// Assume the original smoke field around this point is:
//
//   phi_orig[0] = 0
//   phi_orig[1] = 10
//   phi_orig[2] = 0
//
// The forward SL pass backtraces to b = 1.3, samples the ORIGINAL field there,
// and may produce:
//
//   phi_fwd = 8
//
// That 8 is stored in texSmoke_temp at the current cell.
//
// Also suppose that after the full forward SL pass has updated every cell, the
// local values around the same source region have become:
//
//   texSmoke_temp[0] = 0
//   texSmoke_temp[1] = 1
//   texSmoke_temp[2] = 8
//   texSmoke_temp[3] = 5
//   texSmoke_temp[4] = 1
//
// Notice what changed:
//   around the source region near b, the ORIGINAL field had values 10 and 4,
//   while the FORWARD-SL field now has values 8 and 5 in that same area.
//
// So the spatial region is the same, but the sampled field is not.
//
// Now MacCormack goes to the forward position f = 1.7 and samples the
// ALREADY-ADVECTED field texSmoke_temp there.  Suppose that gives:
//
//   phi_bwd = 5.9
//
// If the round-trip were perfect, we would want:
//
//   phi_bwd ≈ phi_orig = 10
//
// But we got 5.9 instead, so the combination
//
//   original -> forward SL -> reverse estimate
//
// clearly lost amplitude.  That loss is the error MacCormack tries to undo.
//
// The correction formula
// ----------------------
// We read three values:
//
//   phi_orig = texSmoke[i,j,k]       // original smoke at p
//   phi_fwd  = texSmoke_temp[i,j,k]  // forward SL result at p
//   phi_bwd  = sample from texSmoke_temp at f
//
// and compute:
//
//   phi_mc = phi_fwd + 0.5 * (phi_orig - phi_bwd)
//
// Interpretation:
//   phi_fwd is the stable but smoothed result
//   phi_orig - phi_bwd estimates how much the forward+reverse round-trip
//   drifted away from the original
//   adding half of that drift puts back some of the lost detail
//
// Continuing the example above:
//
//   phi_orig = 10
//   phi_fwd  = 8
//   phi_bwd  = 5.9
//
// so:
//
//   phi_mc = 8 + 0.5 * (10 - 5.9)
//          = 8 + 2.05
//          = 10.05
//
// which is much closer to the original sharp peak than the raw SL value 8.
//
// Why the min/max clamp is needed
// -------------------------------
// MacCormack reduces diffusion, but near sharp gradients it can overshoot and
// create non-physical oscillations (ringing), such as:
//
//   negative smoke density
//   values larger than all nearby samples
//
// To prevent this, we compute the min and max of the 8 texels surrounding the
// DEPARTURE position b:
//
//   getMinMaxSmoke(b - 0.5, mn, mx)
//
// getMinMaxSmoke() does not interpolate.  It finds the cell cube containing
// the departure point and reads its 8 corner samples from texSmoke_temp.
//
// Important clarification:
//   getMinMaxSmoke() uses the same SOURCE REGION around b that produced the
//   forward-SL sample, but it reads that region from texSmoke_temp, not from
//   the original texSmoke field.
//
// So it is NOT asking:
//   "what were the exact original samples used by advectSmokeCS?"
//
// It IS asking:
//   "in the stable forward-SL field we have just constructed, what is the
//   locally plausible range around the source region of this sample?"
//
// Example:
//   continuing the 1D example above, around b the original field had:
//     { 10, 4 }
//   but the forward-SL field now has:
//     { 8, 5 }
//   so the local forward-SL range is:
//     mn = 5
//     mx = 8
//
//   In 3D the idea is identical, except that we read 8 neighboring samples
//   instead of 2.
//
//   If the 3D neighborhood values were, for example:
//     { 0, 1, 2, 2, 3, 4, 4, 5 }
//   then:
//     mn = 0
//     mx = 5
//
// After that, we check:
//
//   if (phi_mc < mn || phi_mc > mx)
//       phi_mc = phi_fwd;
//
// meaning:
//   if the corrected value leaves the neighborhood's plausible range in the
//   stable forward-SL field, discard the correction and revert to phi_fwd.
//
// In other words, this limiter asks:
//   "is the correction we are proposing too aggressive relative to the stable
//   field we just built?"
//
// If the answer is yes, we keep the safer forward-SL result.
//
// This is the Selle "revert to SL" limiter:
//   keep the sharpening only when it stays locally reasonable;
//   otherwise fall back to the more diffusive but safer phi_fwd.
//
// Why getMinMaxSmoke() reads texSmoke_temp
// ----------------------------------------
// During this pass, texSmoke_temp is the already-computed forward SL field and
// is read-only.  texSmoke is where we write the corrected result.
//
// Reading neighborhood bounds from texSmoke_temp has two practical benefits:
// 1. it avoids races, because we never read partially-corrected texSmoke data;
// 2. it clamps against the stable forward-SL neighborhood that produced the
//    candidate correction in the first place.
//
// Implementation summary
// ----------------------
// 1. Reconstruct the current cell-center position p.
// 2. Recompute the same local center velocity used for advection.
// 3. Form b = p - vel*dt/h and f = p + vel*dt/h.
// 4. Read phi_orig from texSmoke and phi_fwd from texSmoke_temp.
// 5. Sample texSmoke_temp at f to estimate phi_bwd, the reverse-step result.
// 6. Build the corrected value phi_mc.
// 7. Clamp / revert if phi_mc leaves the local neighborhood range.
// 8. Write the final corrected smoke into texSmoke.
// ============================================================================
#ifdef MACCORMACK_SMOKE
[numthreads(FLUID_THREADS_3D_X, FLUID_THREADS_3D_Y, FLUID_THREADS_3D_Z)]
void macCormackSmokeCS(uint3 dtid : SV_DispatchThreadID)
{
    // Current cell index handled by this thread.
    int i = (int)dtid.x;
    int j = (int)dtid.y;
    int k = (int)dtid.z;

	// Bounds check: if outside the grid, do nothing.
	// Note that the forward SL pass already wrote 0 to texSmoke_temp for
	// out-of-bounds cells, so we don't need to write anything here.
    if (i >= cb.gridX || j >= cb.gridY || k >= cb.gridZ)
        return;

    // Solid cells: _temp was already set to 0 by forward SL. Leave it.
    if (texSolid[int3(i,j,k)] == 0)
        return;

	// Cell size and timestep from constant buffer.
    float h = cb.cellSize;
    float dt = cb.dt;

    // Current smoke sample position p in cell-center coordinates.
    float px = (float)i + 0.5;
    float py = (float)j + 0.5;
    float pz = (float)k + 0.5;

    // Reconstruct the same center velocity used to advect this scalar field.
    float3 vel = avgVelAtCenter(i, j, k);

    // Departure point b: where the smoke now at p came from one timestep ago.
    float bx = px - vel.x * dt / h;
    float by = py - vel.y * dt / h;
    float bz = pz - vel.z * dt / h;
    bx = clamp(bx, 0.5, (float)cb.gridX - 0.5);
    by = clamp(by, 0.5, (float)cb.gridY - 0.5);
    bz = clamp(bz, 0.5, (float)cb.gridZ - 0.5);

    // Forward point f: used to apply the reverse step on the already-advected
    // field stored in texSmoke_temp.
    float fx = px + vel.x * dt / h;
    float fy = py + vel.y * dt / h;
    float fz = pz + vel.z * dt / h;
    fx = clamp(fx, 0.5, (float)cb.gridX - 0.5);
    fy = clamp(fy, 0.5, (float)cb.gridY - 0.5);
    fz = clamp(fz, 0.5, (float)cb.gridZ - 0.5);

    // Original value at p before advection.
    float phi_orig = texSmoke[int3(i,j,k)];
    // Stable forward Semi-Lagrangian result produced by advectSmokeCS.
    float phi_fwd  = texSmoke_temp[int3(i,j,k)];
    // Reverse-step estimate: sample the already-advected field at the forward
    // position to see how well the round-trip reconstructs the original.
    float phi_bwd  = sampleSmokeTemp(float3(fx - 0.5, fy - 0.5, fz - 0.5));

    // MacCormack correction:
    // start from the stable SL value and add back half of the round-trip error.
    float phi_mc = phi_fwd + 0.5 * (phi_orig - phi_bwd);

    // If the MacCormack-corrected value leaves the local min/max range
    // around the departure point (measured on the stable forward-SL field
    // in texSmoke_temp), discard the correction and fall back to phi_fwd.
    float mn, mx;
    getMinMaxSmoke(float3(bx - 0.5, by - 0.5, bz - 0.5), mn, mx);
    if (phi_mc < mn || phi_mc > mx)
        phi_mc = phi_fwd;

    // Write the final corrected smoke back to the original field.
    texSmoke[int3(i,j,k)] = phi_mc;
}
#endif

// ============================================================================
// KERNEL: macCormackTempCS
// MacCormack correction for temperature.
//
// This follows the same scheme as macCormackSmokeCS:
//   forward SL advection writes the stable result into texTemp_temp,
//   then this pass estimates the round-trip error and writes the corrected
//   value back to the original temperature field texTemp.
// ============================================================================
#ifdef MACCORMACK_TEMP
[numthreads(FLUID_THREADS_3D_X, FLUID_THREADS_3D_Y, FLUID_THREADS_3D_Z)]
void macCormackTempCS(uint3 dtid : SV_DispatchThreadID)
{
    int i = (int)dtid.x;
    int j = (int)dtid.y;
    int k = (int)dtid.z;
    if (i >= cb.gridX || j >= cb.gridY || k >= cb.gridZ)
        return;

    if (texSolid[int3(i,j,k)] == 0)
        return;

    float h = cb.cellSize;
    float dt = cb.dt;

    float px = (float)i + 0.5;
    float py = (float)j + 0.5;
    float pz = (float)k + 0.5;

    float3 vel = avgVelAtCenter(i, j, k);

    float bx = px - vel.x * dt / h;
    float by = py - vel.y * dt / h;
    float bz = pz - vel.z * dt / h;
    bx = clamp(bx, 0.5, (float)cb.gridX - 0.5);
    by = clamp(by, 0.5, (float)cb.gridY - 0.5);
    bz = clamp(bz, 0.5, (float)cb.gridZ - 0.5);

    float fx = px + vel.x * dt / h;
    float fy = py + vel.y * dt / h;
    float fz = pz + vel.z * dt / h;
    fx = clamp(fx, 0.5, (float)cb.gridX - 0.5);
    fy = clamp(fy, 0.5, (float)cb.gridY - 0.5);
    fz = clamp(fz, 0.5, (float)cb.gridZ - 0.5);

    float phi_orig = texTemp[int3(i,j,k)];
    float phi_fwd  = texTemp_temp[int3(i,j,k)];
    float phi_bwd  = sampleTemperatureTemp(float3(fx - 0.5, fy - 0.5, fz - 0.5));

    float phi_mc = phi_fwd + 0.5 * (phi_orig - phi_bwd);

    float mn, mx;
    getMinMaxTemp(float3(bx - 0.5, by - 0.5, bz - 0.5), mn, mx);
    if (phi_mc < mn || phi_mc > mx)
        phi_mc = phi_fwd;

    texTemp[int3(i,j,k)] = phi_mc;
}
#endif

