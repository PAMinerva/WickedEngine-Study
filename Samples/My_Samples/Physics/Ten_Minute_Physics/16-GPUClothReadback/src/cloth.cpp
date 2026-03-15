#include <algorithm>
#include <cmath>
#include <vector>
#include <wiScene.h>
#include <wiGraphics.h>
#include <wiRenderer.h>
#include <wiHelper.h>
#include "cloth.h"
#include "simulation_utils.h"

using namespace wi::graphics;

ClothMesh::ClothMesh(const SimulationParams& inParams)
    : params(inParams)
{
    numX = params.numX;
    numZ = params.numZ;

    // Force even (required for correct constraint building; see BuildConstraintPasses)
    if (numX % 2 == 1) numX++;
    if (numZ % 2 == 1) numZ++;

    float spacing = params.spacing;
    float clothY = params.clothY;

    //  Particles
    // numX and numZ define the grid dimensions (number of quads in each direction).
    // Since particles are placed at quad vertices, we need (numX+1) × (numZ+1) particles.
    // Example: a 3×2 grid of quads requires 4×3 = 12 particles.
	// ●---●---●---●
	// |   |   |   |
	// ●---●---●---●
	// |   |   |   |
	// ●---●---●---●
    numParticles = (numX + 1) * (numZ + 1);

	//  Positions and inverse masses (CPU side)
    cpuPos.resize(numParticles);
    cpuInvMass.resize(numParticles);

    // Initialize particle positions in a grid.
    // Grid is in XZ plane, centered at origin.
	// The y position is set to the given yOffset.
	// Initialize positions column by column starting from the bottom-left corner and going upwards,
	// then moving to the next column to the right.
	// Indeed, the particle IDs are computed from their grid coordinates (xi, zi) as
	// id = xi * (numZ + 1) + zj
	// which allows to order particles in the Z direction first, then in the X direction.
	// For example:
	// - xi = 0, zi = 0 -> id = 0 (bottom-left corner)
	// - xi = 0, zi = numZ -> id = numZ (top-left corner)
	// - xi = numX, zi = 0 -> id = numX * (numZ + 1) (bottom-right corner)
	// - xi = numX, zi = numZ -> id = (numX + 1) * (numZ + 1) - 1 (top-right corner)
	// etc.
	// This same layout is used for constraint generation and triangle.
    for (int xi = 0; xi <= numX; xi++)
    {
        for (int zi = 0; zi <= numZ; zi++)
        {
            int id = xi * (numZ + 1) + zi;
            float x = (-numX * 0.5f + xi) * spacing;
            float y = clothY;
            float z = (-numZ * 0.5f + zi) * spacing;
            cpuPos[id] = XMFLOAT4(x, y, z, 0.0f);
            cpuInvMass[id] = 1.0f;
        }
    }

    restPos = cpuPos;

    //  Constraints
    BuildConstraintPasses();

    //  Triangles
    numTris = 2 * numX * numZ;  // 2 triangles per quad
    triIds.resize(numTris * 3); // 3 vertex IDs per triangle

    int ti = 0;
    for (int xi = 0; xi < numX; xi++)
    {
		// Each quad is made of two triangles.
		// Each triangle is defined by 3 vertex IDs, selected in counter-clockwise order.
        for (int zi = 0; zi < numZ; zi++)
        {
            int id0 = xi * (numZ + 1) + zi;
            int id1 = (xi + 1) * (numZ + 1) + zi;
            int id2 = (xi + 1) * (numZ + 1) + zi + 1;
            int id3 = xi * (numZ + 1) + zi + 1;
            triIds[ti++] = id0; triIds[ti++] = id1; triIds[ti++] = id2;
            triIds[ti++] = id0; triIds[ti++] = id2; triIds[ti++] = id3;
        }
    }

    //  Normals (CPU side)
    cpuNormals.resize(numParticles, XMFLOAT3(0, 1, 0));
}

void ClothMesh::BuildConstraintPasses()
{
	// ==============================================================================
	// GRAPH COLORING FOR PARALLEL CONSTRAINT SOLVING
	// ==============================================================================
	// This section implements a GEOMETRICALLY PRECOMPUTED graph coloring strategy
	// to enable parallel resolution of distance constraints on the GPU without race conditions.
	//
	// KEY INSIGHT: Rather than computing graph coloring at runtime (which would be expensive),
	// this approach exploits the REGULAR GRID TOPOLOGY of the cloth mesh to derive the coloring
	// analytically. Since the cloth is a 2D rectangular grid with a predictable particle layout,
	// we can precompute the optimal coloring DURING INITIALIZATION using simple mathematical
	// formulas based on particle grid coordinates (xi, zi).
	//
	// BACKGROUND: THE CONSTRAINT DEPENDENCY PROBLEM
	// -----------------------------------------------
	// In cloth simulation, particles are connected by distance constraints (edges).
	// Each constraint enforces that two particles maintain a specific distance.
	// When solving constraints in parallel on a GPU, we face a critical challenge:
	// if two constraints share a particle, they cannot be solved simultaneously
	// because both GPU threads executing the constraint solver would try to
	// read and write the same particle position, causing a race condition
	// and incorrect results.
	//
	// For example, consider particles P0, P1, P2 in a line:
	//   P0 ---- C1 ---- P1 ---- C2 ---- P2
	// If C1 (constraint between P0 and P1) and C2 (constraint between P1 and P2)
	// are solved in parallel, both GPU threads modify P1 simultaneously, leading to
	// a potential race condition and undefined behaviour.
	//
	// GRAPH COLORING SOLUTION
	// -----------------------------------------------
	// Graph coloring assigns "colors" to constraints such that no two constraints
	// of the same color share a particle. This is exactly the constraint graph
	// coloring problem from graph theory.
	//
	// The solution used here is based on a CHECKERBOARD PATTERN applied to the
	// 2D grid of cloth particles. The cloth is a regular grid where particles are
	// indexed as: id = xi * (numZ + 1) + zi (where xi is the X coordinate, zi is the Z coordinate).
	//
	// GEOMETRIC INTUITION
	// -----------------------------------------------
	// Particle layout (the numbers represent particle IDs, computed as id = xi * (numZ + 1) + zi):
	//
	//   4 -- 9 --14 --19 --24
	//   |    |    |    |    |
	//   3 -- 8 --13 --18 --23
	//   |    |    |    |    |
	//   2 -- 7 --12 --17 --22
	//   |    |    |    |    |
	//   1 -- 6 --11 --16 --21
	//   |    |    |    |    |
	//   0 -- 5 --10 --15 --20
	//
	// Particles are indexed as: id = xi * (numZ + 1) + zi
	// Where:
	//   - xi ranges from 0 to numX (left to right, along X-axis)
	//   - zi ranges from 0 to numZ (bottom to top, along Y-axis)
	//
	// For example, in the 5x5 grid above (numX=4, numZ=4):
	//   - Particle at (xi=0, zi=0) -> id = 0 * 5 + 0 = 0
	//   - Particle at (xi=1, zi=0) -> id = 1 * 5 + 0 = 5
	//   - Particle at (xi=0, zi=1) -> id = 0 * 5 + 1 = 1
	//   - Particle at (xi=1, zi=1) -> id = 1 * 5 + 1 = 6
	//   - Particle at (xi=4, zi=4) -> id = 4 * 5 + 4 = 24
	//
	// CONSTRAINT TYPES IN A CLOTH GRID
	// -----------------------------------------------
	// The cloth has four types of distance constraints:
	//
	// 1. STRETCH CONSTRAINTS (Vertical direction, along Z axis)
	//    Connect adjacent particles in the Z direction: (xi, zi) -- (xi, zi+1)
	//    These form vertical lines when viewing the grid from the Z direction.
	//
	// 2. STRETCH CONSTRAINTS (Horizontal direction, along X axis)
	//    Connect adjacent particles in the X direction: (xi, zi) -- (xi+1, zi)
	//    These form horizontal lines when viewing the grid from the Z direction.
	//
	// 3. SHEAR CONSTRAINTS (Diagonals)
	//    Connect diagonal neighbors: (xi, zi) -- (xi+1, zi+1) and (xi+1, zi) -- (xi, zi+1)
	//    These prevent the cloth from collapsing into a line.
	//
	// 4. BENDING CONSTRAINTS (Long-range)
	//    Connect particles at distance 2: (xi, zi) -- (xi, zi+2) and (xi, zi) -- (xi+2, zi)
	//    These prevent excessive bending.
	//
	// THE 5-PASS CHECKERBOARD STRATEGY
	// -----------------------------------------------
	// This implementation splits constraints into 5 passes:
	//
	// PASS 0: Vertical Stretch Constraints
	//   Connects  0-1, 2-3, 5-6, etc. (see the particle layout above for IDs)
	//   Size: (numX + 1) * (numZ / 2)
	//   Independent: YES - No two constraints in this pass share a particle
	//   Why? Particles are connected by edges that are 2 apart in Z direction,
	//   so no particle is involved in two constraints of this pass.
	//
	// PASS 1: Vertical Stretch Constraints
	//   Connects 1-2, 3-4, 6-7, etc.
	//   Size: (numX + 1) * (numZ / 2)
	//   Independent: YES - Complements Pass 0 to cover all vertical stretch edges
	//
	// PASS 2: Horizontal Stretch Constraints
	//   Connects 0-5, 1-6, 2-7, etc.
	//   Size: (numX / 2) * (numZ + 1)
	//   Independent: YES - Same checkerboard logic as Pass 0, but in X direction
	//
	// PASS 3: Horizontal Stretch Constraints
	//   Connects 5-10, 6-11, 7-12, etc.
	//   Size: (numX / 2) * (numZ + 1)
	//   Independent: YES - Complements Pass 2 to cover all horizontal stretch edges
	//
	// PASS 4: Shear + Bending Constraints
	//   Includes all diagonal and long-range connections.
	//   Size: 2*numX*numZ + (numX+1)*(numZ-1) + (numZ+1)*(numX-1)
	//   ((2*numX*numZ) shear; two diagonal for each grid cell), ((numX+1)*(numZ-1) vertical bending, (numZ+1)*(numX-1) horizontal bending)
	//   Imagine to have a 5x5 grid of particles (numX=4, numZ=4), then:
	//   - Shear constraints: 2 * 4 * 4 = 32 (two diagonals for each of the 16 grid cells)
	//   - Vertical bending: (4+1) * (4-1) = 15 (connects particles 2 apart in Z direction; 3 vertical bending constraints per column, and 5 columns)
	//   - Horizontal bending: (4+1) * (4-1) = 15 (connects particles 2 apart in X direction; 3 horizontal bending constraints per row, and 5 rows)
	//   Shear connects 0-6, 1-5, 1-7, 2-6, etc.
	//   Bending connects 0-2, 1-3, 2-4, 0-10, 5-15, 10-20, etc.
	//   Independent: NO - These constraints are interdependent
	//   Why? See the particle layout below:
	//   diagonals 1-5 and 5-11 share particle 5.
	//   vertical lines 0-2 and 2-4 share particle 2.
	//   We could split these into more passes to achieve independence, but that's not the case here.
	//
	// CHECKERBOARD PATTERN VISUALIZATION
	// -----------------------------------------------
	// Pass 0 (Vertical Stretch):
	// Constraints marked with '|' connect particles with ID=0<->1, ID=2<->3, ID=4<->5, etc.
	//
	//   4    9   14   19   24     (these particles are NOT in Pass 0)
	//
	//   3    8   13   18   23
	//   |    |    |    |    |
	//   2    7   12   17   22
	//
	//   1    6   11   16   21
	//   |    |    |    |    |
	//   0    5   10   15   20
	//
	// No particle appears in two constraints of Pass 0 because each constraint
	// connects particles that are 2 apart in Z direction.
	//
	// Pass 1 (Vertical Stretch):
	// Constraints marked with '|' connect particles with ID=1<->2, ID=3<->4, etc.
	//
	//   4    9   14   19   24
	//   |    |    |    |    |
	//   3    8   13   18   23
	//
	//   2    7   12   17   22
	//   |    |    |    |    |
	//   1    6   11   16   21
	//
	//   0    5   10   15   20     (these particles are NOT in Pass 1)
	//
	// Pass 1 covers the constraints that Pass 0 misses, with the same guarantee
	// of no shared particles within a single pass.
	//
	// Pass 2 (Horizontal Stretch, Even X):
	// Constraints marked with '--' connect particles with ID=0<->5, ID=1<->6, ID=2<->7, etc.
	//
	//   4 -- 9   14 --19   24
	//
	//   3 -- 8   13 --18   23
	//
	//   2 -- 7   12 --17   22
	//
	//   1 -- 6   11 --16   21
	//
	//   0 -- 5   10 --15   20
	//
	// Pass 3 (Horizontal Stretch, Odd X):
	// Constraints marked with '--' connect particles with ID=5<->10, ID=6<->11, etc.
	//   4    9 --14   19 --24
	//
	//   3    8 --13   18 --23
	//
	//   2    7 --12   17 --22
	//
	//   1    6 --11   16 --21
	//
	//   0    5 --10   15 --20
	//
	// EXECUTION STRATEGY IN THE SIMULATOR
	// -----------------------------------------------
	// During cloth simulation (in the SimulateGPU() method):
	//
	// 1. For each each of the indipendent passes 0-3 (independent constraints):
	//    - Use Gauss-Seidel method.
	//    - Execute a compute shader in parallel, where each GPU thread processes one constraint of that pass.
	//    - The compute shader directly modifies particle positions.
	//    - No synchronization needed since constraints in the same pass do not share particles.
	//    - Each constraint adjusts positions based on violated distance constraints.
	//    - Results are immediately available for the next pass.
	//
	// 2. For Pass 4 (dependent constraints):
	//    - Use Jacobi method.
	//    - Execute a compute shader in parallel, where each GPU thread processes one constraint of Pass 4.
	//    - The compute shader CAN'T directly modify particle positions due to shared particles.
	//    - Instead, for each constraint we compute its position corrections and stores\accumulate them in a separate buffer.
	//    - However, since constraints are interdependent, we cannot write in this second buffer without synchronization.
	//    - We need to use an atomic operation (e.g., atomic add).
	//    - After all constraint corrections are accumulated, apply them to positions with a scale factor
	//
	// PERFORMANCE IMPLICATIONS
	// -----------------------------------------------
	// - Passes 0-3: Fully parallel, no synchronization required.
	// - Pass 4: Requires little synchronization (atomic additions).
	//
	// ==============================================================================

	// Initialize pass sizes based on the number of constraints in each pass.
	// The last pass (shear + bending) is not independent because constraints share particles.
    passSizes[0] = (numX + 1) * (numZ / 2);           // Pass 0: vertical stretch
    passSizes[1] = (numX + 1) * (numZ / 2);           // Pass 1: vertical stretch
    passSizes[2] = (numX / 2) * (numZ + 1);           // Pass 2: horizontal stretch
    passSizes[3] = (numX / 2) * (numZ + 1);           // Pass 3: horizontal stretch
    passSizes[4] = 2 * numX * numZ                    // Pass 4: shear (two diagonals per quad) +
                + (numX + 1) * (numZ - 1)             //         bending Y (connects particles 2 apart in Z direction) +
                + (numZ + 1) * (numX - 1);            //         bending X (connects particles 2 apart in X direction)

	// Total number of constraints is the sum of all passes.
    numDistConstraints = 0;
    for (int p = 0; p < NUM_PASSES; p++)
        numDistConstraints += passSizes[p];

	// Each distance constraint is defined by the IDs of the two particles it connects.
	// The two IDs of each constraint are contiguously stored in the constIds vector.
	// Format: [particle_id_0_of_constraint_0, particle_id_1_of_constraint_0,
    //          particle_id_0_of_constraint_1, particle_id_1_of_constraint_1, ...]
    constIds.resize(2 * numDistConstraints);

    int ci = 0;

	// ==============================================================================
	// STRETCH CONSTRAINTS (Vertical Direction, along Z-axis): Pass 0 and Pass 1
	// ==============================================================================
	// Two passes are needed to cover all rows of vertical stretch constraints without shared particles.
    for (int p = 0; p < 2; p++)
    {
        for (int xi = 0; xi <= numX; xi++)
        {
            for (int zi = 0; zi < numZ / 2; zi++)
            {
                constIds[2 * ci]     = xi * (numZ + 1) + 2 * zi + p;
                constIds[2 * ci + 1] = xi * (numZ + 1) + 2 * zi + p + 1;
                ci++;
            }
        }
    }

	// ==============================================================================
	// STRETCH CONSTRAINTS (Horizontal Direction, along X-axis): Pass 2 and Pass 3
	// =============================================================================
	// Two passes are needed to cover all columns of horizontal stretch constraints without shared particles.
    for (int p = 0; p < 2; p++)
    {
        for (int xi = 0; xi < numX / 2; xi++)
        {
            for (int zi = 0; zi <= numZ; zi++)
            {
                constIds[2 * ci]     = (2 * xi + p) * (numZ + 1) + zi;
                constIds[2 * ci + 1] = (2 * xi + p + 1) * (numZ + 1) + zi;
                ci++;
            }
        }
    }

	// ==============================================================================
	// SHEAR CONSTRAINTS (Diagonals): Pass 4
	// =============================================================================
    for (int xi = 0; xi < numX; xi++)
    {
        for (int zi = 0; zi < numZ; zi++)
        {
            constIds[2 * ci]     = xi * (numZ + 1) + zi;
            constIds[2 * ci + 1] = (xi + 1) * (numZ + 1) + zi + 1;
            ci++;
            constIds[2 * ci]     = (xi + 1) * (numZ + 1) + zi;
            constIds[2 * ci + 1] = xi * (numZ + 1) + zi + 1;
            ci++;
        }
    }

	// ==============================================================================
	// BENDING CONSTRAINTS (Long-range): Pass 4 continued
	// =============================================================================
	// Vertical bending constraints (connecting particles 2 units apart in the Z direction)
    for (int xi = 0; xi <= numX; xi++)
    {
        for (int zi = 0; zi < numZ - 1; zi++)
        {
            constIds[2 * ci]     = xi * (numZ + 1) + zi;
            constIds[2 * ci + 1] = xi * (numZ + 1) + zi + 2;
            ci++;
        }
    }

	// Horizontal bending constraints (connecting particles 2 units apart in the X direction)
    for (int xi = 0; xi < numX - 1; xi++)
    {
        for (int zi = 0; zi <= numZ; zi++)
        {
            constIds[2 * ci]     = xi * (numZ + 1) + zi;
            constIds[2 * ci + 1] = (xi + 2) * (numZ + 1) + zi;
            ci++;
        }
    }
}

void ClothMesh::InitGPUBuffers()
{
    CreateGPUBuffers();
    LoadShaders();

    // Rest lengths will be computed on the main thread's first simulation frame (see SimulateGPU) to
    // avoid GPU synchronization at startup (SubmitCommandLists/WaitForGPU) that could causes a VK_TIMEOUT.
    gpuBuffersReady = true;
}

void ClothMesh::CreateGPUBuffers()
{
    GraphicsDevice* device = wi::graphics::GetDevice();

    auto makeRWStructured = [&](GPUBuffer& buf, uint32_t stride, uint32_t count,
                                const void* initData, const char* name)
    {
        GPUBufferDesc desc;
        desc.usage = Usage::DEFAULT;
        desc.bind_flags = BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;
        desc.misc_flags = ResourceMiscFlag::BUFFER_STRUCTURED;
        desc.stride = stride;
        desc.size = stride * count;
        device->CreateBuffer(&desc, initData, &buf);
        device->SetName(&buf, name);
    };

    auto makeReadback = [&](GPUBuffer& buf, uint32_t totalSize, const char* name)
    {
        GPUBufferDesc desc;
        desc.usage = Usage::READBACK;
        desc.bind_flags = BindFlag::NONE;
        desc.misc_flags = ResourceMiscFlag::NONE;
        desc.size = totalSize;
        device->CreateBuffer(&desc, nullptr, &buf);
        device->SetName(&buf, name);
    };

    // Zero-initialized velocity
    std::vector<XMFLOAT4> zeroVel(numParticles, XMFLOAT4(0, 0, 0, 0));

    // Zero-initialized uint buffers (for correctionsBuffer, 3 per particle)
    std::vector<uint32_t> zeroUints(3 * numParticles, 0);

    // Particle buffers
    makeRWStructured(posBuffer, sizeof(XMFLOAT4), numParticles, cpuPos.data(), "cloth::pos");
    makeRWStructured(prevPosBuffer, sizeof(XMFLOAT4), numParticles, cpuPos.data(), "cloth::prevPos");
    makeRWStructured(velBuffer, sizeof(XMFLOAT4), numParticles, zeroVel.data(), "cloth::vel");
    makeRWStructured(invMassBuffer, sizeof(float), numParticles, cpuInvMass.data(), "cloth::invMass");

    // Constraint buffers
    makeRWStructured(constIdsBuffer, sizeof(uint32_t), numDistConstraints * 2, constIds.data(), "cloth::constIds");
    makeRWStructured(restLengthsBuffer, sizeof(float), numDistConstraints, nullptr, "cloth::restLengths");

    // Jacobi corrections buffer
	// (uint, 3 per particle to accumulate x,y,z corrections separately with atomic operations,
	// then convert to float and apply to positions in a second pass;
	// we can't directly accumulate float corrections with atomics in HLSL,
	// but we can do it by encoding the float bits into uints and using atomic operations on the uints,
	// which is legit in HLSL. Then in the second pass we decode the accumulated uints back into floats to apply the corrections).
    makeRWStructured(correctionsBuffer, sizeof(uint32_t), 3 * numParticles, zeroUints.data(), "cloth::corrections");

    // Normals buffer
    makeRWStructured(normalsBuffer, sizeof(XMFLOAT4), numParticles, nullptr, "cloth::normals");

    // Triangle buffer
    makeRWStructured(triIdsBuffer, sizeof(uint32_t), numTris * 3, triIds.data(), "cloth::triIds");

    // Normal accumulation buffer
	// (uint, 3 per particle to accumulate x,y,z normal components separately with atomic operations,
	// then convert to float and normalize in a second pass)
	makeRWStructured(normAccumBuffer, sizeof(uint32_t), 3 * numParticles, zeroUints.data(), "cloth::normAccum");

    // Raycast output
    makeRWStructured(triDistBuffer, sizeof(float), numTris, nullptr, "cloth::triDist");

    // Readback buffers
    makeReadback(posReadbackBuffer, sizeof(XMFLOAT4) * numParticles, "cloth::posReadback");
    makeReadback(normalsReadbackBuffer, sizeof(XMFLOAT4) * numParticles, "cloth::normalsReadback");
    makeReadback(triDistReadbackBuffer, sizeof(float) * numTris, "cloth::triDistReadback");

    // Dispatch sizes
    dispatchParticles = (numParticles + CLOTH_THREAD_GROUP_SIZE - 1) / CLOTH_THREAD_GROUP_SIZE;
    dispatchTris = (numTris + CLOTH_THREAD_GROUP_SIZE - 1) / CLOTH_THREAD_GROUP_SIZE;
}

void ClothMesh::LoadShaders()
{
	// Save original binary and source shader paths to restore after loading cloth shaders
    const std::string originalBinPath = wi::renderer::GetShaderPath();
    const std::string originalSrcPath = wi::renderer::GetShaderSourcePath();

    // .cso outputs go to the shaders folder in the binary directory
    std::string shaderBinPath = wi::helper::GetCurrentPath() + "/shaders/";
    wi::renderer::SetShaderPath(shaderBinPath);

    // .hlsl sources are read directly from the source tree (no copy needed).
    // SHADER_SOURCE_DIR is set by CMake in CMakeLists.txt via target_compile_definitions.
#ifdef SHADER_SOURCE_DIR
    wi::renderer::SetShaderSourcePath(SHADER_SOURCE_DIR);
#else
    wi::renderer::SetShaderSourcePath(shaderBinPath);
#endif

    // All kernels live in cloth_simulation.hlsl. Each LoadShader call uses:
    //   - same base filename "cloth_simulation.cso" -> source = cloth_simulation.hlsl
    //   - a unique permutation define -> unique .cso name (e.g. cloth_simulation_COMPUTE_REST_LENGTHS.cso)
    //   - a named entry point matching the function name in the shader (e.g. cloth_computeRestLengthsCS)
    using SM = ShaderStage;
    wi::renderer::LoadShader(SM::CS, computeRestLengthsCS, "cloth_simulation.cso", ShaderModel::SM_6_0, {"COMPUTE_REST_LENGTHS"}, "cloth_computeRestLengthsCS");
    wi::renderer::LoadShader(SM::CS, integrateCS,          "cloth_simulation.cso", ShaderModel::SM_6_0, {"INTEGRATE"},            "cloth_integrateCS");
    wi::renderer::LoadShader(SM::CS, solveConstraintsCS,   "cloth_simulation.cso", ShaderModel::SM_6_0, {"SOLVE_CONSTRAINTS"},    "cloth_solveConstraintsCS");
    wi::renderer::LoadShader(SM::CS, addCorrectionsCS,     "cloth_simulation.cso", ShaderModel::SM_6_0, {"ADD_CORRECTIONS"},      "cloth_addCorrectionsCS");
    wi::renderer::LoadShader(SM::CS, updateVelCS,          "cloth_simulation.cso", ShaderModel::SM_6_0, {"UPDATE_VEL"},           "cloth_updateVelCS");
    wi::renderer::LoadShader(SM::CS, clearNormalsCS,       "cloth_simulation.cso", ShaderModel::SM_6_0, {"CLEAR_NORMALS"},        "cloth_clearNormalsCS");
    wi::renderer::LoadShader(SM::CS, addNormalsCS,         "cloth_simulation.cso", ShaderModel::SM_6_0, {"ADD_NORMALS"},          "cloth_addNormalsCS");
    wi::renderer::LoadShader(SM::CS, normalizeNormalsCS,   "cloth_simulation.cso", ShaderModel::SM_6_0, {"NORMALIZE_NORMALS"},    "cloth_normalizeNormalsCS");
    wi::renderer::LoadShader(SM::CS, raycastTriangleCS,    "cloth_simulation.cso", ShaderModel::SM_6_0, {"RAYCAST_TRIANGLE"},     "cloth_raycastTriangleCS");

	// Old LoadShader calls without permutation defines (left here for reference, but not used since we switched to a single .cso with permutations)
    // wi::renderer::LoadShader(ShaderStage::CS, computeRestLengthsCS, "cloth_computeRestLengthsCS.cso");
    // wi::renderer::LoadShader(ShaderStage::CS, integrateCS, "cloth_integrateCS.cso");
    // wi::renderer::LoadShader(ShaderStage::CS, solveConstraintsCS, "cloth_solveConstraintsCS.cso");
    // wi::renderer::LoadShader(ShaderStage::CS, addCorrectionsCS, "cloth_addCorrectionsCS.cso");
    // wi::renderer::LoadShader(ShaderStage::CS, updateVelCS, "cloth_updateVelCS.cso");
    // wi::renderer::LoadShader(ShaderStage::CS, clearNormalsCS, "cloth_clearNormalsCS.cso");
    // wi::renderer::LoadShader(ShaderStage::CS, addNormalsCS, "cloth_addNormalsCS.cso");
    // wi::renderer::LoadShader(ShaderStage::CS, normalizeNormalsCS, "cloth_normalizeNormalsCS.cso");
    // wi::renderer::LoadShader(ShaderStage::CS, raycastTriangleCS, "cloth_raycastTriangleCS.cso");

	// Restore original shader paths after loading cloth shaders
    wi::renderer::SetShaderPath(originalBinPath);
    wi::renderer::SetShaderSourcePath(originalSrcPath);
}

void ClothMesh::ComputeRestLengthsGPU(wi::graphics::CommandList cmd)
{
    GraphicsDevice* device = wi::graphics::GetDevice();

	// Total number of distance constraints across all passes (rest lengths are computed for all constraints in one compute shader dispatch)
    ClothSimConstants cb = {};
    cb.numConstraintsInPass = numDistConstraints;

	// Bind the compute shader for rest length computation and set up the constant buffer with the total number of constraints.
	// The compute shader will use this information to determine if the current thread ID corresponds to a valid constraint.
	// See the compute shader "cloth_computeRestLengthsCS" in "cloth_simulation.hlsl"
	// for details on how the rest lengths are computed.
    device->BindComputeShader(&computeRestLengthsCS, cmd);
    device->BindDynamicConstantBuffer(cb, CLOTH_CB_SLOT, cmd);

	// Bind the necessary buffers as UAVs for the compute shader:
	// - posBuffer: contains the current positions of particles, used to compute the rest lengths based on initial particle positions.
	// - constIdsBuffer: contains the pairs of particle IDs for each constraint, used to identify which particles are connected by each constraint.
	// - restLengthsBuffer: the output buffer where the computed rest lengths for each constraint will be stored by the compute shader.
	// Note: We bind posBuffer as a UAV because the compute shader might read from it to calculate rest lengths based on initial positions,
	// but it won't write to it. The constIdsBuffer is also read-only for the compute shader, while the restLengthsBuffer is written by the
	// compute shader to store the results.
    const GPUResource* uav_pos[] = { &posBuffer };
    device->BindUAVs(uav_pos, CLOTH_SLOT_POS, 1, cmd);
    const GPUResource* uav_const[] = { &constIdsBuffer, &restLengthsBuffer };
    device->BindUAVs(uav_const, CLOTH_SLOT_CONST_IDS, 2, cmd);

	// Dispatch the compute shader with enough thread groups to cover all distance constraints.
    uint32_t dispatchCount = (numDistConstraints + CLOTH_THREAD_GROUP_SIZE - 1) / CLOTH_THREAD_GROUP_SIZE;
    device->Dispatch(dispatchCount, 1, 1, cmd);

	// Insert a memory barrier to ensure that the rest lengths are fully computed and visible in the
	// restLengthsBuffer before any subsequent compute shader dispatches that might read from it.
    GPUBarrier barrier = GPUBarrier::Memory();
    device->Barrier(&barrier, 1, cmd);
}

void ClothMesh::SimulateGPU(float dt, wi::graphics::CommandList cmd, int solveType)
{
    if (!gpuBuffersReady)
        return;

    // Compute rest lengths on first simulation frame (deferred from init to avoid VK_TIMEOUT)
    if (!restLengthsComputed)
    {
        ComputeRestLengthsGPU(cmd);
        restLengthsComputed = true;
    }

    GraphicsDevice* device = wi::graphics::GetDevice();
    const float substepDt = gPhysicsScene.dt / static_cast<float>(numSubSteps);
    // const float substepDt = dt / static_cast<float>(numSubSteps);

	// Memory barrier that will be used between passes to ensure correct ordering and
	// visibility of UAV writes and reads across different compute shader dispatches.
    GPUBarrier barrier = GPUBarrier::Memory();

    // Helper lambdas for binding UAV groups
    auto bindIntegrateUAVs = [&]() {
        const GPUResource* uavs[] = { &posBuffer, &prevPosBuffer, &velBuffer, &invMassBuffer };
        device->BindUAVs(uavs, CLOTH_SLOT_POS, 4, cmd);
    };

    auto bindSolveConstraintUAVs = [&]() {
        const GPUResource* uav_pos[] = { &posBuffer };
        device->BindUAVs(uav_pos, CLOTH_SLOT_POS, 1, cmd);
        const GPUResource* uav_mass[] = { &invMassBuffer, &constIdsBuffer, &restLengthsBuffer, &correctionsBuffer };
        device->BindUAVs(uav_mass, CLOTH_SLOT_INV_MASS, 4, cmd);
    };

    auto bindAddCorrUAVs = [&]() {
        const GPUResource* uav_pos[] = { &posBuffer };
        device->BindUAVs(uav_pos, CLOTH_SLOT_POS, 1, cmd);
        const GPUResource* uav_corr[] = { &correctionsBuffer };
        device->BindUAVs(uav_corr, CLOTH_SLOT_CORRECTIONS, 1, cmd);
    };

    auto bindUpdateVelUAVs = [&]() {
        const GPUResource* uavs[] = { &posBuffer, &prevPosBuffer, &velBuffer };
        device->BindUAVs(uavs, CLOTH_SLOT_POS, 3, cmd);
    };

	// ---------------------------------------------------------------------
	// THE SIMULATION LOOP (PBD-like workflow)
	// ---------------------------------------------------------------------
	// The simulation step in simulate() follows a common Position-Based Dynamics pattern:
	//
	//   (1) integrate():       "predict" motion under external forces
	//       - update vel using gravity
	//       - update pos using vel
	//       - handle collisions by directly projecting pos out of penetrations
	//
	//   (2) solve constraints: enforce cloth distance constraints by directly correcting positions
	//       - this is why constraints operate on pos (not on forces)
	//       - the solver may run multiple passes / substeps to converge
	//
	//   (3) updateVel():       recompute velocities from the final corrected positions
	//       vel = (pos - prevPos) / dt
	//
	// In other words: integrate() produces a tentative motion, constraints/collisions fix it,
	// and then velocity is derived from what actually happened.
	// ---------------------------------------------------------------------
	// WHY SUBSTEPPING?
	// ---------------------------------------------------------------------
	// This simulator uses a Position-Based Dynamics (PBD-like) loop:
	//   1) integrate()        -> predict positions under forces (gravity) for a time step dt
	//   2) constraint solve   -> project positions to satisfy distance constraints
	//   3) updateVel()        -> recompute velocity from final corrected positions
	//
	// IMPORTANT: the constraint solve is ITERATIVE / approximate.
	// - We do not solve a global system exactly in one shot.
	// - We apply many local projections (constraints) and stop after a limited amount of work.
	//
	// If we used only ONE step per frame with dt = timeStep:
	// - integration would move particles farther in one go (Δx ~ v*dt)
	// - constraint violations would be larger (springs/distances are more "broken")
	// - collision penetrations can become deeper (bigger push-outs)
	// - the required corrections become large and can lead to:
	//     * overshoot (corrections "over-compensate" (exceeds the constaint) and create new violations with opposite sign)
	//     * jitter / oscillations (constraints keep fighting each other)
	//     * visible stretching (the solver cannot propagate corrections through the cloth fast enough)
	//
	// Substepping splits one frame into multiple smaller steps:
	//   dt = timeStep / numSubsteps
	// and repeats "integrate + solve" numSubsteps times.
	//
	// Benefits:
	// - Smaller per-step motion -> smaller constraint violations -> smaller corrections
	// - Fewer overshoot/oscillations because we correct gently and frequently
	// - Better propagation through a network of constraints:
	//   In a cloth, each correction changes neighboring constraints. With large dt, a big local
	//   correction can disrupt neighbors, and with limited iterations you may not reach a
	//   globally consistent state in one frame. With substeps, corrections spread gradually.
	// - More stable collision handling:
	//   Deep penetrations in one large step produce large "push-out" projections -> noisy motion.
	//   With substeps, penetrations are shallower and projections are smoother.
	// - Less tunneling through colliders:
	//   With large dt, fast particles can tunnel through thin colliders (e.g., the ground plane) because
	//   they move from one side to the other without ever being detected as penetrating.
	//   Substepping reduces this risk by checking collisions more frequently with smaller motions.
	//
	// Trade-off:
	// - More substeps means more kernel launches and more compute per rendered frame.
	//
	// Practical intuition:
	// - numSubsteps acts like a "stability / stiffness budget" for PBD.
	// - If you reduce numSubsteps you may see more stretch and jitter; increasing it typically
	//   makes the cloth appear stiffer and more stable (at higher cost).
	// ---------------------------------------------------------------------
    for (int step = 0; step < numSubSteps; step++)
    {
        // ===== 1. Integrate =====
		// See cloth_integrateCS in cloth_simulation.hlsl for details on the integration step.
        {
            ClothSimConstants cb = {};
            cb.dt = substepDt;
            cb.gravX = gravity.x; cb.gravY = gravity.y; cb.gravZ = gravity.z;
            cb.sphereCX = params.sphereCenter[0];
            cb.sphereCY = params.sphereCenter[1];
            cb.sphereCZ = params.sphereCenter[2];
            cb.sphereR = params.sphereRadius;
            cb.numParticles = numParticles;
            cb.dragParticleNr = dragParticleNr;
            cb.dragPosX = currentDragPos.x;
            cb.dragPosY = currentDragPos.y;
            cb.dragPosZ = currentDragPos.z;

			// Bind the compute shader responsible for integrating particle positions and velocities based on forces, collisions, and drag
			// (see cloth_integrateCS in cloth_simulation.hlsl).
			// Bind the constant buffer with integration parameters (time step, gravity, sphere collision parameters, drag parameters).
			// Bind the necessary buffers for integration (current and previous positions, velocities, inverse masses).
			// Dispatch the compute shader with enough thread groups to cover all particles.
			// Insert a memory barrier to ensure that all position updates from integration are visible before constraint solving.
            device->BindComputeShader(&integrateCS, cmd);
            device->BindDynamicConstantBuffer(cb, CLOTH_CB_SLOT, cmd);
            bindIntegrateUAVs();
            device->Dispatch(dispatchParticles, 1, 1, cmd);
            device->Barrier(&barrier, 1, cmd);
        }

        // ===== 2. Solve constraints =====
		// See cloth_solveConstraintsCS in cloth_simulation.hlsl and cloth_addCorrectionsCS for details on
		// the hybrid constraint solving approach (coloring + Jacobi)
        if (solveType == 0)
        {
            // Coloring hybrid
			// Solve distance constraints in multiple passes.
			// The first 4 passes can be solved in parallel with direct position updates thanks to a careful
			// construction of independent constraint sets using an approach that looks like the graph coloring one.
			// The 5th pass (shear + bending) is solved using the Jacobi method, which requires accumulating
			// corrections and applying them in a separate step at the end of the Jacobi pass.
            int first = 0;
            for (int pNr = 0; pNr < NUM_PASSES; pNr++)
            {
                int n = passSizes[pNr];
                if (n == 0) { first += n; continue; }

                ClothSimConstants cb = {};
                cb.dt = substepDt;
                // jacobiScale=0 signals coloring (direct write), >0 signals Jacobi (atomic accumulate)
                cb.jacobiScale = passIndependent[pNr] ? 0.0f : jacobiScale;
                cb.numParticles = numParticles;
                cb.firstConstraint = first;
                cb.numConstraintsInPass = n;
                cb.dragParticleNr = dragParticleNr;

                device->BindComputeShader(&solveConstraintsCS, cmd);
                device->BindDynamicConstantBuffer(cb, CLOTH_CB_SLOT, cmd);
                bindSolveConstraintUAVs();

                uint32_t dc = (n + CLOTH_THREAD_GROUP_SIZE - 1) / CLOTH_THREAD_GROUP_SIZE;
                device->Dispatch(dc, 1, 1, cmd);
                device->Barrier(&barrier, 1, cmd);

                if (!passIndependent[pNr])
                {
                    // Jacobi pass: apply accumulated corrections
                    device->BindComputeShader(&addCorrectionsCS, cmd);
                    device->BindDynamicConstantBuffer(cb, CLOTH_CB_SLOT, cmd);
                    bindAddCorrUAVs();
                    device->Dispatch(dispatchParticles, 1, 1, cmd);
                    device->Barrier(&barrier, 1, cmd);
                }

                first += n;
            }
        }
        else
        {
            // Full Jacobi
            // We solve all constraints using the Jacobi method in a single pass.
            // This is less efficient than solving independent passes in parallel, but it demonstrates the
			// flexibility of the Jacobi approach to handle all constraints uniformly, including dependent ones.
            // In this mode, we dispatch one compute shader to compute corrections for all constraints, and then
			// apply those corrections to positions in a separate step at the end.
            // This approach is more computationally expensive due to the iterative nature and the fact that all
			// constraints are treated as dependent, but it can be useful for scenarios where all constraints
			// have complex dependencies.
            ClothSimConstants cb = {};
            cb.dt = substepDt;
            cb.jacobiScale = jacobiScale;
            cb.numParticles = numParticles;
            cb.firstConstraint = 0;
            cb.numConstraintsInPass = numDistConstraints;
            cb.dragParticleNr = dragParticleNr;

            device->BindComputeShader(&solveConstraintsCS, cmd);
            device->BindDynamicConstantBuffer(cb, CLOTH_CB_SLOT, cmd);
            bindSolveConstraintUAVs();

            uint32_t dc = (numDistConstraints + CLOTH_THREAD_GROUP_SIZE - 1) / CLOTH_THREAD_GROUP_SIZE;
            device->Dispatch(dc, 1, 1, cmd);
            device->Barrier(&barrier, 1, cmd);

            device->BindComputeShader(&addCorrectionsCS, cmd);
            device->BindDynamicConstantBuffer(cb, CLOTH_CB_SLOT, cmd);
            bindAddCorrUAVs();
            device->Dispatch(dispatchParticles, 1, 1, cmd);
            device->Barrier(&barrier, 1, cmd);
        }

        // ===== 3. Update velocities =====
        {
            ClothSimConstants cb = {};
            cb.dt = substepDt;
            cb.numParticles = numParticles;
            cb.dragParticleNr = dragParticleNr;

            device->BindComputeShader(&updateVelCS, cmd);
            device->BindDynamicConstantBuffer(cb, CLOTH_CB_SLOT, cmd);
            bindUpdateVelUAVs();
            device->Dispatch(dispatchParticles, 1, 1, cmd);
            device->Barrier(&barrier, 1, cmd);
        }
    }
}

void ClothMesh::UpdateMeshNormalsGPU(wi::graphics::CommandList cmd)
{
    if (!gpuBuffersReady)
        return;

    GraphicsDevice* device = wi::graphics::GetDevice();
    GPUBarrier barrier = GPUBarrier::Memory();

    // 1. Clear normAccum
    {
        ClothSimConstants cb = {};
        cb.numParticles = numParticles;

        device->BindComputeShader(&clearNormalsCS, cmd);
        device->BindDynamicConstantBuffer(cb, CLOTH_CB_SLOT, cmd);
        const GPUResource* uav_norm[] = { &normAccumBuffer };
        device->BindUAVs(uav_norm, CLOTH_SLOT_NORM_ACCUM, 1, cmd);
        device->Dispatch(dispatchParticles, 1, 1, cmd);
        device->Barrier(&barrier, 1, cmd);
    }

    // 2. Add normals (one thread per triangle)
    {
        ClothSimConstants cb = {};
        cb.numParticles = numParticles;
        cb.numConstraintsInPass = numTris;  // repurposed for numTris

        device->BindComputeShader(&addNormalsCS, cmd);
        device->BindDynamicConstantBuffer(cb, CLOTH_CB_SLOT, cmd);
        const GPUResource* uav_pos[] = { &posBuffer };
        device->BindUAVs(uav_pos, CLOTH_SLOT_POS, 1, cmd);
        const GPUResource* uav_tri[] = { &triIdsBuffer, &normAccumBuffer };
        device->BindUAVs(uav_tri, CLOTH_SLOT_TRI_IDS, 2, cmd);
        device->Dispatch(dispatchTris, 1, 1, cmd);
        device->Barrier(&barrier, 1, cmd);
    }

    // 3. Normalize normals
    {
        ClothSimConstants cb = {};
        cb.numParticles = numParticles;

        device->BindComputeShader(&normalizeNormalsCS, cmd);
        device->BindDynamicConstantBuffer(cb, CLOTH_CB_SLOT, cmd);
        const GPUResource* uav_normals[] = { &normalsBuffer };
        device->BindUAVs(uav_normals, CLOTH_SLOT_NORMALS, 1, cmd);
        const GPUResource* uav_normAccum[] = { &normAccumBuffer };
        device->BindUAVs(uav_normAccum, CLOTH_SLOT_NORM_ACCUM, 1, cmd);
        device->Dispatch(dispatchParticles, 1, 1, cmd);
        device->Barrier(&barrier, 1, cmd);
    }
}

void ClothMesh::Reset()
{
    if (!gpuBuffersReady)
        return;

    GraphicsDevice* device = wi::graphics::GetDevice();
    CommandList cmd = device->BeginCommandList();

    // Restore positions
    device->UpdateBuffer(&posBuffer, restPos.data(), cmd,
        sizeof(XMFLOAT4) * numParticles);
    device->UpdateBuffer(&prevPosBuffer, restPos.data(), cmd,
        sizeof(XMFLOAT4) * numParticles);

    // Zero velocities
    std::vector<XMFLOAT4> zeroVel(numParticles, XMFLOAT4(0, 0, 0, 0));
    device->UpdateBuffer(&velBuffer, zeroVel.data(), cmd,
        sizeof(XMFLOAT4) * numParticles);

    // Restore invMass
    device->UpdateBuffer(&invMassBuffer, cpuInvMass.data(), cmd,
        sizeof(float) * numParticles);

    device->SubmitCommandLists();

    dragParticleNr = -1;

    // Reset CPU pos
    cpuPos = restPos;
}

// ==============================================================================
// GPU RAYCAST FOR CLOTH GRABBING
// ==============================================================================
//
// WHY GPU INSTEAD OF CPU?
// -----------------------
// With a 500×500 cloth grid we have ~498,002 triangles. A CPU raycast must
// iterate over every triangle sequentially, performing a Möller-Trumbore
// intersection test for each one. This is extremely slow and causes visible
// stuttering on every click.
// Moreover, a CPU raycast would operate on the CPU copy of the cloth vertex data,
// which may be out of sync with the GPU simulation results, leading to inaccurate picking.
// To ensure accurate picking, we would need to read back the entire vertex buffer from the
// GPU to the CPU before raycasting, which would introduce additional latency and performance issues.
//
// The GPU version dispatches one compute-shader thread per triangle
// (cloth_raycastTriangleCS), so all ~500k intersection tests run in parallel
// in a single dispatch. The CPU then only needs to scan the small readback
// array of per-triangle distances to find the minimum — a trivial O(N) pass
// on already-computed results.
//
// HOW IT WORKS
// ------------------------
// Start a "drag" interaction using a picking ray.
//
// INPUTS:
// - rayOrig: ray origin (3D point), typically the camera position in the same space as the cloth vertices (object space)
// - rayDir:  ray direction (unit vector), pointing from the camera through the mouse cursor pixel
//
// HIGH-LEVEL IDEA:
// 1) Raycast the ray against ALL cloth triangles on the GPU (parallel).
// 2) Copy the per-triangle hit distances back to CPU.
// 3) Pick the triangle with the smallest distance (closest hit).
// 4) Choose one vertex of that triangle as the "dragged particle".
// 5) Temporarily set that particle invMass to 0 so the solver won't move it (see cloth_solveConstraintsCS).
// 6) Move that particle to the ray point at the hit depth (orig + depth * dir) (see cloth_integrateCS).
// ------------------------
bool ClothMesh::StartGrabGPU(const XMFLOAT3& rayOrigin, const XMFLOAT3& rayDir)
{
    if (!gpuBuffersReady)
        return false;

    GraphicsDevice* device = wi::graphics::GetDevice();
    CommandList cmd = device->BeginCommandList();

    // Retrieve ray parameters
    RaycastConstants rayCB = {};
    rayCB.origX = rayOrigin.x; rayCB.origY = rayOrigin.y; rayCB.origZ = rayOrigin.z;
    rayCB.dirX = rayDir.x; rayCB.dirY = rayDir.y; rayCB.dirZ = rayDir.z;

	// Bind the compute shader responsible for raycasting against cloth triangles (cloth_raycastTriangleCS)
	// and the constant buffer with the ray parameters (origin and direction).
    device->BindComputeShader(&raycastTriangleCS, cmd);
    device->BindDynamicConstantBuffer(rayCB, CLOTH_CB_RAYCAST_SLOT, cmd);

	// Bind the necessary buffers for raycasting and dispatch the compute shader,
	// which computes the distance from the ray to each triangle in parallel on the GPU.
	// At the end triDistBuffer[triNr] will contain:
	// - a finite t value if the ray hits triangle triNr
	// - a large sentinel value (1e9) if it doesn't hit
    const GPUResource* uav_pos[] = { &posBuffer };
    device->BindUAVs(uav_pos, CLOTH_SLOT_POS, 1, cmd);
    const GPUResource* uav_tri[] = { &triIdsBuffer };
    device->BindUAVs(uav_tri, CLOTH_SLOT_TRI_IDS, 1, cmd);
    const GPUResource* uav_dist[] = { &triDistBuffer };
    device->BindUAVs(uav_dist, CLOTH_SLOT_TRI_DIST, 1, cmd);
    device->Dispatch(dispatchTris, 1, 1, cmd);

	// Insert a memory barrier to ensure that all raycast results (distances) are
	// fully computed and visible in the triDistBuffer
    GPUBarrier barrier = GPUBarrier::Memory();
    device->Barrier(&barrier, 1, cmd);

	// Copy the triDistBuffer back to CPU for picking the closest hit triangle.
    GPUBarrier copyBarrier = GPUBarrier::Buffer(&triDistBuffer,
        ResourceState::UNORDERED_ACCESS, ResourceState::COPY_SRC);
    device->Barrier(&copyBarrier, 1, cmd);
    device->CopyResource(&triDistReadbackBuffer, &triDistBuffer, cmd);
    copyBarrier = GPUBarrier::Buffer(&triDistBuffer,
        ResourceState::COPY_SRC, ResourceState::UNORDERED_ACCESS);
    device->Barrier(&copyBarrier, 1, cmd);

	// Submit GPU work and wait for completion before reading back results.
    device->SubmitCommandLists();
    device->WaitForGPU();

	// Check if the readback buffer is mapped and contains valid data before accessing it.
    if (!triDistReadbackBuffer.mapped_data)
        return false;

	// Find the triangle with the minimum distance (closest hit).
    const float* dists = static_cast<const float*>(triDistReadbackBuffer.mapped_data);
    float minDist = 1e9f;
    int minTri = -1;
    for (int t = 0; t < numTris; t++)
    {
        if (dists[t] < minDist)
        {
            minDist = dists[t];
            minTri = t;
        }
    }

	// If no triangle was hit (minDist is still the large sentinel value),
	// or if the minimum distance is negative (which shouldn't happen with proper ray-triangle intersection),
	// we consider that the raycast failed and reset the drag state.
    if (minTri < 0 || minDist >= 1e8f)
    {
        dragParticleNr = -1;
        return false;
    }

    // Get first vertex of hit triangle as the dragged particle
	// (could be any of the three vertices, we just pick the first one for simplicity)
	// and store the hit depth that represents how far along the ray the hit occurred.
    dragParticleNr = triIds[3 * minTri];
    dragDepth = minDist;

	// Compute the current drag position in world space using the ray parameters and the hit depth.
	// This will be passed to the compute shader responsible for integrating particle positions via constant buffer,
	// so that the dragged particle can be moved to this position during integration without modifying its velocity
	// (see cloth_integrateCS for details on how drag is handled).
    currentDragPos.x = rayOrigin.x + dragDepth * rayDir.x;
    currentDragPos.y = rayOrigin.y + dragDepth * rayDir.y;
    currentDragPos.z = rayOrigin.z + dragDepth * rayDir.z;

	// Return true to indicate that the grab was successfully initiated with a valid hit triangle and dragged particle.
    return true;
}

void ClothMesh::DragGPU(const XMFLOAT3& rayOrigin, const XMFLOAT3& rayDir, wi::graphics::CommandList cmd)
{
    if (dragParticleNr < 0 || !gpuBuffersReady)
        return;

	// Update the current drag position based on the ray parameters and the hit depth.
	// This allows the dragged particle to follow the mouse cursor at the same depth as the initial hit,
	// creating an intuitive dragging interaction.
    currentDragPos.x = rayOrigin.x + dragDepth * rayDir.x;
    currentDragPos.y = rayOrigin.y + dragDepth * rayDir.y;
    currentDragPos.z = rayOrigin.z + dragDepth * rayDir.z;
}

void ClothMesh::EndGrabGPU(wi::graphics::CommandList cmd)
{
    if (dragParticleNr < 0 || !gpuBuffersReady)
        return;

    // Just reset the drag state; next frame shaders see dragParticleNr = -1.
    dragParticleNr = -1;
}

void ClothMesh::RequestPositionsReadback(wi::graphics::CommandList cmd)
{
    if (!gpuBuffersReady)
        return;

    GraphicsDevice* device = wi::graphics::GetDevice();

	// Insert a barrier to ensure that all GPU writes to the position buffer are
	// completed and visible before we copy it for readback.
    GPUBarrier barrier = GPUBarrier::Buffer(&posBuffer,
        ResourceState::UNORDERED_ACCESS, ResourceState::COPY_SRC);
    device->Barrier(&barrier, 1, cmd);

	// Copy the position buffer to the readback buffer so that we can access it on the CPU.
    device->CopyResource(&posReadbackBuffer, &posBuffer, cmd);

	// Insert another barrier to transition the position buffer back to
	// unordered access for the next simulation step.
    barrier = GPUBarrier::Buffer(&posBuffer,
        ResourceState::COPY_SRC, ResourceState::UNORDERED_ACCESS);
    device->Barrier(&barrier, 1, cmd);

	// Mark that we have a pending readback for positions, so that the CPU knows to process it when ready.
    posReadbackPending = true;
}

void ClothMesh::RequestNormalsReadback(wi::graphics::CommandList cmd)
{
    if (!gpuBuffersReady)
        return;

    GraphicsDevice* device = wi::graphics::GetDevice();

    GPUBarrier barrier = GPUBarrier::Buffer(&normalsBuffer,
        ResourceState::UNORDERED_ACCESS, ResourceState::COPY_SRC);
    device->Barrier(&barrier, 1, cmd);

    device->CopyResource(&normalsReadbackBuffer, &normalsBuffer, cmd);

    barrier = GPUBarrier::Buffer(&normalsBuffer,
        ResourceState::COPY_SRC, ResourceState::UNORDERED_ACCESS);
    device->Barrier(&barrier, 1, cmd);

    normalsReadbackPending = true;
}

void ClothMesh::ProcessPositionsReadback()
{
	// If the GPU buffers are not ready or if there is no pending readback for positions,
	// we cannot process the readback data.
    if (!gpuBuffersReady || !posReadbackPending)
        return;

	// Check if the readback buffer for positions is mapped and contains valid data before accessing it.
    if (posReadbackBuffer.mapped_data)
    {
		// Copy the position data from the mapped readback buffer to the CPU array (cpuPos)
		// for use in rendering and other possible CPU-side operations.
        const XMFLOAT4* src = static_cast<const XMFLOAT4*>(posReadbackBuffer.mapped_data);
        for (int i = 0; i < numParticles; i++)
            cpuPos[i] = src[i];

		// Mark the readback as processed by setting posReadbackPending to false,
		// so that we don't process the same data again.
        posReadbackPending = false;
    }
}

void ClothMesh::ProcessNormalsReadback()
{
    if (!gpuBuffersReady || !normalsReadbackPending)
        return;

    if (normalsReadbackBuffer.mapped_data)
    {
        const XMFLOAT4* src = static_cast<const XMFLOAT4*>(normalsReadbackBuffer.mapped_data);
        for (int i = 0; i < numParticles; i++)
            cpuNormals[i] = XMFLOAT3(src[i].x, src[i].y, src[i].z);

        normalsReadbackPending = false;
    }
}
