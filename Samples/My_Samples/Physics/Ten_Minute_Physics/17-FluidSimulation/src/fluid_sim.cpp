#include "fluid_sim.h"
#include <wiRenderer.h>
#include <wiHelper.h>
#include <wiGraphics.h>
#include <wiBacklog.h>
#include <cmath>
#include <thread>

using namespace wi::graphics;

// ============================================================================
// Construction / Initialization
// ============================================================================

FluidSimulation3D::FluidSimulation3D() {}

void FluidSimulation3D::Initialize(const FluidParams& inParams)
{
    params = inParams;
    SetupScene();
    CreateTextures();
    LoadShaders();
    gpuReady = true;
	resetPending = true;
}

void FluidSimulation3D::SetupScene()
{
    // Compute grid dimensions.
	// Domain is a cube of size params.domainSize which represents the physical space where the grid lives.
    // Each cell has size computed as domainSize / gridX, so the number of cells is just gridX in each dimension.
    gridX = params.resolution;
    gridY = params.resolution;
    gridZ = params.resolution;
    cellSize = params.domainSize / (float)gridX;

    // Domain is centered at (0, domainSize/2, 0) so the bottom sits on the ground plane
    float halfSize = params.domainSize * 0.5f;
    domainMin = XMFLOAT3(-halfSize, 0.0f, -halfSize);
    domainMax = XMFLOAT3(halfSize, params.domainSize, halfSize);
}

// ============================================================================
// GPU Texture Creation
// ============================================================================
//
// All simulation fields (velocity components, pressure, smoke density, solid mask, temperature, curl)
// are stored in 3D textures. We use UAVs for read/write access in compute shaders,
// and SRVs for read-only access when needed (e.g., advection).
void FluidSimulation3D::CreateTextures()
{
    GraphicsDevice* device = GetDevice();

    // Helper: create a 3D texture with UAV + SRV
    auto make3D = [&](Texture& tex, Format format, int w, int h, int d, const char* name)
    {
        TextureDesc desc;
        desc.type = TextureDesc::Type::TEXTURE_3D;
        desc.width = (uint32_t)w;
        desc.height = (uint32_t)h;
        desc.depth = (uint32_t)d;
        desc.mip_levels = 1;
        desc.format = format;
        desc.bind_flags = BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;
        desc.usage = Usage::DEFAULT;
        device->CreateTexture(&desc, nullptr, &tex);
        device->SetName(&tex, name);
    };

    // Velocity fields (staggered MAC) store components at cell boundaries:
    // u follows the x direction, so it has one extra column in that direction to store all u velocities.
	// Immagine the grid having 2 cells in x, then we need 3 u values to store:
	//  - the u velocity at the left boundary of the first cell
	//  - the u velocity at the interface between the first and second cell
	//  - the u velocity at the right boundary of the second cell
	// Similarly for v and w.
	// temp textures are used for ping-ponging during advection (see advectVelCS and advectSmokeCS).
    make3D(texU,      Format::R32_FLOAT, gridX + 1, gridY,     gridZ,     "fluid::u");
    make3D(texV,      Format::R32_FLOAT, gridX,     gridY + 1, gridZ,     "fluid::v");
    make3D(texW,      Format::R32_FLOAT, gridX,     gridY,     gridZ + 1, "fluid::w");
    make3D(texU_temp, Format::R32_FLOAT, gridX + 1, gridY,     gridZ,     "fluid::u_temp");
    make3D(texV_temp, Format::R32_FLOAT, gridX,     gridY + 1, gridZ,     "fluid::v_temp");
    make3D(texW_temp, Format::R32_FLOAT, gridX,     gridY,     gridZ + 1, "fluid::w_temp");

    // Cell-centered scalar fields
    make3D(texPressure,   Format::R32_FLOAT, gridX, gridY, gridZ, "fluid::pressure");
    make3D(texDivergence, Format::R32_FLOAT, gridX, gridY, gridZ, "fluid::divergence");
    make3D(texSmoke,      Format::R16_FLOAT, gridX, gridY, gridZ, "fluid::smoke");
    make3D(texSmoke_temp, Format::R16_FLOAT, gridX, gridY, gridZ, "fluid::smoke_temp");

    // Solid mask: 0 = solid, 1 = fluid
    // R8_UINT is sufficient (only stores 0 or 1, no mathematical operations)
    make3D(texSolid, Format::R8_UINT, gridX, gridY, gridZ, "fluid::solid");

    // Temperature field (cell-centered, for buoyancy)
    make3D(texTemperature,      Format::R16_FLOAT, gridX, gridY, gridZ, "fluid::temperature");
    make3D(texTemperature_temp, Format::R16_FLOAT, gridX, gridY, gridZ, "fluid::temperature_temp");

    // Curl / vorticity vector (cell-centered, xyz = curl components)
    make3D(texCurl, Format::R16G16B16A16_FLOAT, gridX, gridY, gridZ, "fluid::curl");
}

// ============================================================================
// Shader Loading
// ============================================================================

void FluidSimulation3D::LoadShaders()
{
    // Wait for all engine background shader compilation jobs to finish before
    // temporarily changing global shader paths.  Engine jobs (e.g. objectPS
    // permutations) read SHADERSOURCEPATH at execution time — if we change it
    // while they're still in the queue, they'll look for engine .hlsl files in
    // our sample's shaders/ directory and fail.
    while (wi::renderer::IsPipelineCreationActive() > 0)
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }

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

    using SM = ShaderStage;

    // Simulation kernels — all from fluid_simulation.hlsl with permutation defines
    wi::renderer::LoadShader(SM::CS, applyForcesCS,   "fluid_simulation.cso", ShaderModel::SM_6_0, {"APPLY_FORCES"},     "applyForcesCS");
    wi::renderer::LoadShader(SM::CS, setObstacleCS,   "fluid_simulation.cso", ShaderModel::SM_6_0, {"SET_OBSTACLE"},     "setObstacleCS");
    wi::renderer::LoadShader(SM::CS, divergenceCS,    "fluid_simulation.cso", ShaderModel::SM_6_0, {"DIVERGENCE"},       "divergenceCS");
    wi::renderer::LoadShader(SM::CS, pressureRedCS,   "fluid_simulation.cso", ShaderModel::SM_6_0, {"PRESSURE_RED"},     "pressureRedCS");
    wi::renderer::LoadShader(SM::CS, pressureBlackCS, "fluid_simulation.cso", ShaderModel::SM_6_0, {"PRESSURE_BLACK"},   "pressureBlackCS");
    wi::renderer::LoadShader(SM::CS, projectCS,       "fluid_simulation.cso", ShaderModel::SM_6_0, {"PROJECT"},          "projectCS");
    wi::renderer::LoadShader(SM::CS, advectVelCS,     "fluid_simulation.cso", ShaderModel::SM_6_0, {"ADVECT_VELOCITY"},  "advectVelCS");
    wi::renderer::LoadShader(SM::CS, advectSmokeCS,   "fluid_simulation.cso", ShaderModel::SM_6_0, {"ADVECT_SMOKE"},     "advectSmokeCS");
    wi::renderer::LoadShader(SM::CS, boundaryCS,      "fluid_simulation.cso", ShaderModel::SM_6_0, {"BOUNDARY"},         "boundaryCS");
    wi::renderer::LoadShader(SM::CS, computeCurlCS,   "fluid_simulation.cso", ShaderModel::SM_6_0, {"COMPUTE_CURL"},     "computeCurlCS");
    wi::renderer::LoadShader(SM::CS, applyVorticityCS,"fluid_simulation.cso", ShaderModel::SM_6_0, {"APPLY_VORTICITY"},  "applyVorticityCS");
    wi::renderer::LoadShader(SM::CS, advectTempCS,    "fluid_simulation.cso", ShaderModel::SM_6_0, {"ADVECT_TEMPERATURE"},"advectTemperatureCS");

    // MacCormack correction kernels
    wi::renderer::LoadShader(SM::CS, macCormackVelCS,   "fluid_simulation.cso", ShaderModel::SM_6_0, {"MACCORMACK_VEL"},   "macCormackVelCS");
    wi::renderer::LoadShader(SM::CS, macCormackSmokeCS, "fluid_simulation.cso", ShaderModel::SM_6_0, {"MACCORMACK_SMOKE"}, "macCormackSmokeCS");
    wi::renderer::LoadShader(SM::CS, macCormackTempCS,  "fluid_simulation.cso", ShaderModel::SM_6_0, {"MACCORMACK_TEMP"},  "macCormackTempCS");

    // Ray-march render shader (separate file)
    wi::renderer::LoadShader(SM::CS, raymarchCS, "fluid_raymarching.cso");

    // Restore original paths
    wi::renderer::SetShaderPath(originalBinPath);
    wi::renderer::SetShaderSourcePath(originalSrcPath);
}

// ============================================================================
// Simulation Dispatch
// ============================================================================

void FluidSimulation3D::DispatchSim3D(CommandList cmd, int cellsX, int cellsY, int cellsZ)
{
    GraphicsDevice* device = GetDevice();
    uint32_t gx = ((uint32_t)cellsX + FLUID_THREADS_3D_X - 1) / FLUID_THREADS_3D_X;
    uint32_t gy = ((uint32_t)cellsY + FLUID_THREADS_3D_Y - 1) / FLUID_THREADS_3D_Y;
    uint32_t gz = ((uint32_t)cellsZ + FLUID_THREADS_3D_Z - 1) / FLUID_THREADS_3D_Z;
    device->Dispatch(gx, gy, gz, cmd);
}

void FluidSimulation3D::SetObstacle(const XMFLOAT3& center, const XMFLOAT3& velocity)
{
    params.obsCenter = center;
    obsVelocity = velocity;
}

void FluidSimulation3D::SetSmokeSrc(const XMFLOAT3& center)
{
    params.smokeSrcCenter = center;
}

void FluidSimulation3D::Simulate(CommandList cmd)
{
    if (!gpuReady)
        return;
    if (resetPending)
    {
        ResetState(cmd);
    }
    if (paused)
        return;

    GraphicsDevice* device = GetDevice();
    GPUBarrier barrier = GPUBarrier::Memory();

    // Build the constant buffer shared by all simulation kernels
    FluidConstants cb = {};
    cb.gridX = gridX;
    cb.gridY = gridY;
    cb.gridZ = gridZ;
    cb.cellSize = cellSize;
    cb.dt = params.dt;
    cb.density = params.density;
    cb.overRelaxation = params.overRelaxation;
    cb.numPressureIters = params.numPressureIters;
    cb.gravX = params.gravity.x;
    cb.gravY = params.gravity.y;
    cb.gravZ = params.gravity.z;
    cb.obsCenterX = params.obsCenter.x;
    cb.obsCenterY = params.obsCenter.y;
    cb.obsCenterZ = params.obsCenter.z;
    cb.obsRadius = params.obsRadius;
    cb.obsVelX = obsVelocity.x;
    cb.obsVelY = obsVelocity.y;
    cb.obsVelZ = obsVelocity.z;
    cb.buoyancy = params.buoyancy;
    cb.enableSmoke = 1;

    // Smoke source: draggable emitter sphere
    cb.smokeSrcX = params.smokeSrcCenter.x;
    cb.smokeSrcY = params.smokeSrcCenter.y;
    cb.smokeSrcZ = params.smokeSrcCenter.z;
    cb.smokeSrcRadius = params.smokeSrcRadius;
    cb.smokeInflowDensity = params.smokeInflowDensity;
    cb.smokeDissipation = params.dissipation;
    cb.vorticityStrength = params.vorticityStrength;
    cb.tempInjection = params.tempInjection;
    cb.ambientTemp = params.ambientTemp;
    cb.fluidWeight = params.fluidWeight;
    cb.velDissipation = params.velDissipation;

    // Helper: bind base simulation UAVs (u0-u6)
    auto bindBaseUAVs = [&]()
    {
        const GPUResource* uavs[] = {
            &texU, &texV, &texW, &texPressure,
            &texDivergence, &texSmoke, &texSolid
        };
        device->BindUAVs(uavs, FLUID_SLOT_VELOCITY_U, 7, cmd);
    };

    // Helper: bind temperature UAVs (u11-u12)
    auto bindTempUAVs = [&]()
    {
        const GPUResource* uavs[] = { &texTemperature, &texTemperature_temp };
        device->BindUAVs(uavs, FLUID_SLOT_TEMPERATURE, 2, cmd);
    };

    // Helper: bind curl UAV (u13)
    auto bindCurlUAV = [&]()
    {
        const GPUResource* uavs[] = { &texCurl };
        device->BindUAVs(uavs, FLUID_SLOT_CURL, 1, cmd);
    };

    // --- STEP 1: Set obstacle (mark solid cells, inject velocity) ---
    {
        device->BindComputeShader(&setObstacleCS, cmd);
        device->BindDynamicConstantBuffer(cb, FLUID_CB_SLOT, cmd);
        bindBaseUAVs();
        DispatchSim3D(cmd, gridX, gridY, gridZ);
        device->Barrier(&barrier, 1, cmd);
    }

    // --- STEP 2: Apply forces (gravity + temp buoyancy + weight + injection) ---
    {
        device->BindComputeShader(&applyForcesCS, cmd);
        device->BindDynamicConstantBuffer(cb, FLUID_CB_SLOT, cmd);
        bindBaseUAVs();
        bindTempUAVs();
        DispatchSim3D(cmd, gridX, gridY, gridZ);
        device->Barrier(&barrier, 1, cmd);
    }

    // --- STEP 3: Compute curl (vorticity ∇×v) ---
    {
        device->BindComputeShader(&computeCurlCS, cmd);
        device->BindDynamicConstantBuffer(cb, FLUID_CB_SLOT, cmd);
        bindBaseUAVs();
        bindCurlUAV();
        DispatchSim3D(cmd, gridX, gridY, gridZ);
        device->Barrier(&barrier, 1, cmd);
    }

    // --- STEP 4: Apply vorticity confinement (velocity-only swirl reinforcement, not MacCormack) ---
    {
        device->BindComputeShader(&applyVorticityCS, cmd);
        device->BindDynamicConstantBuffer(cb, FLUID_CB_SLOT, cmd);
        bindBaseUAVs();
        bindCurlUAV();
        DispatchSim3D(cmd, gridX, gridY, gridZ);
        device->Barrier(&barrier, 1, cmd);
    }

    // --- STEP 5: Boundary conditions (pre-projection) ---
    {
        device->BindComputeShader(&boundaryCS, cmd);
        device->BindDynamicConstantBuffer(cb, FLUID_CB_SLOT, cmd);
        bindBaseUAVs();
        DispatchSim3D(cmd, gridX, gridY, gridZ);
        device->Barrier(&barrier, 1, cmd);
    }

    // --- STEP 6: Compute divergence ---
    {
        device->BindComputeShader(&divergenceCS, cmd);
        device->BindDynamicConstantBuffer(cb, FLUID_CB_SLOT, cmd);
        bindBaseUAVs();
        DispatchSim3D(cmd, gridX, gridY, gridZ);
        device->Barrier(&barrier, 1, cmd);
    }

    // --- STEP 7: Pressure solve (Red-Black Gauss-Seidel + SOR) ---
    for (int iter = 0; iter < params.numPressureIters; iter++)
    {
        cb.redBlackPass = 0;
        device->BindComputeShader(&pressureRedCS, cmd);
        device->BindDynamicConstantBuffer(cb, FLUID_CB_SLOT, cmd);
        bindBaseUAVs();
        DispatchSim3D(cmd, gridX, gridY, gridZ);
        device->Barrier(&barrier, 1, cmd);

        cb.redBlackPass = 1;
        device->BindComputeShader(&pressureBlackCS, cmd);
        device->BindDynamicConstantBuffer(cb, FLUID_CB_SLOT, cmd);
        bindBaseUAVs();
        DispatchSim3D(cmd, gridX, gridY, gridZ);
        device->Barrier(&barrier, 1, cmd);
    }

    // --- STEP 8: Project (subtract pressure gradient from velocity to make it divergence-free)  ---
    {
        device->BindComputeShader(&projectCS, cmd);
        device->BindDynamicConstantBuffer(cb, FLUID_CB_SLOT, cmd);
        bindBaseUAVs();
        DispatchSim3D(cmd, gridX, gridY, gridZ);
        device->Barrier(&barrier, 1, cmd);
    }

    // --- STEP 9: Boundary conditions (post-projection) ---
    {
        device->BindComputeShader(&boundaryCS, cmd);
        device->BindDynamicConstantBuffer(cb, FLUID_CB_SLOT, cmd);
        bindBaseUAVs();
        DispatchSim3D(cmd, gridX, gridY, gridZ);
        device->Barrier(&barrier, 1, cmd);
    }

    // --- STEP 10: Advect velocity (forward Semi-Lagrangian) ---
    {
        device->BindComputeShader(&advectVelCS, cmd);
        device->BindDynamicConstantBuffer(cb, FLUID_CB_SLOT, cmd);
        bindBaseUAVs();
        const GPUResource* uavs_dst[] = {
            &texU_temp, &texV_temp, &texW_temp, &texSmoke_temp
        };
        device->BindUAVs(uavs_dst, FLUID_SLOT_VELOCITY_U_TEMP, 4, cmd);
        int maxDim = std::max({gridX + 1, gridY + 1, gridZ + 1});
        DispatchSim3D(cmd, maxDim, maxDim, maxDim);
        device->Barrier(&barrier, 1, cmd);
    }

    // --- STEP 11: Advect smoke (forward Semi-Lagrangian) ---
    {
        device->BindComputeShader(&advectSmokeCS, cmd);
        device->BindDynamicConstantBuffer(cb, FLUID_CB_SLOT, cmd);
        bindBaseUAVs();
        const GPUResource* uavs_dst[] = {
            &texU_temp, &texV_temp, &texW_temp, &texSmoke_temp
        };
        device->BindUAVs(uavs_dst, FLUID_SLOT_VELOCITY_U_TEMP, 4, cmd);
        DispatchSim3D(cmd, gridX, gridY, gridZ);
        device->Barrier(&barrier, 1, cmd);
    }

    // --- STEP 12: Advect temperature (forward Semi-Lagrangian) ---
    {
        device->BindComputeShader(&advectTempCS, cmd);
        device->BindDynamicConstantBuffer(cb, FLUID_CB_SLOT, cmd);
        bindBaseUAVs();
        bindTempUAVs();
        DispatchSim3D(cmd, gridX, gridY, gridZ);
        device->Barrier(&barrier, 1, cmd);
    }

    // --- STEP 13: MacCormack correction (optional) ---
    // Corrections read the stable forward-SL fields in _temp (read-only) and
    // write corrected results back to the ORIGINAL textures.
    // When MacCormack is ON, skip CopyResource (originals already correct).
    // When MacCormack is OFF, CopyResource moves forward SL results from _temp → originals.
    if (useMacCormack)
    {
        // 13a: Smoke correction — same logic as macCormackSmokeCS:
        // reads original smoke + forward-SL smoke in texSmoke_temp,
        // then writes corrected smoke to texSmoke.
        {
            device->BindComputeShader(&macCormackSmokeCS, cmd);
            device->BindDynamicConstantBuffer(cb, FLUID_CB_SLOT, cmd);
            bindBaseUAVs();
            const GPUResource* uavs_dst[] = {
                &texU_temp, &texV_temp, &texW_temp, &texSmoke_temp
            };
            device->BindUAVs(uavs_dst, FLUID_SLOT_VELOCITY_U_TEMP, 4, cmd);
            DispatchSim3D(cmd, gridX, gridY, gridZ);
            device->Barrier(&barrier, 1, cmd);
        }

        // 13b: Temperature correction — same scheme as smoke, but without a
        // separate per-step dissipation term during temperature advection.
        {
            device->BindComputeShader(&macCormackTempCS, cmd);
            device->BindDynamicConstantBuffer(cb, FLUID_CB_SLOT, cmd);
            bindBaseUAVs();
            bindTempUAVs();
            DispatchSim3D(cmd, gridX, gridY, gridZ);
            device->Barrier(&barrier, 1, cmd);
        }

        // 13c: Velocity correction — same MacCormack idea, but on the
        // staggered velocity field. Traces and phi_bwd are reconstructed from
        // the already-advected _temp velocity field, then corrected values are
        // written to texU / texV / texW.
        {
            device->BindComputeShader(&macCormackVelCS, cmd);
            device->BindDynamicConstantBuffer(cb, FLUID_CB_SLOT, cmd);
            bindBaseUAVs();
            const GPUResource* uavs_dst[] = {
                &texU_temp, &texV_temp, &texW_temp, &texSmoke_temp
            };
            device->BindUAVs(uavs_dst, FLUID_SLOT_VELOCITY_U_TEMP, 4, cmd);
            int maxDim = std::max({gridX + 1, gridY + 1, gridZ + 1});
            DispatchSim3D(cmd, maxDim, maxDim, maxDim);
            device->Barrier(&barrier, 1, cmd);
        }
        // No CopyResource needed — corrections wrote directly to originals.
    }
    else
    {
        // --- STEP 14: Copy _temp → originals (pure SL path) ---
        device->CopyResource(&texU, &texU_temp, cmd);
        device->CopyResource(&texV, &texV_temp, cmd);
        device->CopyResource(&texW, &texW_temp, cmd);
        device->CopyResource(&texSmoke, &texSmoke_temp, cmd);
        device->CopyResource(&texTemperature, &texTemperature_temp, cmd);
        device->Barrier(&barrier, 1, cmd);
    }

    // --- STEP 15: Boundary conditions (post-advection) ---
    // Advection can pollute velocities at solid faces; re-enforce BCs.
    {
        device->BindComputeShader(&boundaryCS, cmd);
        device->BindDynamicConstantBuffer(cb, FLUID_CB_SLOT, cmd);
        bindBaseUAVs();
        DispatchSim3D(cmd, gridX, gridY, gridZ);
        device->Barrier(&barrier, 1, cmd);
    }
}

// ============================================================================
// Volume Rendering (Ray Marching)
//
// The simulation step updates 3D scalar/vector fields (smoke, pressure,
// velocity), but those textures "only" hold raw simulation data and don't have a direct
// visual representation. They represent the physical state of the fluid,
// but to visualize it we need to convert that data into something we can draw on the screen.
// The 3D textures themselves are not meant to be rendered directly; they are more like
// "physics data containers" for the simulation. They are not directly drawable like
// ordinary scene meshes. To visualize the smoke from the current camera view, we run a
// canvas-sized compute pass that ray-marches the simulation domain and
// converts the sampled volume into a 2D image (renderOutput).
//
// In other words:
//   simulation 3D textures -> per-pixel ray marching -> 2D overlay texture
//
// This is a standard volumetric rendering approach for smoke/fog because the
// fluid is a participating medium distributed through space, not a closed
// triangle surface that could be rendered with the usual mesh pipeline.
//
// The resulting renderOutput is composited later in SampleRenderPath::Compose()
// as a 2D overlay on top of the already rendered scene.
//
// Current limitation:
// this sample does not make the ray-march pass read the scene depth buffer.
// Because of that, the fluid overlay has no per-pixel knowledge of where the
// already rendered opaque scene should stop the volume along the current view
// ray.
// Practical example:
// suppose a wall is in front of part of the smoke volume. The normal 3D scene
// pass renders that wall first. Later, the ray-march pass still computes smoke
// color for the same screen pixels and writes it into renderOutput. During
// Compose(), that 2D fluid image is simply blended on top of the scene image.
// Since the ray march never checked the wall depth, it cannot know that the
// volume contribution for those pixels should have been clipped by the wall.
// Workaround:
// The current workaround is to skip the overlay entirely when the camera is
// inside the domain. A more robust solution would be a depth-aware ray march:
// sample the scene depth buffer for each pixel, convert that depth to a ray
// distance, and stop marching when the ray reaches the nearest opaque surface.
// That would let the volume respect scene occlusion instead of relying on the
// coarse workaround of disabling the fluid when the camera is inside the
// domain.
// ============================================================================

void FluidSimulation3D::RenderVolume(CommandList cmd, const wi::scene::CameraComponent& camera,
                                     uint32_t screenWidth, uint32_t screenHeight)
{
    renderOverlay = false;  // default: don't draw overlay

    if (!gpuReady)
        return;

    GraphicsDevice* device = GetDevice();

    uint32_t screenW = screenWidth;
    uint32_t screenH = screenHeight;

    if (screenW == 0 || screenH == 0)
        return;

    // Skip the overlay when the camera is inside the domain.
    // See the note above: without depth-aware ray marching, composing the
    // volume as a 2D overlay produces incorrect in-front-of-camera artifacts.
    XMFLOAT3 camPos;
    XMStoreFloat3(&camPos, camera.GetEye());
    if (camPos.x > domainMin.x && camPos.x < domainMax.x &&
        camPos.y > domainMin.y && camPos.y < domainMax.y &&
        camPos.z > domainMin.z && camPos.z < domainMax.z)
    {
        return;
    }

    renderOverlay = true;  // camera is outside, we will render

    if (!renderOutput.IsValid() ||
        renderOutput.GetDesc().width != screenW ||
        renderOutput.GetDesc().height != screenH)
    {
        TextureDesc desc;
        desc.type = TextureDesc::Type::TEXTURE_2D;
        desc.width = screenW;
        desc.height = screenH;
        desc.mip_levels = 1;
        desc.format = Format::R16G16B16A16_FLOAT;
        desc.bind_flags = BindFlag::SHADER_RESOURCE | BindFlag::UNORDERED_ACCESS;
        desc.usage = Usage::DEFAULT;
        device->CreateTexture(&desc, nullptr, &renderOutput);
        device->SetName(&renderOutput, "fluid::renderOutput");
    }

    // Build render constants for the canvas-sized ray-march pass:
    // domain extents, camera data for ray reconstruction, current visualization
    // mode, and a few helper parameters for obstacle/source visualization.
    FluidRenderConstants rc = {};
    rc.gridX = gridX;
    rc.gridY = gridY;
    rc.gridZ = gridZ;
    rc.cellSize = cellSize;
    rc.domainMinX = domainMin.x;
    rc.domainMinY = domainMin.y;
    rc.domainMinZ = domainMin.z;
    rc.domainMaxX = domainMax.x;
    rc.domainMaxY = domainMax.y;
    rc.domainMaxZ = domainMax.z;
    rc.renderMode = renderMode;
    rc.maxSteps = 256.0f;
    rc.screenWidth = (int)screenW;
    rc.screenHeight = (int)screenH;
    rc.smokeAbsorption = smokeAbsorption;
    rc.obsCenterX = params.obsCenter.x;
    rc.obsCenterY = params.obsCenter.y;
    rc.obsCenterZ = params.obsCenter.z;
    rc.obsRadius = params.obsRadius;

    // Camera world position
    XMFLOAT3 camEye;
    XMStoreFloat3(&camEye, camera.GetEye());
    rc.cameraPosX = camEye.x;
    rc.cameraPosY = camEye.y;
    rc.cameraPosZ = camEye.z;

    // Smoke source sphere (for rendering a faint indicator)
    rc.smokeSrcRenderX = params.smokeSrcCenter.x;
    rc.smokeSrcRenderY = params.smokeSrcCenter.y;
    rc.smokeSrcRenderZ = params.smokeSrcCenter.z;
    rc.smokeSrcRenderRadius = params.smokeSrcRadius;

    // Configurable smoke color
    rc.smokeColorR = params.smokeColorR;
    rc.smokeColorG = params.smokeColorG;
    rc.smokeColorB = params.smokeColorB;

    // Store inverse VP matrix
    XMMATRIX VP = camera.GetViewProjection();
    XMMATRIX invVP = XMMatrixInverse(nullptr, VP);
    XMStoreFloat4x4(reinterpret_cast<XMFLOAT4X4*>(rc.invVP), invVP);

    device->BindComputeShader(&raymarchCS, cmd);
    device->BindDynamicConstantBuffer(rc, FLUID_CB_SLOT, cmd);

    // Bind the simulation 3D textures as read-only inputs for the render pass.
    // The compute shader samples these volumes along each camera ray and writes
    // a single 2D RGBA result per screen pixel.
    const GPUResource* srvs[] = {
        &texSmoke, &texPressure, &texU, &texV, &texW
    };
    device->BindResources(srvs, FLUID_SRV_SMOKE, 5, cmd);

    // Bind render output as UAV
    const GPUResource* uavs[] = { &renderOutput };
    device->BindUAVs(uavs, FLUID_UAV_RENDER_OUTPUT, 1, cmd);

    // Dispatch screen-sized thread groups for ray marching.
	// Each thread computes one pixel's ray march, so we need enough threads to cover the entire screen.
    uint32_t gx = (screenW + FLUID_THREADS_2D - 1) / FLUID_THREADS_2D;
    uint32_t gy = (screenH + FLUID_THREADS_2D - 1) / FLUID_THREADS_2D;
    device->Dispatch(gx, gy, 1, cmd);

    GPUBarrier barrier = GPUBarrier::Memory();
    device->Barrier(&barrier, 1, cmd);
}

// ============================================================================
// Reset
// ============================================================================

void FluidSimulation3D::Reset()
{
	resetPending = true;
}

void FluidSimulation3D::ResetState(CommandList cmd)
{
    GraphicsDevice* device = GetDevice();

    // Clear every simulation volume explicitly so the paused renderer and the
    // next simulation step both see a deterministic all-zero state.
    device->ClearUAV(&texU, 0, cmd);
    device->ClearUAV(&texV, 0, cmd);
    device->ClearUAV(&texW, 0, cmd);
    device->ClearUAV(&texU_temp, 0, cmd);
    device->ClearUAV(&texV_temp, 0, cmd);
    device->ClearUAV(&texW_temp, 0, cmd);
    device->ClearUAV(&texPressure, 0, cmd);
    device->ClearUAV(&texDivergence, 0, cmd);
    device->ClearUAV(&texSmoke, 0, cmd);
    device->ClearUAV(&texSmoke_temp, 0, cmd);
    device->ClearUAV(&texSolid, 0, cmd);
    device->ClearUAV(&texTemperature, 0, cmd);
    device->ClearUAV(&texTemperature_temp, 0, cmd);
    device->ClearUAV(&texCurl, 0, cmd);

    if (renderOutput.IsValid())
    {
        device->ClearUAV(&renderOutput, 0, cmd);
    }

    GPUBarrier barrier = GPUBarrier::Memory();
    device->Barrier(&barrier, 1, cmd);

    obsVelocity = {0, 0, 0};
    renderOverlay = false;
    resetPending = false;
}
