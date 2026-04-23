#pragma once

#include <wiGraphics.h>
#include <wiScene.h>
#include "shader_interop_fluid.h"

// ============================================================================
// FluidSimulation3D
//
// 3D Eulerian fluid simulator on GPU using compute shaders.
// Single mode: buoyant smoke rising from a draggable spherical emitter,
// with a draggable solid obstacle sphere.
// ============================================================================

struct FluidParams
{
    int   resolution    = 128;
    float domainSize    = 4.0f;
    float dt            = 1.0f / 60.0f;
    float density       = 1000.0f;
    float overRelaxation = 1.9f;
    int   numPressureIters = 40;
    XMFLOAT3 gravity    = {0.0f, 0.0f, 0.0f};

    // Solid obstacle sphere
    XMFLOAT3 obsCenter  = {0.0f, 2.0f, 0.0f};
    float    obsRadius  = 0.3f;

    // Smoke emitter sphere (draggable)
    XMFLOAT3 smokeSrcCenter = {0.0f, 0.5f, 0.0f};
    float    smokeSrcRadius = 0.35f;
    float    smokeInflowDensity = 1.0f;

    // Buoyancy & dissipation (tuned for realistic smoke plume)
    float buoyancy      = 1.0f;     // upward force per unit temperature difference
    float dissipation   = 0.999f;   // smoke kept per advection step (multiplicative)

    // Temperature-based buoyancy
    float tempInjection  = 12.0f;    // temperature injected at source per frame
    float ambientTemp    = 0.0f;     // reference temperature
    float fluidWeight    = 0.125f;   // weight pulling dense smoke down

    // Vorticity confinement (stronger = more turbulent detail)
    float vorticityStrength = 0.35f;

    // Velocity dissipation
    float velDissipation = 0.9995f;  // velocity kept per step (1.0 = no loss)

    // Smoke color for ray-march rendering
    float smokeColorR = 0.92f;
    float smokeColorG = 0.90f;
    float smokeColorB = 0.85f;
};

class FluidSimulation3D
{
public:
    FluidSimulation3D();

    void Initialize(const FluidParams& params);
    void Simulate(wi::graphics::CommandList cmd);
    void RenderVolume(wi::graphics::CommandList cmd, const wi::scene::CameraComponent& camera,
                      uint32_t screenWidth, uint32_t screenHeight);
    void Reset();

    // Call when obstacle moves (e.g., mouse drag)
    void SetObstacle(const XMFLOAT3& center, const XMFLOAT3& velocity);
    // Call when smoke source moves (e.g., mouse drag)
    void SetSmokeSrc(const XMFLOAT3& center);

    FluidParams params;

    // Render mode: 0=smoke, 1=pressure, 2=velocity
    int renderMode = 0;
    float smokeAbsorption = 80.0f;
    bool paused = true;
    bool useMacCormack = true;  // MacCormack advection (sharper, less diffusive)

    // Grid dimensions (computed from resolution)
    int gridX = 0, gridY = 0, gridZ = 0;
    float cellSize = 0.0f;
    XMFLOAT3 domainMin = {0, 0, 0};
    XMFLOAT3 domainMax = {0, 0, 0};

    // Render output texture (ray-marched result)
    wi::graphics::Texture renderOutput;
    bool renderOverlay = false;   // set by RenderVolume, checked by Compose

    bool IsReady() const { return gpuReady; }

private:
    void CreateTextures();
    void LoadShaders();
    void SetupScene();
    void ResetState(wi::graphics::CommandList cmd);


    // Dispatch helpers
    void DispatchSim3D(wi::graphics::CommandList cmd, int cellsX, int cellsY, int cellsZ);

    // Obstacle state
    XMFLOAT3 obsVelocity = {0, 0, 0};

    // =========================================================================
    // 3D Textures (simulation fields)
    // =========================================================================
    // Velocity components (staggered MAC grid)
    wi::graphics::Texture texU, texV, texW;             // current
    wi::graphics::Texture texU_temp, texV_temp, texW_temp; // ping-pong for advection
    // Scalar fields (cell-centered)
    wi::graphics::Texture texPressure;
    wi::graphics::Texture texDivergence;
    wi::graphics::Texture texSmoke, texSmoke_temp;
    wi::graphics::Texture texSolid;                     // R8_UINT: 0=solid, 1=fluid
    // Temperature (cell-centered scalar for buoyancy)
    wi::graphics::Texture texTemperature, texTemperature_temp;
    // Curl / vorticity (cell-centered float4: xyz = curl)
    wi::graphics::Texture texCurl;

    // =========================================================================
    // Compute shaders (loaded from fluid_simulation.hlsl with permutation defines)
    // =========================================================================
    wi::graphics::Shader applyForcesCS;
    wi::graphics::Shader setObstacleCS;
    wi::graphics::Shader divergenceCS;
    wi::graphics::Shader pressureRedCS;
    wi::graphics::Shader pressureBlackCS;
    wi::graphics::Shader projectCS;
    wi::graphics::Shader advectVelCS;
    wi::graphics::Shader advectSmokeCS;
    wi::graphics::Shader boundaryCS;
    wi::graphics::Shader computeCurlCS;
    wi::graphics::Shader applyVorticityCS;
    wi::graphics::Shader advectTempCS;

    // MacCormack correction shaders
    wi::graphics::Shader macCormackVelCS;
    wi::graphics::Shader macCormackSmokeCS;
    wi::graphics::Shader macCormackTempCS;

    // Ray-march render shader
    wi::graphics::Shader raymarchCS;

    bool gpuReady = false;
	bool resetPending = false;
};
