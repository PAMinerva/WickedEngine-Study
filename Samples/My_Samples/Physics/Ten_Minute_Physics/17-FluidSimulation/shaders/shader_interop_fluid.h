#pragma once

// ============================================================================
// shader_interop_fluid.h
// Shared between C++ and HLSL — keep types compatible with both languages.
// ============================================================================

#ifdef __cplusplus
#include <cstdint>
#endif

// ---------------------------------------------------------------------------
// Thread group sizes
// ---------------------------------------------------------------------------
#define FLUID_THREADS_3D_X  8
#define FLUID_THREADS_3D_Y  8
#define FLUID_THREADS_3D_Z  8

// For 2D fullscreen ray-march dispatch
#define FLUID_THREADS_2D    16

// ---------------------------------------------------------------------------
// UAV/SRV slot assignments for simulation kernels
// ---------------------------------------------------------------------------
#define FLUID_SLOT_VELOCITY_U       0
#define FLUID_SLOT_VELOCITY_V       1
#define FLUID_SLOT_VELOCITY_W       2
#define FLUID_SLOT_PRESSURE         3
#define FLUID_SLOT_DIVERGENCE       4
#define FLUID_SLOT_SMOKE            5
#define FLUID_SLOT_SOLID            6

// Temp buffers for ping-pong advection (written as UAV)
#define FLUID_SLOT_VELOCITY_U_TEMP  7
#define FLUID_SLOT_VELOCITY_V_TEMP  8
#define FLUID_SLOT_VELOCITY_W_TEMP  9
#define FLUID_SLOT_SMOKE_TEMP       10

// Temperature field (cell-centered scalar)
#define FLUID_SLOT_TEMPERATURE      11
#define FLUID_SLOT_TEMPERATURE_TEMP 12

// Curl / vorticity (cell-centered float4: xyz = curl, w unused)
#define FLUID_SLOT_CURL             13

// ---------------------------------------------------------------------------
// SRV slot assignments for ray-march render shader
// ---------------------------------------------------------------------------
#define FLUID_SRV_SMOKE             0
#define FLUID_SRV_PRESSURE          1
#define FLUID_SRV_VELOCITY_U        2
#define FLUID_SRV_VELOCITY_V        3
#define FLUID_SRV_VELOCITY_W        4

// UAV for render output
#define FLUID_UAV_RENDER_OUTPUT     0

// ---------------------------------------------------------------------------
// Constant buffer slot
// ---------------------------------------------------------------------------
#define FLUID_CB_SLOT               0

// ---------------------------------------------------------------------------
// Simulation constants — pushed every frame via dynamic constant buffer
// ---------------------------------------------------------------------------
struct FluidConstants
{
    // Grid dimensions
    int gridX;
    int gridY;
    int gridZ;
    float cellSize;

    float dt;
    float density;
    float overRelaxation;
    int   numPressureIters;

    // Gravity
    float gravX;
    float gravY;
    float gravZ;
    int   redBlackPass;

    // Obstacle (solid sphere)
    float obsCenterX;
    float obsCenterY;
    float obsCenterZ;
    float obsRadius;

    float obsVelX;
    float obsVelY;
    float obsVelZ;
    float buoyancy;             // upward force per unit temperature difference

    // Smoke source (draggable sphere emitter)
    float smokeSrcX;
    float smokeSrcY;
    float smokeSrcZ;
    float smokeSrcRadius;

    float smokeInflowDensity;
    float smokeDissipation;     // smoke fade rate per second (0 = no fade)
    int   enableSmoke;
    float vorticityStrength;    // epsilon for vorticity confinement (0.35 typical)

    // Temperature-based buoyancy
    float tempInjection;        // temperature value injected at source (12.0 typical)
    float ambientTemp;          // reference temperature (0.0 typical)
    float fluidWeight;          // weight pulling dense smoke down (0.125 typical)
    float velDissipation;       // velocity kept per step (0.9995 typical, 1.0 = no loss)
};

// ---------------------------------------------------------------------------
// Ray-march render constants
// ---------------------------------------------------------------------------
struct FluidRenderConstants
{
    // Grid info
    int gridX;
    int gridY;
    int gridZ;
    float cellSize;

    // Domain AABB in world space
    float domainMinX;
    float domainMinY;
    float domainMinZ;
    float domainMaxX;

    float domainMaxY;
    float domainMaxZ;
    int   renderMode;
    float maxSteps;

    // Camera inverse VP for ray reconstruction
#ifdef __cplusplus
    float invVP[16];
#else
    row_major float4x4 invVP;
#endif

    // Screen resolution
    int screenWidth;
    int screenHeight;
    float smokeAbsorption;
    float _pad0;

    // Obstacle sphere (for solid rendering in ray-march)
    float obsCenterX;
    float obsCenterY;
    float obsCenterZ;
    float obsRadius;

    // Camera world position
    float cameraPosX;
    float cameraPosY;
    float cameraPosZ;

    // Smoke source sphere (rendered as faint glow in ray-march)
    float smokeSrcRenderX;
    float smokeSrcRenderY;
    float smokeSrcRenderZ;
    float smokeSrcRenderRadius;

    // Configurable smoke color (default: warm white-gray)
    float smokeColorR;
    float smokeColorG;
    float smokeColorB;
    float _pad1;
};
