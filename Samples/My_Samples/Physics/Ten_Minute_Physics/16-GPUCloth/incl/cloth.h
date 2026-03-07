#pragma once

#include <memory>
#include <vector>
#include <cstdint>
#include <wiScene.h>
#include <wiGraphics.h>
#include "shader_interop_cloth.h"

// Simulation parameters
struct SimulationParams
{
    int numX = 500;
    int numZ = 500;
    float spacing = 0.01f;
    float clothY = 2.2f;

    // Sphere collision
    float sphereCenter[3] = {0.0f, 1.5f, 0.0f};
    float sphereRadius = 0.5f;
};

class ClothMesh
{
public:
    ClothMesh(const SimulationParams& params);

    // Initialize GPU buffers and shaders (call once after construction)
    void InitGPUBuffers();

    // GPU simulation (called every frame)
    void SimulateGPU(float dt, wi::graphics::CommandList cmd, int solveType);

    // Update mesh normals on GPU
    void UpdateMeshNormalsGPU(wi::graphics::CommandList cmd);

    // GPU raycast for grabbing (returns true if hit found)
    bool StartGrabGPU(const XMFLOAT3& rayOrigin, const XMFLOAT3& rayDir);
    void DragGPU(const XMFLOAT3& rayOrigin, const XMFLOAT3& rayDir, wi::graphics::CommandList cmd);
    void EndGrabGPU(wi::graphics::CommandList cmd);

    // Reset to initial state
    void Reset();

    // Readback management
    void RequestPositionsReadback(wi::graphics::CommandList cmd);
    void RequestNormalsReadback(wi::graphics::CommandList cmd);
    void ProcessPositionsReadback();
    void ProcessNormalsReadback();

    // =====================================================================
    // Simulation parameters
    // =====================================================================
    SimulationParams params;

    int numSubSteps = 30;
    float dt = 1.0f / 60.0f;
    float jacobiScale = 0.15f;
    XMFLOAT3 gravity = {0.0f, -10.0f, 0.0f};

    // =====================================================================
    // CPU particle data
    // =====================================================================
    int numParticles = 0;
    int numX = 0;
    int numZ = 0;

    std::vector<XMFLOAT4> cpuPos;        // Current positions
    std::vector<float> cpuInvMass;       // Inverse mass
    std::vector<XMFLOAT4> restPos;       // Rest positions
    std::vector<XMFLOAT3> cpuNormals;    // Normals

    // =====================================================================
    // Constraint data — 5 passes
    // =====================================================================
    static constexpr int NUM_PASSES = 5;
    int passSizes[NUM_PASSES] = {};
    bool passIndependent[NUM_PASSES] = {true, true, true, true, false};

    int numDistConstraints = 0;
    std::vector<uint32_t> constIds;

    // =====================================================================
    // Triangle data
    // =====================================================================
    std::vector<uint32_t> triIds;
    int numTris = 0;

    // =====================================================================
    // GPU buffers
    // =====================================================================
    wi::graphics::GPUBuffer posBuffer;          // float4
    wi::graphics::GPUBuffer prevPosBuffer;      // float4
    wi::graphics::GPUBuffer velBuffer;          // float4
    wi::graphics::GPUBuffer invMassBuffer;      // float
    wi::graphics::GPUBuffer constIdsBuffer;     // uint
    wi::graphics::GPUBuffer restLengthsBuffer;  // float
    wi::graphics::GPUBuffer correctionsBuffer;  // uint (3*N, stores float bits via CAS atomics)
    wi::graphics::GPUBuffer normalsBuffer;      // float4
    wi::graphics::GPUBuffer triIdsBuffer;       // uint
    wi::graphics::GPUBuffer normAccumBuffer;    // uint (3*N, stores float bits via CAS atomics)
    wi::graphics::GPUBuffer triDistBuffer;      // float (raycast output)

    // Readback buffers
    wi::graphics::GPUBuffer posReadbackBuffer;
    wi::graphics::GPUBuffer normalsReadbackBuffer;
    wi::graphics::GPUBuffer triDistReadbackBuffer;

    bool posReadbackPending = false;
    bool normalsReadbackPending = false;

    // =====================================================================
    // Compute shaders
    // =====================================================================
    wi::graphics::Shader computeRestLengthsCS;
    wi::graphics::Shader integrateCS;
    wi::graphics::Shader solveConstraintsCS;
    wi::graphics::Shader addCorrectionsCS;
    wi::graphics::Shader updateVelCS;
    wi::graphics::Shader clearNormalsCS;
    wi::graphics::Shader addNormalsCS;
    wi::graphics::Shader normalizeNormalsCS;
    wi::graphics::Shader raycastTriangleCS;

    bool gpuBuffersReady = false;
    bool restLengthsComputed = false;

    // =====================================================================
    // Dispatch sizes
    // =====================================================================
    uint32_t dispatchParticles = 0;
    uint32_t dispatchTris = 0;

    // =====================================================================
    // Grabbing state
    // =====================================================================
    int dragParticleNr = -1;
    float dragDepth = 0.0f;
    XMFLOAT3 currentDragPos = {0, 0, 0};

private:
    void BuildConstraintPasses();
    void CreateGPUBuffers();
    void LoadShaders();
    void ComputeRestLengthsGPU(wi::graphics::CommandList cmd);
};

// SimulationObject: cloth with associated wireframe mesh entity
struct SimulationObject
{
    std::unique_ptr<ClothMesh> cloth;
    uint64_t wireEntity = 0;
};
