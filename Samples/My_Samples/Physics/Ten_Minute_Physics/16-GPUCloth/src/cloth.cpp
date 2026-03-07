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

// ═══════════════════════════════════════════════════════════════════
// Constructor — matches HTML ClothSim constructor exactly
// ═══════════════════════════════════════════════════════════════════
ClothMesh::ClothMesh(const SimulationParams& inParams)
    : params(inParams)
{
    numX = params.numX;
    numY = params.numY;

    // Force even (required for correct constraint coloring)
    if (numX % 2 == 1) numX++;
    if (numY % 2 == 1) numY++;

    float spacing = params.spacing;
    float clothY = params.clothY;

    // ── Particles ──
    numParticles = (numX + 1) * (numY + 1);

    cpuPos.resize(numParticles);
    cpuInvMass.resize(numParticles);

    for (int xi = 0; xi <= numX; xi++)
    {
        for (int yi = 0; yi <= numY; yi++)
        {
            int id = xi * (numY + 1) + yi;
            float x = (-numX * 0.5f + xi) * spacing;
            float y = clothY;
            float z = (-numY * 0.5f + yi) * spacing;
            cpuPos[id] = XMFLOAT4(x, y, z, 0.0f);
            cpuInvMass[id] = 1.0f;
        }
    }

    restPos = cpuPos;

    // ── Constraints ──
    BuildConstraintPasses();

    // ── Triangles ──
    numTris = 2 * numX * numY;
    triIds.resize(numTris * 3);

    int ti = 0;
    for (int xi = 0; xi < numX; xi++)
    {
        for (int yi = 0; yi < numY; yi++)
        {
            int id0 = xi * (numY + 1) + yi;
            int id1 = (xi + 1) * (numY + 1) + yi;
            int id2 = (xi + 1) * (numY + 1) + yi + 1;
            int id3 = xi * (numY + 1) + yi + 1;
            triIds[ti++] = id0; triIds[ti++] = id1; triIds[ti++] = id2;
            triIds[ti++] = id0; triIds[ti++] = id2; triIds[ti++] = id3;
        }
    }

    // ── Normals (CPU side for readback) ──
    cpuNormals.resize(numParticles, XMFLOAT3(0, 1, 0));
}

// ═══════════════════════════════════════════════════════════════════
// Constraint generation — matches HTML exactly
// ═══════════════════════════════════════════════════════════════════
void ClothMesh::BuildConstraintPasses()
{
    // Pass sizes matching HTML:
    passSizes[0] = (numX + 1) * (numY / 2);           // stretch Y even
    passSizes[1] = (numX + 1) * (numY / 2);           // stretch Y odd
    passSizes[2] = (numX / 2) * (numY + 1);           // stretch X even
    passSizes[3] = (numX / 2) * (numY + 1);           // stretch X odd
    passSizes[4] = 2 * numX * numY                    // shear
                 + (numX + 1) * (numY - 1)            // bending Y
                 + (numY + 1) * (numX - 1);           // bending X

    numDistConstraints = 0;
    for (int p = 0; p < NUM_PASSES; p++)
        numDistConstraints += passSizes[p];

    constIds.resize(2 * numDistConstraints);

    int ci = 0;

    // Pass 0 & 1: stretch Y (along Y axis, even/odd)
    for (int p = 0; p < 2; p++)
    {
        for (int xi = 0; xi <= numX; xi++)
        {
            for (int yi = 0; yi < numY / 2; yi++)
            {
                constIds[2 * ci]     = xi * (numY + 1) + 2 * yi + p;
                constIds[2 * ci + 1] = xi * (numY + 1) + 2 * yi + p + 1;
                ci++;
            }
        }
    }

    // Pass 2 & 3: stretch X (along X axis, even/odd)
    for (int p = 0; p < 2; p++)
    {
        for (int xi = 0; xi < numX / 2; xi++)
        {
            for (int yi = 0; yi <= numY; yi++)
            {
                constIds[2 * ci]     = (2 * xi + p) * (numY + 1) + yi;
                constIds[2 * ci + 1] = (2 * xi + p + 1) * (numY + 1) + yi;
                ci++;
            }
        }
    }

    // Pass 4: shear
    for (int xi = 0; xi < numX; xi++)
    {
        for (int yi = 0; yi < numY; yi++)
        {
            constIds[2 * ci]     = xi * (numY + 1) + yi;
            constIds[2 * ci + 1] = (xi + 1) * (numY + 1) + yi + 1;
            ci++;
            constIds[2 * ci]     = (xi + 1) * (numY + 1) + yi;
            constIds[2 * ci + 1] = xi * (numY + 1) + yi + 1;
            ci++;
        }
    }

    // Pass 4 continued: bending Y
    for (int xi = 0; xi <= numX; xi++)
    {
        for (int yi = 0; yi < numY - 1; yi++)
        {
            constIds[2 * ci]     = xi * (numY + 1) + yi;
            constIds[2 * ci + 1] = xi * (numY + 1) + yi + 2;
            ci++;
        }
    }

    // Pass 4 continued: bending X
    for (int xi = 0; xi < numX - 1; xi++)
    {
        for (int yi = 0; yi <= numY; yi++)
        {
            constIds[2 * ci]     = xi * (numY + 1) + yi;
            constIds[2 * ci + 1] = (xi + 2) * (numY + 1) + yi;
            ci++;
        }
    }
}

// ═══════════════════════════════════════════════════════════════════
// GPU initialization
// ═══════════════════════════════════════════════════════════════════
void ClothMesh::InitGPUBuffers()
{
    CreateGPUBuffers();
    LoadShaders();

    // Rest lengths will be computed on the main thread's first simulation frame
    // to avoid SubmitCommandLists/WaitForGPU from a background job (causes VK_TIMEOUT).
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

    // Zero-initialized uint buffers (for CAS-loop float atomic accumulators, 3 per particle)
    std::vector<uint32_t> zeroUints(3 * numParticles, 0);

    // Particle buffers (float4)
    makeRWStructured(posBuffer, sizeof(XMFLOAT4), numParticles, cpuPos.data(), "cloth::pos");
    makeRWStructured(prevPosBuffer, sizeof(XMFLOAT4), numParticles, cpuPos.data(), "cloth::prevPos");
    makeRWStructured(velBuffer, sizeof(XMFLOAT4), numParticles, zeroVel.data(), "cloth::vel");
    makeRWStructured(invMassBuffer, sizeof(float), numParticles, cpuInvMass.data(), "cloth::invMass");

    // Constraint buffers
    makeRWStructured(constIdsBuffer, sizeof(uint32_t), numDistConstraints * 2, constIds.data(), "cloth::constIds");
    makeRWStructured(restLengthsBuffer, sizeof(float), numDistConstraints, nullptr, "cloth::restLengths");

    // Jacobi corrections buffer (uint, 3 per particle — stores float bits via CAS atomics)
    makeRWStructured(correctionsBuffer, sizeof(uint32_t), 3 * numParticles, zeroUints.data(), "cloth::corrections");

    // Normals buffer (float4)
    makeRWStructured(normalsBuffer, sizeof(XMFLOAT4), numParticles, nullptr, "cloth::normals");

    // Triangle buffer
    makeRWStructured(triIdsBuffer, sizeof(uint32_t), numTris * 3, triIds.data(), "cloth::triIds");

    // Normal accumulation buffer (uint, 3 per particle — stores float bits via CAS atomics)
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
    const std::string originalBinPath = wi::renderer::GetShaderPath();
    const std::string originalSrcPath = wi::renderer::GetShaderSourcePath();

    // .cso outputs go to the binary directory
    std::string shaderBinPath = wi::helper::GetCurrentPath() + "/shaders/";
    wi::renderer::SetShaderPath(shaderBinPath);

    // .hlsl sources are read directly from the source tree (no copy needed).
    // SHADER_SOURCE_DIR is set by CMake via target_compile_definitions.
#ifdef SHADER_SOURCE_DIR
    wi::renderer::SetShaderSourcePath(SHADER_SOURCE_DIR);
#else
    wi::renderer::SetShaderSourcePath(shaderBinPath);
#endif

    wi::renderer::LoadShader(ShaderStage::CS, computeRestLengthsCS, "cloth_computeRestLengthsCS.cso");
    wi::renderer::LoadShader(ShaderStage::CS, integrateCS, "cloth_integrateCS.cso");
    wi::renderer::LoadShader(ShaderStage::CS, solveConstraintsCS, "cloth_solveConstraintsCS.cso");
    wi::renderer::LoadShader(ShaderStage::CS, addCorrectionsCS, "cloth_addCorrectionsCS.cso");
    wi::renderer::LoadShader(ShaderStage::CS, updateVelCS, "cloth_updateVelCS.cso");
    wi::renderer::LoadShader(ShaderStage::CS, clearNormalsCS, "cloth_clearNormalsCS.cso");
    wi::renderer::LoadShader(ShaderStage::CS, addNormalsCS, "cloth_addNormalsCS.cso");
    wi::renderer::LoadShader(ShaderStage::CS, normalizeNormalsCS, "cloth_normalizeNormalsCS.cso");
    wi::renderer::LoadShader(ShaderStage::CS, raycastTriangleCS, "cloth_raycastTriangleCS.cso");

    wi::renderer::SetShaderPath(originalBinPath);
    wi::renderer::SetShaderSourcePath(originalSrcPath);
}

void ClothMesh::ComputeRestLengthsGPU(wi::graphics::CommandList cmd)
{
    GraphicsDevice* device = wi::graphics::GetDevice();

    // Use numConstraintsInPass to pass total constraint count
    ClothSimConstants cb = {};
    cb.numConstraintsInPass = numDistConstraints;

    device->BindComputeShader(&computeRestLengthsCS, cmd);
    device->BindDynamicConstantBuffer(cb, CLOTH_CB_SLOT, cmd);

    const GPUResource* uav_pos[] = { &posBuffer };
    device->BindUAVs(uav_pos, CLOTH_SLOT_POS, 1, cmd);
    const GPUResource* uav_const[] = { &constIdsBuffer, &restLengthsBuffer };
    device->BindUAVs(uav_const, CLOTH_SLOT_CONST_IDS, 2, cmd);

    uint32_t dispatchCount = (numDistConstraints + CLOTH_THREAD_GROUP_SIZE - 1) / CLOTH_THREAD_GROUP_SIZE;
    device->Dispatch(dispatchCount, 1, 1, cmd);

    GPUBarrier barrier = GPUBarrier::Memory();
    device->Barrier(&barrier, 1, cmd);
}

// ═══════════════════════════════════════════════════════════════════
// SimulateGPU — matches HTML simulate() exactly
// ═══════════════════════════════════════════════════════════════════
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

    GPUBarrier barrier = GPUBarrier::Memory();

    // Helper lambdas for binding UAV groups (avoids nullptr in arrays)
    auto bindIntegrateUAVs = [&]() {
        const GPUResource* uavs[] = { &posBuffer, &prevPosBuffer, &velBuffer, &invMassBuffer };
        device->BindUAVs(uavs, CLOTH_SLOT_POS, 4, cmd);
    };

    auto bindSolveConstraintUAVs = [&]() {
        const GPUResource* uav_pos[] = { &posBuffer };
        device->BindUAVs(uav_pos, CLOTH_SLOT_POS, 1, cmd);
        const GPUResource* uav_mass[] = { &invMassBuffer, &constIdsBuffer, &restLengthsBuffer,
                                           &correctionsBuffer };
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

    for (int step = 0; step < numSubSteps; step++)
    {
        // ===== 1. Integrate =====
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

            device->BindComputeShader(&integrateCS, cmd);
            device->BindDynamicConstantBuffer(cb, CLOTH_CB_SLOT, cmd);
            bindIntegrateUAVs();
            device->Dispatch(dispatchParticles, 1, 1, cmd);
            device->Barrier(&barrier, 1, cmd);
        }

        // ===== 2. Solve constraints =====
        if (solveType == 0)
        {
            // Coloring hybrid — each pass dispatched separately
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
            // Full Jacobi — all constraints at once
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

// ═══════════════════════════════════════════════════════════════════
// UpdateMeshGPU — matches HTML updateNormals() exactly
// ═══════════════════════════════════════════════════════════════════
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

// ═══════════════════════════════════════════════════════════════════
// Reset — matches HTML reset()
// ═══════════════════════════════════════════════════════════════════
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

// ═══════════════════════════════════════════════════════════════════
// Grab — matches HTML startDrag/drag/endDrag
// ═══════════════════════════════════════════════════════════════════
bool ClothMesh::StartGrabGPU(const XMFLOAT3& rayOrigin, const XMFLOAT3& rayDir)
{
    if (!gpuBuffersReady)
        return false;

    GraphicsDevice* device = wi::graphics::GetDevice();
    CommandList cmd = device->BeginCommandList();

    // Upload ray parameters
    RaycastConstants rayCB = {};
    rayCB.origX = rayOrigin.x; rayCB.origY = rayOrigin.y; rayCB.origZ = rayOrigin.z;
    rayCB.dirX = rayDir.x; rayCB.dirY = rayDir.y; rayCB.dirZ = rayDir.z;

    device->BindComputeShader(&raycastTriangleCS, cmd);
    device->BindDynamicConstantBuffer(rayCB, CLOTH_CB_RAYCAST_SLOT, cmd);

    const GPUResource* uav_pos[] = { &posBuffer };
    device->BindUAVs(uav_pos, CLOTH_SLOT_POS, 1, cmd);
    const GPUResource* uav_tri[] = { &triIdsBuffer };
    device->BindUAVs(uav_tri, CLOTH_SLOT_TRI_IDS, 1, cmd);
    const GPUResource* uav_dist[] = { &triDistBuffer };
    device->BindUAVs(uav_dist, CLOTH_SLOT_TRI_DIST, 1, cmd);
    device->Dispatch(dispatchTris, 1, 1, cmd);

    GPUBarrier barrier = GPUBarrier::Memory();
    device->Barrier(&barrier, 1, cmd);

    // Copy to readback
    GPUBarrier copyBarrier = GPUBarrier::Buffer(&triDistBuffer,
        ResourceState::UNORDERED_ACCESS, ResourceState::COPY_SRC);
    device->Barrier(&copyBarrier, 1, cmd);
    device->CopyResource(&triDistReadbackBuffer, &triDistBuffer, cmd);
    copyBarrier = GPUBarrier::Buffer(&triDistBuffer,
        ResourceState::COPY_SRC, ResourceState::UNORDERED_ACCESS);
    device->Barrier(&copyBarrier, 1, cmd);

    device->SubmitCommandLists();
    device->WaitForGPU();

    // Read results
    if (!triDistReadbackBuffer.mapped_data)
        return false;

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

    if (minTri < 0 || minDist >= 1e8f)
    {
        dragParticleNr = -1;
        return false;
    }

    // Get first vertex of hit triangle
    dragParticleNr = triIds[3 * minTri];
    dragDepth = minDist;
    dragInvMass = cpuInvMass[dragParticleNr];

    // Store drag position (integrate shader will set it via CB)
    currentDragPos.x = rayOrigin.x + dragDepth * rayDir.x;
    currentDragPos.y = rayOrigin.y + dragDepth * rayDir.y;
    currentDragPos.z = rayOrigin.z + dragDepth * rayDir.z;

    // invMass pinning is handled in shaders via CB (dragParticleNr check)
    return true;
}

void ClothMesh::DragGPU(const XMFLOAT3& rayOrigin, const XMFLOAT3& rayDir, wi::graphics::CommandList cmd)
{
    if (dragParticleNr < 0 || !gpuBuffersReady)
        return;

    // Store drag position — will be passed to integrate shader via constant buffer
    currentDragPos.x = rayOrigin.x + dragDepth * rayDir.x;
    currentDragPos.y = rayOrigin.y + dragDepth * rayDir.y;
    currentDragPos.z = rayOrigin.z + dragDepth * rayDir.z;
}

void ClothMesh::EndGrabGPU(wi::graphics::CommandList cmd)
{
    if (dragParticleNr < 0 || !gpuBuffersReady)
        return;

    // invMass is never modified on GPU — shaders use CB dragParticleNr check.
    // Just reset the drag state; next frame shaders see dragParticleNr = -1.
    dragParticleNr = -1;
}

// ═══════════════════════════════════════════════════════════════════
// Readback
// ═══════════════════════════════════════════════════════════════════
void ClothMesh::RequestPositionsReadback(wi::graphics::CommandList cmd)
{
    if (!gpuBuffersReady)
        return;

    GraphicsDevice* device = wi::graphics::GetDevice();

    GPUBarrier barrier = GPUBarrier::Buffer(&posBuffer,
        ResourceState::UNORDERED_ACCESS, ResourceState::COPY_SRC);
    device->Barrier(&barrier, 1, cmd);

    device->CopyResource(&posReadbackBuffer, &posBuffer, cmd);

    barrier = GPUBarrier::Buffer(&posBuffer,
        ResourceState::COPY_SRC, ResourceState::UNORDERED_ACCESS);
    device->Barrier(&barrier, 1, cmd);

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
    if (!gpuBuffersReady || !posReadbackPending)
        return;

    if (posReadbackBuffer.mapped_data)
    {
        const XMFLOAT4* src = static_cast<const XMFLOAT4*>(posReadbackBuffer.mapped_data);
        for (int i = 0; i < numParticles; i++)
            cpuPos[i] = src[i];
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
