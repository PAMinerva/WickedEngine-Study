// cloth_updateStreamoutCS.hlsl
// Copies simulation positions and normals into renderer-compatible buffers.
// Positions: float4 → float4 (direct copy)
// Normals:   float4 → R8G8B8A8_SNORM packed as uint (matches Vertex_NOR encoding)

#include "shader_interop_cloth.h"

cbuffer UpdateStreamoutCB : register(b0)
{
    uint numParticles;
    uint _pad0;
    uint _pad1;
    uint _pad2;
};

// Input: simulation buffers (structured float4)
RWStructuredBuffer<float4> posIn  : register(u0);
RWStructuredBuffer<float4> norIn  : register(u1);

// Output: renderer buffers
RWStructuredBuffer<float4> posOut : register(u2);   // typed SRV as R32G32B32A32_FLOAT during rendering
RWStructuredBuffer<uint>   norOut : register(u3);   // typed SRV as R8G8B8A8_SNORM during rendering

[numthreads(CLOTH_THREAD_GROUP_SIZE, 1, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
    uint pNr = DTid.x;
    if (pNr >= numParticles)
        return;

    // Position: direct copy (float4 → float4)
    posOut[pNr] = posIn[pNr];

    // Normal: pack float3 -> R8G8B8A8_SNORM as uint
    // Matches Wicked Engine Vertex_NOR::FromFULL encoding (multiply by 127.5, truncate)
    float3 n = normalize(norIn[pNr].xyz);
    int snx = (int)(n.x * 127.5f);
    int sny = (int)(n.y * 127.5f);
    int snz = (int)(n.z * 127.5f);

	// Move the 3 signed bytes into a single uint (R8G8B8A8_SNORM format, alpha unused)
    norOut[pNr] = (uint(snx & 0xFF))
               | (uint(sny & 0xFF) << 8)
               | (uint(snz & 0xFF) << 16);
}
