// cloth_normalizeNormalsCS.hlsl
// Read accumulated normals from CAS-loop float atomic buffer, normalize, write to float4 normals.

#include "cloth_common.hlsli"

[numthreads(CLOTH_THREAD_GROUP_SIZE, 1, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
    uint pNr = DTid.x;
    if (pNr >= cb.numParticles)
        return;

    // Read accumulated normals (float bits stored as uint)
    uint base = 3u * pNr;
    float nx = asfloat(normAccum[base + 0u]);
    float ny = asfloat(normAccum[base + 1u]);
    float nz = asfloat(normAccum[base + 2u]);

    normals[pNr] = float4(normalize(float3(nx, ny, nz)), 0.0f);
}
