// cloth_addCorrectionsCS.hlsl
// Apply accumulated Jacobi corrections: read from corrections buffer (float bits stored as uint),
// apply to positions with jacobiScale, then zero the buffer for next pass.

#include "cloth_common.hlsli"

[numthreads(CLOTH_THREAD_GROUP_SIZE, 1, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
    uint pNr = DTid.x;
    if (pNr >= cb.numParticles)
        return;

    // Read accumulated corrections (stored as float bits in uint buffer)
    uint base = 3u * pNr;
    float cx = asfloat(corrections[base + 0u]);
    float cy = asfloat(corrections[base + 1u]);
    float cz = asfloat(corrections[base + 2u]);

    // Zero for next pass (safe: one thread per particle, barrier before next Jacobi dispatch)
    corrections[base + 0u] = 0;
    corrections[base + 1u] = 0;
    corrections[base + 2u] = 0;

    float scale = cb.jacobiScale;
    float3 correction = float3(cx, cy, cz) * scale;

    // Clamp total applied correction to prevent Jacobi divergence during grab
    float corrLen = length(correction);
    const float maxCorr = 0.1f;
    if (corrLen > maxCorr)
        correction = correction * (maxCorr / corrLen);

    float3 p = pos[pNr].xyz;
    pos[pNr] = float4(p + correction, 0.0f);
}
