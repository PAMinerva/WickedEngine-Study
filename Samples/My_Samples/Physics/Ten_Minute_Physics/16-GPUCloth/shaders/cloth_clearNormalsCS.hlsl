// cloth_clearNormalsCS.hlsl
// Zero the normAccum buffer before normal accumulation.

#include "cloth_common.hlsli"

[numthreads(CLOTH_THREAD_GROUP_SIZE, 1, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
    uint pNr = DTid.x;
    if (pNr >= cb.numParticles)
        return;

    uint base = 3u * pNr;
    normAccum[base + 0u] = 0;
    normAccum[base + 1u] = 0;
    normAccum[base + 2u] = 0;
}
