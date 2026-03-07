// cloth_computeRestLengthsCS.hlsl
// One-time initialization: compute rest lengths from initial positions.
// One thread per constraint.

#include "cloth_common.hlsli"

[numthreads(CLOTH_THREAD_GROUP_SIZE, 1, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
    uint cNr = DTid.x;
    uint n = cb.numConstraintsInPass;  // total constraints passed via this field
    if (cNr >= n)
        return;

    uint id0 = constIds[2u * cNr];
    uint id1 = constIds[2u * cNr + 1u];
    float3 p0 = pos[id0].xyz;
    float3 p1 = pos[id1].xyz;
    restLengths[cNr] = length(p1 - p0);
}
