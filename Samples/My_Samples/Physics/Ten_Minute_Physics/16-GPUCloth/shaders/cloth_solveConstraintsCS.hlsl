// cloth_solveConstraintsCS.hlsl
// Unified constraint solver for both coloring and Jacobi modes.
// When cb.jacobiScale == 0: direct position write (coloring, no race by construction).
// When cb.jacobiScale >  0: CAS-loop float atomic accumulation (Jacobi).

#include "cloth_common.hlsli"

[numthreads(CLOTH_THREAD_GROUP_SIZE, 1, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
    uint localIdx = DTid.x;
    if (localIdx >= cb.numConstraintsInPass)
        return;

    uint cNr = cb.firstConstraint + localIdx;

    uint id0 = constIds[2u * cNr];
    uint id1 = constIds[2u * cNr + 1u];
    float w0 = ((int)id0 == cb.dragParticleNr) ? 0.0f : invMass[id0];
    float w1 = ((int)id1 == cb.dragParticleNr) ? 0.0f : invMass[id1];
    float w = w0 + w1;
    if (w == 0.0f)
        return;

    float3 p0 = pos[id0].xyz;
    float3 p1 = pos[id1].xyz;
    float3 dv = p1 - p0;
    float l = length(dv);
    float l0 = restLengths[cNr];
    if (l < 1e-9f)
        return;
    float3 n = dv / l;
    float3 dP = n * (l - l0) / w;

    float3 c0 =  w0 * dP;
    float3 c1 = -w1 * dP;

    if (cb.jacobiScale > 0.0f)
    {
        // Jacobi: CAS-loop float atomic accumulation
        ATOMIC_ADD_FLOAT(corrections, 3u * id0 + 0u, c0.x);
        ATOMIC_ADD_FLOAT(corrections, 3u * id0 + 1u, c0.y);
        ATOMIC_ADD_FLOAT(corrections, 3u * id0 + 2u, c0.z);
        ATOMIC_ADD_FLOAT(corrections, 3u * id1 + 0u, c1.x);
        ATOMIC_ADD_FLOAT(corrections, 3u * id1 + 1u, c1.y);
        ATOMIC_ADD_FLOAT(corrections, 3u * id1 + 2u, c1.z);
    }
    else
    {
        // Coloring: direct position write (no race by construction)
        pos[id0] = float4(p0 + c0, 0.0f);
        pos[id1] = float4(p1 + c1, 0.0f);
    }
}
