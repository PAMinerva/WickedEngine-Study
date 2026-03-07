// cloth_addNormalsCS.hlsl
// Accumulate face normals per vertex using CAS-loop float atomic addition.
// One thread per triangle.

#include "cloth_common.hlsli"

[numthreads(CLOTH_THREAD_GROUP_SIZE, 1, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
    uint triNr = DTid.x;
    // numConstraintsInPass is repurposed as numTris for this dispatch
    if (triNr >= cb.numConstraintsInPass)
        return;

    uint id0 = triIds[3u * triNr];
    uint id1 = triIds[3u * triNr + 1u];
    uint id2 = triIds[3u * triNr + 2u];

    float3 e1 = pos[id1].xyz - pos[id0].xyz;
    float3 e2 = pos[id2].xyz - pos[id0].xyz;
    float3 n = cross(e2, e1);

    // CAS-loop float atomic accumulation for normals
    ATOMIC_ADD_FLOAT(normAccum, 3u * id0 + 0u, n.x);
    ATOMIC_ADD_FLOAT(normAccum, 3u * id0 + 1u, n.y);
    ATOMIC_ADD_FLOAT(normAccum, 3u * id0 + 2u, n.z);
    ATOMIC_ADD_FLOAT(normAccum, 3u * id1 + 0u, n.x);
    ATOMIC_ADD_FLOAT(normAccum, 3u * id1 + 1u, n.y);
    ATOMIC_ADD_FLOAT(normAccum, 3u * id1 + 2u, n.z);
    ATOMIC_ADD_FLOAT(normAccum, 3u * id2 + 0u, n.x);
    ATOMIC_ADD_FLOAT(normAccum, 3u * id2 + 1u, n.y);
    ATOMIC_ADD_FLOAT(normAccum, 3u * id2 + 2u, n.z);
}
