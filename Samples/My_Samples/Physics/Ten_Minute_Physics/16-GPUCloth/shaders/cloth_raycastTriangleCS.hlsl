// cloth_raycastTriangleCS.hlsl
// GPU raycast: Möller-Trumbore per-triangle intersection.

#include "cloth_common.hlsli"

cbuffer RaycastCB : register(b1)
{
    RaycastConstants rayCB;
};

[numthreads(CLOTH_THREAD_GROUP_SIZE, 1, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
    uint triNr = DTid.x;
    // numTris passed via constant buffer _pad0 field or we use buffer size
    // We'll use a simple bound check with the triDist buffer
    uint numTris;
    uint stride;
    triDist.GetDimensions(numTris, stride);
    if (triNr >= numTris)
        return;

    const float noHit = 1e9f;
    float3 orig = float3(rayCB.origX, rayCB.origY, rayCB.origZ);
    float3 dir  = float3(rayCB.dirX,  rayCB.dirY,  rayCB.dirZ);

    uint id0 = triIds[3u * triNr];
    uint id1 = triIds[3u * triNr + 1u];
    uint id2 = triIds[3u * triNr + 2u];

    float3 p0 = pos[id0].xyz;
    float3 p1 = pos[id1].xyz;
    float3 p2 = pos[id2].xyz;

    float3 e1 = p1 - p0;
    float3 e2 = p2 - p0;
    float3 pv = cross(dir, e2);
    float det = dot(e1, pv);

    if (abs(det) < 1e-9f)
    {
        triDist[triNr] = noHit;
        return;
    }

    float inv = 1.0f / det;
    float3 tv = orig - p0;
    float u = dot(tv, pv) * inv;
    if (u < 0.0f || u > 1.0f)
    {
        triDist[triNr] = noHit;
        return;
    }

    float3 qv = cross(tv, e1);
    float v = dot(dir, qv) * inv;
    if (v < 0.0f || u + v > 1.0f)
    {
        triDist[triNr] = noHit;
        return;
    }

    float d = dot(e2, qv) * inv;
    if (d <= 0.0f)
    {
        triDist[triNr] = noHit;
        return;
    }
    triDist[triNr] = d;
}
