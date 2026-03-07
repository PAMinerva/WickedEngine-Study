// cloth_integrateCS.hlsl
// Integration: save prevPos, apply gravity, integrate, sphere + ground collision.
// Matches HTML sh-integrate exactly.

#include "cloth_common.hlsli"

[numthreads(CLOTH_THREAD_GROUP_SIZE, 1, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
    uint pNr = DTid.x;
    if (pNr >= cb.numParticles)
        return;

    prevPos[pNr] = pos[pNr]; // test

    // Dragged particle: set position from constant buffer, skip integration
    if ((int)pNr == cb.dragParticleNr)
    {
        float4 dp = float4(cb.dragPosX, cb.dragPosY, cb.dragPosZ, 0.0f);
        pos[pNr] = dp;
        prevPos[pNr] = dp;
        return;
    }

    float w = invMass[pNr];
    if (w == 0.0f)
        return;

    float3 grav = float3(cb.gravX, cb.gravY, cb.gravZ);
    float dt = cb.dt;
    float3 v = vel[pNr].xyz;
    float3 p = pos[pNr].xyz;

    v = v + grav * dt;
    p = p + v * dt;

    // Sphere collision
    float3 sc = float3(cb.sphereCX, cb.sphereCY, cb.sphereCZ);
    float sr = cb.sphereR;
    float thickness = 0.001f;
    float friction = 0.01f;

    float d = length(p - sc);
    if (d < (sr + thickness))
    {
        float3 pp = p * (1.0f - friction) + prevPos[pNr].xyz * friction;
        float3 r = pp - sc;
        d = length(r);
        p = sc + r * ((sr + thickness) / d);
    }

    // Ground collision
    if (p.y < thickness)
    {
        float3 pp = p * (1.0f - friction) + prevPos[pNr].xyz * friction;
        p = float3(pp.x, thickness, pp.z);
    }

    pos[pNr] = float4(p, 0.0f);
    vel[pNr] = float4(v, 0.0f);
}
