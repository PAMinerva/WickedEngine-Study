// cloth_updateVelCS.hlsl
// Update velocity from position change: vel = (pos - prevPos) / dt
// Matches HTML sh-updatevel exactly.

#include "cloth_common.hlsli"

[numthreads(CLOTH_THREAD_GROUP_SIZE, 1, 1)]
void main(uint3 DTid : SV_DispatchThreadID)
{
    uint pNr = DTid.x;
    if (pNr >= cb.numParticles)
        return;

    float invDt = 1.0f / cb.dt;

    // Dragged particle: force zero velocity
    if ((int)pNr == cb.dragParticleNr)
    {
        vel[pNr] = float4(0.0f, 0.0f, 0.0f, 0.0f);
        return;
    }

    vel[pNr] = float4((pos[pNr].xyz - prevPos[pNr].xyz) * invDt, 0.0f);
}
