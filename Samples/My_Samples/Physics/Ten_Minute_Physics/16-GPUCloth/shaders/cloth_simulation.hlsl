// cloth_simulation.hlsl
// All cloth simulation compute kernels in a single file.
// Each kernel is a named entry point. A unique permutation define per kernel
// ensures LoadShader produces distinct .cso files from the same source.

#include "cloth_common.hlsli"

// ---------------------------------------------------------------------------
// computeRestLengthsCS - One-time init: compute rest lengths from initial positions.
// ---------------------------------------------------------------------------
[numthreads(CLOTH_THREAD_GROUP_SIZE, 1, 1)]
void cloth_computeRestLengthsCS(uint3 DTid : SV_DispatchThreadID)
{
	// Each GPU thread processes one constraint in parallel, identified by the global thread index DTid.x.
	// This shader is dispatched as a 1D grid. See ClothMesh::ComputeRestLengthsGPU() for dispatch size calculation.
    uint cNr = DTid.x;
    uint n = cb.numConstraintsInPass;
    if (cNr >= n)
        return;

	// cNr is the global constraint number. Each constraint connects two particles, whose indices are stored in constIds buffer.
    uint id0 = constIds[2u * cNr];
    uint id1 = constIds[2u * cNr + 1u];
    float3 p0 = pos[id0].xyz;
    float3 p1 = pos[id1].xyz;
    restLengths[cNr] = length(p1 - p0);
}

// ---------------------------------------------------------------------------
// integrateCS - Save prevPos, apply gravity, integrate, sphere + ground collision.
// ---------------------------------------------------------------------------
[numthreads(CLOTH_THREAD_GROUP_SIZE, 1, 1)]
void cloth_integrateCS(uint3 DTid : SV_DispatchThreadID)
{
	// Each GPU thread processes one particle in parallel, identified by global thread index DTid.x.
	// See ClothMesh::SimulateGPU() for dispatch size calculation.
    uint pNr = DTid.x;
    if (pNr >= cb.numParticles)
        return;

	// Save the current position before integrating.
	// This is used later for:
	// - velocity update: vel = (pos - prevPos) / dt
	// - simple collision "friction" trick that blends current/previous positions (see below)
    prevPos[pNr] = pos[pNr];

	// If this particle is being dragged by the mouse,
	// override its position with the drag position from the constant buffer.
    if ((int)pNr == cb.dragParticleNr)
    {
        float4 dp = float4(cb.dragPosX, cb.dragPosY, cb.dragPosZ, 0.0f);
        pos[pNr] = dp;
        prevPos[pNr] = dp;
        return;
    }

	// If invMass == 0 the particle is pinned:
	// - invMass = 1/mass
	// - invMass = 0 means infinite mass -> it does not move under forces or constraints
	// This is used for fixed corners or for the particle being dragged by the mouse.
    float w = invMass[pNr];
    if (w == 0.0f)
        return;

	// ---------------------------------------------------------------------
	// SEMI-IMPLICIT EULER INTEGRATION
	//
	// Update velocity from acceleration (gravity):
	//   v_{t+dt} = v_t + g * dt
	//
	// Then update position using the new velocity:
	//   x_{t+dt} = x_t + v_{t+dt} * dt
	//
    float3 grav = float3(cb.gravX, cb.gravY, cb.gravZ);
    float dt = cb.dt;
    float3 v = vel[pNr].xyz;
    float3 p = pos[pNr].xyz;

    v = v + grav * dt;
    p = p + v * dt;

	// ---------------------------------------------------------------------
	// COLLISION PARAMETERS AND COLLIDER DEFINITIONS
	//
	// This demo resolves collisions directly by projecting particles out of
	// penetration (position correction). Velocities are later recomputed from
	// position differences, so this "position push-out" implicitly affects velocity.
	//
	// Two colliders:
	// 1) A sphere at cs=(sphereCX, sphereCY, sphereCZ) with radius sr=sphereR, defined in the constant buffer.
	// 2) A ground plane at y = 0
	//
    float3 sc = float3(cb.sphereCX, cb.sphereCY, cb.sphereCZ);
    float sr = cb.sphereR;

	// thickness:
	// - expands the colliders a bit so cloth doesn't numerically "buzz" at exact contact
	// - also acts like a minimum separation
	// friction:
	// In this context friction is NOT intended as a physically accurate surface friction model.
	// Instead, it is a simple "friction-like" damping factor applied during contact resolution to
	// reduce sliding and improve visual stability:
	//
	//  During contact, we LERP between the newly integrated position (pos)
	//  and the previous position (prevPos) before projecting out to reduce the effective velocity at contacts.
	//  This can help mitigate sliding and improve visual stability.
	//
	//   bp = lerp(pos, prevPos, friction)
	//     = pos*(1-friction) + prevPos*friction
	//
	// If friction = 0:
	//   bp == pos (no damping; collision uses the fully integrated current position)
	// If friction = 1:
	//   bp == prevPos (max damping; tends to "stick" to the previous position)
	//
	// With friction = 0.01:
	//   bp moves only 1% toward prevPos: a very small damping effect each frame.
	//
	// Why this is called "friction-like":
	// - True surface friction primarily reduces tangential motion (sliding) while permitting
	//   normal separation. This hack does not explicitly decompose into normal/tangent.
	// - However, because we keep re-projecting onto the collider surface, the normal component
	//   is mostly handled by the projection step, and the remaining visually dominant motion at
	//   the contact is often tangential sliding along the surface.
	// - Pulling the contact point slightly toward prevPos tends to reduce the net displacement,
	//   and thus reduces the reconstructed velocity (pos_final - prevPos)/dt, which often appears
	//   as reduced sliding (tangential damping). More on this in the detailed comment below.
    float thickness = 0.001f;
    float friction = 0.01f;

	// =====================================================================
	// SPHERE COLLISION
	// =====================================================================
	// Distance from particle p to sphere center sc:
    float d = length(p - sc);

	// If inside the sphere (or within a small thickness margin), resolve penetration.
    if (d < (sr + thickness))
    {
		// -----------------------------------------------------------------
		// STEP 1: BLEND POSITION (CONTACT DAMPING)
		//
		// bp is ONLY used to compute the contact direction for projection.
		// Think of this as: "during contact, don't fully trust the newly integrated
		// position p; bias it slightly toward where we were last frame."
		//
		// This reduces abrupt changes in the contact direction and reduces the net
		// displacement after projection, which in turn reduces the final velocity
		// reconstructed later.
		//
		// Why blend BEFORE projection (ordering matters):
		// - Projection enforces the geometric constraint "stay outside the sphere":
		//     |p - sc| = sr (or sr + thickness)
		// - The projection direction is radial: from center to the point we project.
		// - If we project using the raw integrated position, we may "accept" the full
		//   tangential drift that occurred during integration, causing more sliding.
		// - If we first blend toward prevPos, the direction used for projection is slightly
		//   biased toward last frame's contact configuration, which reduces the apparent
		//   tangential drift along the surface.
		//
		// In short:
		//   blend -> choose a slightly "more conservative" contact point
		//   project -> enforce non-penetration onto the sphere surface
		//
		// See PRACTICAL EXAMPLE and VISUALIZATION OF BLENDING EFFECT ON PROJECTION below for more insight.
        float3 bp = p * (1.0f - friction) + prevPos[pNr].xyz * friction; // lerp(pos, prevPos, friction)

		// -----------------------------------------------------------------
		// STEP 2: PROJECT OUT OF THE SPHERE (GEOMETRIC PUSH-OUT)
		//
		// r is the vector from sphere center sc to the blended point bp.
		// We rescale r so that its length becomes (sphereRadius + thickness),
		// i.e. we move the particle onto the sphere surface (expanded by thickness).
        float3 r = bp - sc;
        d = length(r);

		// Projection formula:
		//   p_new = sc + r * ((sr + thickness) / |r|)
		//
		// (This is equivalent to p_new = sc + normalize(r) * (sr + thickness), but avoids computing normalize separately)
        p = sc + r * ((sr + thickness) / d);

		// -----------------------------------------------------------------
		// PRACTICAL EXAMPLE (why blending reduces apparent tangential sliding)
		//
		// Important: the collision branch (the current if block) runs only if we are INSIDE the sphere:
		//   d = |pos - sc| < (sr + thickness)
		// So the "pos" used in this example must satisfy d < (sr + thickness) to enter the collision block.
		//
		// Consider a unit sphere centered at the origin:
		//   sc = (0, 0, 0), sr = 1
		// and ignore thickness for simplicity (thickness = 0).
		//
		// Suppose last frame the particle ended up on the sphere surface at:
		//   prevPos = (1.0, 0.0, 0.0)
		//
		// After integration (gravity/velocity), assume the particle moved to a point
		// slightly INSIDE the sphere (so we do enter the collision block), e.g.:
		//   pos = (0.9, 0.3, 0.0)
		// Check:
		//   |pos| = sqrt(0.9^2 + 0.3^2) = sqrt(0.90) ≈ 0.9499 < 1  -> penetration
		//
		// A) WITHOUT blending (friction = 0):
		//   bp = pos = (0.9, 0.3, 0.0)
		//   Project to sphere surface:
		//     newPos = normalize(bp) * sr
		//     normalize(bp) = bp / |bp| ≈ (0.9, 0.3, 0) / 0.9499 ≈ (0.947, 0.316, 0)
		//   => newPos ≈ (0.947, 0.316, 0.0)
		//
		// B) WITH blending (use friction = 0.2 to make the effect clearly visible):
		//   bp = 0.8*pos + 0.2*prevPos
		//     = 0.8*(0.9, 0.3, 0) + 0.2*(1, 0, 0)
		//     = (0.92, 0.24, 0.0)
		//   Project to sphere surface:
		//     |bp| = sqrt(0.92^2 + 0.24^2) = sqrt(0.904) ≈ 0.9508
		//     normalize(bp) ≈ (0.92, 0.24, 0) / 0.9508 ≈ (0.968, 0.252, 0)
		//   => newPos ≈ (0.968, 0.252, 0.0)
		//
		// Compare A vs B:
		// - Both results are on the sphere surface (penetration solved).
		// - With blending, the contact point advanced LESS along the surface.
		//   A quick way to see this is by comparing the polar angles (in XY plane):
		//     theta_A ≈ atan(0.316/0.947) ≈ 18.4°
		//     theta_B ≈ atan(0.252/0.968) ≈ 14.6°
		//   Blending produced a smaller angle change -> less "along-the-surface" motion.
		//
		// That reduced along-the-surface advance often looks like reduced tangential motion
		// (less sliding) during contact.
		//
		// IMPORTANT CAVEAT:
		// - This is not physically accurate Coulomb friction.
		// - It is a damping/biasing trick that *often* reduces tangential sliding visually.

		// -----------------------------------------------------------------
		// VISUALIZATION OF THE BLENDING EFFECT ON PROJECTION
		//
		// A = prevPos
		// B = pos (integrated position, inside the sphere)
		// C = bp (blended position, biased toward prevPos)
		// P_B = projection of B onto the sphere surface (without blending)
		// P_C = projection of C onto the sphere surface (with blending)

		//
		//             P_B
		//      P_C   , - ~ ~ ~ - ,
		//        , '   '           ' ,
		//  A - , - - C - B             ,
		//     ,        '  '              ,
		//    ,           ' '             ,
		//    ,             °             ,
		//    ,                           ,
		//     ,                         ,
		//      ,                       ,
		//        ,                  , '
		//          ' - , _ _ _ ,  '
		//
		// As you can see in the diagram below, P_C is closer to A than P_B is,
		// which means the blended projection results in less tangential advance along
		// the surface compared to the non-blended projection.
		// This results in a smaller reconstructed velocity
		//   (pos_final - prevPos) / dt
		// which often appears as reduced sliding (tangential damping) during contact.
    }

	// =====================================================================
	// GROUND COLLISION
	// =====================================================================
	// After sphere resolution, clamp to ground plane.
	// We enter the collision block if we are below the ground, plus a small thickness margin.
	// Use the same LERP trick (as the one used in sphere collision) for a friction-like damping on ground contact.
    if (p.y < thickness)
    {
        float3 pp = p * (1.0f - friction) + prevPos[pNr].xyz * friction;
        p = float3(pp.x, thickness, pp.z);
    }

	// Write the final integrated (and collided) position back to the pos buffer.
	// Write the raw integrated velocity to the vel buffer.
	// The velocity will be updated later (after all constraints are solved) in a separate
	// compute shader, based on the position change from prevPos to pos.
    pos[pNr] = float4(p, 0.0f);
    vel[pNr] = float4(v, 0.0f);
}

// ---------------------------------------------------------------------------
// solveConstraintsCS - Unified solver for coloring (direct write) and Jacobi
//                      ( atomic accumulation) modes.
// ---------------------------------------------------------------------------
// --------------------------------------------------------------------------
// NOTE ON "POSITION CORRECTIONS" (Δx) AND WHY JACOBI NEEDS AN ACCUMULATOR
// --------------------------------------------------------------------------
// This compute shader enforces one distance constraint per thread. For each constraint,
// it computes a POSITION CORRECTION that should be applied to the two endpoint particles
// to reduce the constraint violation (current distance - rest distance).
//
// You can think of each constraint producing two per-particle correction vectors:
//
//   Δx_id0^(c)  and  Δx_id1^(c)
//
// where:
// - "c" is the current constraint
// - id0 and id1 are the two particles connected by that constraint
// - Δx is a small displacement vector that moves a particle toward satisfying the constraint
//
// In the code below these are:
//   dP = n * (l - l0) / (w0 + w1)
//   Δx_id0^(c) = c0 = + w0 * dP
//   Δx_id1^(c) = c1 = - w1 * dP
//
// w0 and w1 are inverse masses, so heavier particles move less, since the inverse mass is smaller.
// n is the normalized direction of the edge, from particle id0 to id1, so corrections are along the line connecting the particles.
// (l - l0) is the constraint violation: how much the current distance (l) exceeds the rest length (l0).
// The division by (w0 + w1) ensures that the total correction is distributed according to the particles' masses.
// The sign convention is such that if the current distance l is greater than the rest length l0, the correction
// will pull the particles together. If l is less than l0, the correction will push them apart.
// dP is the base correction vector that represents how much we need to adjust the positions along the edge direction to satisfy the constraint.
// The individual particle corrections are then scaled by their inverse masses to determine how much each particle should move.
//
// ---------------------------------------------------------------------
// Gauss-Seidel-like mode (solveType==0): update positions immediately (in-place)
// ---------------------------------------------------------------------
// If solveType==0 in C++, the jacobi scale is zero in the first 4 passes, which means
// we apply each constraint's correction directly to pos[]:
//
//   pos[id0] += Δx_id0^(c)
//   pos[id1] += Δx_id1^(c)
//
// This is similar in spirit to Gauss-Seidel iteration: as soon as a correction is computed,
// it is written back, so subsequent constraints may see updated positions.
//
// This is safe only when the constraints in the current pass are independent (do not share
// particles), otherwise many constraints would try to update the same pos[] entry in parallel.
//
// Subsequent passes can then build on the updated positions from previous passes,
// allowing for a more "immediate" propagation of corrections.
// In this case, the convergence to the final solution is often faster (fewer iterations needed)
// because each constraint can see the latest positions.
//
// ---------------------------------------------------------------------
// Jacobi mode (solveType==1): accumulate per-particle corrections into corrections[]
// ---------------------------------------------------------------------
// If solveType==1 in C++, the jacobi scale is NOT zero, which means we CANNOT directly modify pos[] in this compute shader.
// Instead, we accumulate all per-constraint contributions in a separate array corrections[]:
//
//   corcorrections[id0] += Δx_id0^(c)
//   cocorrectionsr[id1] += Δx_id1^(c)
//
// Then a second compute shader (addCorrections) applies the sum once per particle:
//
//   pos[p] += corrections[p] * jacobiScale
//
// This accumulator is necessary because a single particle usually participates in MANY
// constraints at once (especially shear/bending in this sample).
// In Jacobi we want the update for particle p to be the SUM of ALL contributions from all constraints that involve p,
// computed from the SAME original/old positions.
//
//   pos_new[p] = pos_old[p] + Σ_c Δx_p^(c)
//
// In contrast to Gauss-Seidel, Jacobi's corrections are based on the same original positions,
// so they can be computed in parallel without read/write conflicts.
// However, Jacobi typically requires more iterations to converge compared to Gauss-Seidel because it
// does not immediately see the effects of other constraints' corrections.
// The choice between Gauss-Seidel and Jacobi often depends on the structure of the constraints and
// the desired convergence properties.
//
// ----------------------------------------------------------------------
// HYBRID APPROACH IN THIS SIMULATOR
// ----------------------------------------------------------------------
// - if solveType==0 in C++, the jacobi scale is zero for passes 0-3, which are independent, so
//   we use the Gauss-Seidel-like approach to solve each of those passes with direct updates to pos[],
//   in parallel across constraints in the same pass, but without conflicts since they are independent.
//   For pass 4, which includes shear and bending constraints that create dependencies, the jacobi scale
//   is non-zero, which means we switch to the Jacobi approach and accumulate corrections into corrections[]
//   instead of writing to pos[] to avoid conflicts. Then we apply those corrections in a separate step.
//   However, note that constraints in pass 4 are NOT independent, so we CANNOT accumulate corrections in-place
//   without conflicts: we need atomic adds to safely sum contributions from multiple constraints that share particles.
//   This hybrid approach allows us to maximize performance for the independent passes while
//   ensuring correctness and stability for the dependent pass.
//
// - if solveType==1
//   we use Jacobi for all passes, which is simpler but may require more iterations overall.
//   This is useful if the constraints are not easily partitioned into independent sets.
//
// ---------------------------------------------------------------------
// Concrete 1D example (intuitive sign/meaning of Δx terms)
// ---------------------------------------------------------------------
// Consider three particles on a line: A -- B -- C (1D positions).
// Particle B participates in two constraints: (A,B) and (B,C).
//
// Suppose during one iteration (computed from the SAME old positions) we get:
//   ΔB^(AB) = -1.0   (constraint AB wants to move B left by 1 unit)
//   ΔB^(BC) = +3.0   (constraint BC wants to move B right by 3 units)
//
// In Jacobi, both contributions must be combined:
//   corrections[B] = ΔB^(AB) + ΔB^(BC) = (-1) + (+3) = +2
//   B_new = B_old + corrections[B]
//
// So B moves right by 2 units overall, representing the net effect of both constraints.
//
// Note: In Jacobi, constraints are NOT immediately satisfied each iteration.
// Instead, corrections are accumulated and applied together at the end of the iteration.
// The next iteration will recompute corrections from the new positions.
// All corrections for a given particle are computed from the SAME original positions.
// This is why we need atomic operations: many constraints may update corrections[B]
// in parallel, and we must safely accumulate all contributions before the apply step.
// In THIS code, those Δ contributions correspond to the terms written via atomic adds:
// corrections[B] gets +wB*dP from one constraint and -wB*dP from the other, and we sum them.
//
// ---------------------------------------------------------------------
[numthreads(CLOTH_THREAD_GROUP_SIZE, 1, 1)]
void cloth_solveConstraintsCS(uint3 DTid : SV_DispatchThreadID)
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
		// Many constraints may touch the same particle -> atomic accumulation is required.
		// See cloth_common.hlsli for the ATOMIC_ADD_FLOAT macro, which implements a CAS-loop to
		// atomically add a float value to a uint buffer storing float bits.
		// Note that we store 3 floats per particle in the corrections buffer for x, y, z components,
		// so the index is 3*id + component.
        ATOMIC_ADD_FLOAT(corrections, 3u * id0 + 0u, c0.x);
        ATOMIC_ADD_FLOAT(corrections, 3u * id0 + 1u, c0.y);
        ATOMIC_ADD_FLOAT(corrections, 3u * id0 + 2u, c0.z);
        ATOMIC_ADD_FLOAT(corrections, 3u * id1 + 0u, c1.x);
        ATOMIC_ADD_FLOAT(corrections, 3u * id1 + 1u, c1.y);
        ATOMIC_ADD_FLOAT(corrections, 3u * id1 + 2u, c1.z);
    }
    else
    {
		// Constraints in this pass are independent, so we can directly update positions without conflicts.
        pos[id0] = float4(p0 + c0, 0.0f);
        pos[id1] = float4(p1 + c1, 0.0f);
    }
}

// ---------------------------------------------------------------------------
// addCorrectionsCS - Apply accumulated Jacobi corrections then zero the buffer.
// ---------------------------------------------------------------------------
//
// JACOBI APPLICATION STEP: pos += corr * jacobiScale
//
// The parameter "scale" (jacobiScale in the demo) is an UNDER-RELAXATION factor:
// - scale = 1.0 would apply the full Jacobi correction each iteration.
// - scale < 1.0 applies only a fraction of the correction.
//
// Why use scale < 1?
// - Jacobi updates are computed from the same "old" positions and then applied together.
//   In strongly coupled constraint networks (cloth), the summed correction can be large,
//   and applying it fully may overshoot, causing oscillations/jitter.
// - Under-relaxation reduces overshoot and improves stability, at the cost of needing
//   more iterations/substeps to converge.
//
// Practical effect when tuning:
// - Higher scale (closer to 1): faster convergence but more risk of instability/jitter.
// - Lower scale (e.g. 0.1..0.3): more stable but "softer" unless you increase iterations.
// ---------------------------------------------------------------------
[numthreads(CLOTH_THREAD_GROUP_SIZE, 1, 1)]
void cloth_addCorrectionsCS(uint3 DTid : SV_DispatchThreadID)
{
	// Each thread processes one particle, applying the sum of all corrections accumulated for that particle.
    uint pNr = DTid.x;
    if (pNr >= cb.numParticles)
        return;

	// Read the accumulated correction for this particle.
	// Since we used atomic adds to accumulate into a uint buffer, we need to reinterpret the bits back to float.
	// Also, we stored 3 floats per particle (x, y, z) in the corrections buffer, so we read 3 consecutive uints and reinterpret them as floats.
    uint base = 3u * pNr;
    float cx = asfloat(corrections[base + 0u]);
    float cy = asfloat(corrections[base + 1u]);
    float cz = asfloat(corrections[base + 2u]);

	// Zero the correction after reading, so it's ready for the next iteration.
    corrections[base + 0u] = 0;
    corrections[base + 1u] = 0;
    corrections[base + 2u] = 0;

	// Multiply the accumulated correction by the jacobi scale factor before applying.
    float scale = cb.jacobiScale;
    float3 correction = float3(cx, cy, cz) * scale;

	// Optional: clamp the correction to prevent extreme jumps that could cause instability.
    // float corrLen = length(correction);
    // const float maxCorr = 0.1f;
    // if (corrLen > maxCorr)
    //     correction = correction * (maxCorr / corrLen);

	// Apply the scaled correction to the position.
	// We don't need to worry about conflicts here because each particle is processed by only one thread,
	// and we've already safely accumulated all contributions into the corrections buffer.
    float3 p = pos[pNr].xyz;
    pos[pNr] = float4(p + correction, 0.0f);
}

// ---------------------------------------------------------------------------
// updateVelCS - Update velocity from position change: vel = (pos - prevPos) / dt
// ---------------------------------------------------------------------------
[numthreads(CLOTH_THREAD_GROUP_SIZE, 1, 1)]
void cloth_updateVelCS(uint3 DTid : SV_DispatchThreadID)
{
	// Each thread processes one particle, updating its velocity based on the change in position from prevPos to pos.
    uint pNr = DTid.x;
    if (pNr >= cb.numParticles)
        return;

    float invDt = 1.0f / cb.dt;

	// If this particle is being dragged by the mouse,
	// we set its velocity to zero to avoid a sudden jump when released.
    if ((int)pNr == cb.dragParticleNr)
    {
        vel[pNr] = float4(0.0f, 0.0f, 0.0f, 0.0f);
        return;
    }

	// Compute velocity as the change in position divided by the time step.
	// This is a simple finite difference approximation of velocity based on the new and previous positions.
	// Note that this overwrites the velocity computed during integration, which is intentional because after solving constraints,
	// the actual velocity should reflect the final position change, not just the integrated velocity before constraints.
	// We don't need to worry about conflicts here because each particle is processed by only one thread.
    vel[pNr] = float4((pos[pNr].xyz - prevPos[pNr].xyz) * invDt, 0.0f);
}

// ---------------------------------------------------------------------------
// clearNormalsCS - Zero the normAccum buffer before normal accumulation.
// ---------------------------------------------------------------------------
[numthreads(CLOTH_THREAD_GROUP_SIZE, 1, 1)]
void cloth_clearNormalsCS(uint3 DTid : SV_DispatchThreadID)
{
    uint pNr = DTid.x;
    if (pNr >= cb.numParticles)
        return;

    uint base = 3u * pNr;
    normAccum[base + 0u] = 0;
    normAccum[base + 1u] = 0;
    normAccum[base + 2u] = 0;
}

// ---------------------------------------------------------------------------
// addNormalsCS - Accumulate face normals per vertex via CAS-loop float atomics.
//                One thread per triangle.
// ---------------------------------------------------------------------------
[numthreads(CLOTH_THREAD_GROUP_SIZE, 1, 1)]
void cloth_addNormalsCS(uint3 DTid : SV_DispatchThreadID)
{
	// Each thread processes one triangle, identified by the global thread index DTid.x.
	// Here numConstraintsInPass is repurposed to store the number of triangles for normal accumulation,
	// which is set in C++ before dispatch.
    uint triNr = DTid.x;
    if (triNr >= cb.numConstraintsInPass)
        return;

	// Read the vertex indices of the triangle from the triIds buffer.
    uint id0 = triIds[3u * triNr];
    uint id1 = triIds[3u * triNr + 1u];
    uint id2 = triIds[3u * triNr + 2u];

	// Compute the face normal using the cross product of two edges of the triangle.
	// We are in a left-handed coordinate system, where the resulting vector from the cross product "see"
	// the first vector as rotating toward the second vector in a clockwise manner.
	// However, the winding order of the triangle vertices (id0, id1, id2) is counter-clockwise,
	// so the cross product becomes n = cross(e2, e1) instead of n = cross(e1, e2) to get the correct normal direction.
    float3 e1 = pos[id1].xyz - pos[id0].xyz;
    float3 e2 = pos[id2].xyz - pos[id0].xyz;
    float3 n = cross(e2, e1);

	// Accumulate the face normal n into the normAccum buffer for each of the three vertices.
	// Since multiple triangles can share the same vertex, we must use atomic operations to
	// safely accumulate the normals without race conditions.
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

// ---------------------------------------------------------------------------
// normalizeNormalsCS - Read accumulated normals, normalize, write to normals buffer.
// ---------------------------------------------------------------------------
[numthreads(CLOTH_THREAD_GROUP_SIZE, 1, 1)]
void cloth_normalizeNormalsCS(uint3 DTid : SV_DispatchThreadID)
{
	// Each thread processes one particle.
    uint pNr = DTid.x;
    if (pNr >= cb.numParticles)
        return;

	// Read the accumulated normal for this particle from the normAccum buffer.
	// Since we used atomic adds to accumulate into a uint buffer, we need to reinterpret the bits back to float.
    uint base = 3u * pNr;
    float nx = asfloat(normAccum[base + 0u]);
    float ny = asfloat(normAccum[base + 1u]);
    float nz = asfloat(normAccum[base + 2u]);

	// Normalize the accumulated normal vector to get the final vertex normal, storing it in the normals buffer.
	// We don't need to worry about conflicts here because each particle is processed by only one thread.
    normals[pNr] = float4(normalize(float3(nx, ny, nz)), 0.0f);
}

// ---------------------------------------------------------------------------
// raycastTriangleCS - GPU raycast: Moller-Trumbore per-triangle intersection.
// ---------------------------------------------------------------------------
//
// GENERIC DESCRIPTION:
// This compute shader performs ray-triangle intersection testing on the GPU in parallel.
// For each triangle in the cloth mesh, it determines if the ray (from the mouse click)
// intersects that triangle, and if so, computes the distance from the ray origin to
// the intersection point. The result allows us to find which triangle was clicked.
//
// ALGORITHM: Moller-Trumbore Ray-Triangle Intersection
// The steps are as follows:
// 1. Check if ray is parallel to triangle plane (det == 0) -> NO HIT
// 2. Compute barycentric u: is intersection "within edge1 range"? (0 ≤ u ≤ 1)
// 3. Compute barycentric v: is intersection "within edge2 range"? (0 ≤ v ≤ 1)
// 4. Check if u + v ≤ 1: is intersection within triangle boundaries?
// 5. If all checks pass, compute t = distance along ray to intersection point
// 6. Return t (distance) to use in ClothMesh::StartGrabGPU for selecting the closest triangle
//
// WHY THIS ALGORITHM?
// - Very efficient: only one division per test (vs. multiple divisions in other methods)
// - Naturally gives us barycentric coordinates (u, v) to check if point is inside triangle
// - Perfect for GPU parallelization: each thread tests one triangle independently
//
// ---------------------------------------------------------------------
// THE GEOMETRY PROBLEM (from first principles):
//
// Ray equation:
//   R(t) = O + t * D
// where:
//   O = ray origin (orig)
//   D = ray direction (dir), typically normalized
//   t = scalar distance along the ray
//
// Triangle vertices:
//   V0, V1, V2
//
// Define triangle edges (a local 2D basis on the triangle plane):
//   E1 = V1 - V0
//   E2 = V2 - V0
//
// The parametric equation of a plane is defined by a point and two spanning vectors.
// So, we can use V0 as the point and E1, E2 as the spanning vectors.
// Any point P on the triangle plane can be written as:
//   P(u, v) = V0 + u * E1 + v * E2
//
// For P to be INSIDE the triangle (barycentric constraints):
//   u >= 0
//   v >= 0
//   u + v <= 1
//
// Ray hits the triangle if ray equation and triangle plane equation have a common solution.
// That is, there exist scalars (t, u, v) such that:
//   O + t * D = V0 + u * E1 + v * E2
//
// Rearranging:
//   O - V0 = u * E1 + v * E2 - t * D
//
// Let:
//   T = O - V0
//
// Then:
//   T = u * E1 + v * E2 - t * D
//
// This is a 3D vector equation (3 components) with 3 unknown scalars (u, v, t).
// It is equivalent to solving a 3x3 linear system
//
//   | E1  E2  -D |  *  | u  v  t |^T = | T |
//
// or in expanded form:
//
//   | E1.x  E2.x  -D.x |   | u | = | T.x |
//   | E1.y  E2.y  -D.y | * | v | = | T.y |
//   | E1.z  E2.z  -D.z |   | t | = | T.z |
//
// The determinant of the coefficient matrix tells us if the system has a
// unique solution (det != 0) or not (det == 0).
// Möller–Trumbore algirithm avoids explicitly solving this system using matrices and
// instead directly computes u, v, t using dot and cross products,
// which is more efficient and numerically stable.
//
// ---------------------------------------------------------------------
// WHY DOT AND CROSS PRODUCTS APPEAR (key idea):
//
// The scalar triple product:
//   a . (b x c)
//
// equals the oriented volume of the parallelepiped spanned by (a, b, c).
// If it is 0, the vectors are coplanar -> the corresponding 3x3 system is singular
// (no unique solution). In ray/triangle terms, this corresponds to the ray being
// parallel (or nearly parallel) to the triangle plane.
//
// In this algorithm we compute (derived and explained below):
//   det = E1 . (D x E2)
//
// This "det" plays the role of the determinant of the implicit 3x3 system.
// - det == 0  -> ray parallel to triangle plane (or degenerate triangle) -> no hit
// - |det| big -> stable intersection computation
//
// We then compute inv_det = 1/det once, because u, v, t are all
//
// (some_value) * inv_det
//
// (explained in the next section).
//
// ---------------------------------------------------------------------
// HOW u, v, t FORMULAS ARE OBTAINED:
//
// Starting from:
//   T = uE1 + vE2 - tD
//
// If we dot both sides with (D x E2), two terms vanish automatically:
// - E2 . (D x E2) = 0 because (D x E2) is perpendicular to E2
// - D  . (D x E2) = 0 because (D x E2) is perpendicular to D
//
// So we isolate u:
//   T . (D x E2) = u * E1 . (D x E2)
//   u = [T . (D x E2)] / [E1 . (D x E2)]
//
// That is exactly (used in the code):
//   pv = D x E2
//   det  = E1 . pvec
//   u    = (T . pvec) / det
//
// Similarly, we compute:
//   qv = T x E1
// and obtain:
//   v = [D . (T x E1)] / det  = (D . qvec) / det
//   t = [E2 . (T x E1)] / det = (E2 . qvec) / det
//
// ---------------------------------------------------------------------
// INSIDE-TRIANGLE TEST (why these checks):
//
// The barycentric constraints:
//   u >= 0, v >= 0, u + v <= 1
// are exactly the condition that the intersection point lies inside the triangle.
//
// The code checks:
//   if u < 0 or u > 1 -> outside
//   if v < 0 or u + v > 1 -> outside
//
// ---------------------------------------------------------------------
// NUMERICAL NOTE
//
// In the current explanation we check "det == 0.0" to determine if the ray is parallel to the triangle plane.
// However, in practice, due to floating point precision issues, we should check if "abs(det) < eps" for some
// small epsilon value (e.g., 1e-9f) instead of checking for exact equality to zero because floating point
// rarely hits exactly 0.0 and near-parallel cases can be unstable.
//
// ---------------------------------------------------------------------
// CONCRETE NUMERICAL EXAMPLE (sanity check):
//
// Triangle on z=0 plane:
//   V0=(0,0,0), V1=(1,0,0), V2=(0,1,0)
//   E1=(1,0,0), E2=(0,1,0)
//
// Ray:
//   O=(0.25,0.25,1), D=(0,0,-1)
//
// pvec = DxE2 = (0,0,-1)x(0,1,0) = (1,0,0)
// det  = E1.pvec = (1,0,0).(1,0,0) = 1
// T    = O-V0 = (0.25,0.25,1)
// u    = T.pvec / det = 0.25
// qvec = TxE1 = (0.25,0.25,1)x(1,0,0) = (0,1,-0.25)
// v    = D.qvec / det = (0,0,-1).(0,1,-0.25) = 0.25
// u+v  = 0.5 <= 1 -> inside
// t    = E2.qvec / det = (0,1,0).(0,1,-0.25) = 1
// Intersection point:
//   P = O + tD = (0.25,0.25,1) + 1*(0,0,-1) = (0.25,0.25,0)
[numthreads(CLOTH_THREAD_GROUP_SIZE, 1, 1)]
void cloth_raycastTriangleCS(uint3 DTid : SV_DispatchThreadID)
{
	// Each thread processes one triangle, identified by the global thread index DTid.x.
    uint triNr = DTid.x;
    uint numTris;
    uint stride;

	// Get the number of triangles from the triDist buffer,
	// which is set in C++ code before dispatch (see ClothMesh::CreateGPUBuffers).
	// If triNr exceeds the number of triangles, we return early.
    triDist.GetDimensions(numTris, stride);
    if (triNr >= numTris)
        return;

	// Set a default "no hit" distance value, which is a large number that indicates the ray did not intersect the triangle.
    const float noHit = 1e9f;
    float3 orig = float3(rayCB.origX, rayCB.origY, rayCB.origZ);
    float3 dir  = float3(rayCB.dirX,  rayCB.dirY,  rayCB.dirZ);

	// Read the vertex indices of the triangle from the triIds buffer.
    uint id0 = triIds[3u * triNr];
    uint id1 = triIds[3u * triNr + 1u];
    uint id2 = triIds[3u * triNr + 2u];

	// Read the vertex positions of the triangle from the pos buffer.
    float3 p0 = pos[id0].xyz;
    float3 p1 = pos[id1].xyz;
    float3 p2 = pos[id2].xyz;

	// Compute the edges of the triangle.
    float3 e1 = p1 - p0;
    float3 e2 = p2 - p0;

	// Compute pv and det (as explained in the algorithm description).
    float3 pv = cross(dir, e2);
    float det = dot(e1, pv);

	// Check if the ray is parallel to the triangle plane (det == 0) or nearly parallel (abs(det) < epsilon).
	// If so, we consider it a "no hit" case and write the noHit value to the triDist buffer for this triangle.
    if (abs(det) < 1e-9f)
    {
        triDist[triNr] = noHit;
        return;
    }

	// Compute the inverse of det once, since we will use it to compute u, v, and t.
    float inv = 1.0f / det;

	// Compute tv (T in the algorithm description) and then compute u.
    float3 tv = orig - p0;
    float u = dot(tv, pv) * inv;

	// Check if u is outside the range [0, 1].
	// If so, the intersection point is outside the triangle, and we write noHit.
    if (u < 0.0f || u > 1.0f)
    {
        triDist[triNr] = noHit;
        return;
    }

	// Compute qv and then compute v.
    float3 qv = cross(tv, e1);
    float v = dot(dir, qv) * inv;

	// Check if v is outside the range [0, 1] or if u + v > 1.
	// If so, the intersection point is outside the triangle, and we write noHit.
    if (v < 0.0f || u + v > 1.0f)
    {
        triDist[triNr] = noHit;
        return;
    }

	// Compute d (the distance along the ray to the intersection poin; t in the algorithm description).
    float d = dot(e2, qv) * inv;

	// If d is negative or zero, the intersection point is behind the ray origin, so we consider it a "no hit".
    if (d <= 0.0f)
    {
        triDist[triNr] = noHit;
        return;
    }

	// If we reach this point, we have a valid intersection.
	// We write the distance d to the triDist buffer for this triangle.
    triDist[triNr] = d;
}
