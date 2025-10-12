#include <algorithm>
#include <cmath>
#include <vector>
#include <wiScene.h>
#include "softBody.h"
#include "bunnyMesh.h"
#include "simulation_utils.h"

SoftBody::SoftBody(
    const TetMesh& tetMesh,
    float edgeCompliance_,
    float volCompliance_
)
    : edgeCompliance(edgeCompliance_)
    , volCompliance(volCompliance_)
    , numVerts(static_cast<int>(tetMesh.verts.size() / 3))
    , numTets(static_cast<int>(tetMesh.tetIds.size() / 4))
    , grads(4 * 3, 0.0f)
	, invMass(numVerts, 0.0f)
    , restVol(numTets, 0.0f)
    , edgeLengths(tetMesh.tetEdgeIds.size() / 2, 0.0f)
{
    pos.reserve(numVerts);
    prevPos.reserve(numVerts);
    for (int i = 0; i < numVerts; ++i) {
        XMFLOAT3 v(
            tetMesh.verts[3 * i + 0],
            tetMesh.verts[3 * i + 1],
            tetMesh.verts[3 * i + 2]
        );
        pos.push_back(v);
        prevPos.push_back(v);
    }

	tetIds.reserve(numTets);
	for (int i = 0; i < numTets; ++i) {
		XMINT4 t(
			tetMesh.tetIds[4 * i + 0],
			tetMesh.tetIds[4 * i + 1],
			tetMesh.tetIds[4 * i + 2],
			tetMesh.tetIds[4 * i + 3]
		);
		tetIds.push_back(t);
	}

	edgeIds.reserve(tetMesh.tetEdgeIds.size() / 2);
	for (size_t i = 0; i < tetMesh.tetEdgeIds.size() / 2; ++i) {
		XMINT2 e(
			tetMesh.tetEdgeIds[2 * i + 0],
			tetMesh.tetEdgeIds[2 * i + 1]
		);
		edgeIds.push_back(e);
	}

	surfaceTriIds.reserve(tetMesh.tetSurfaceTriIds.size() / 3);
	for (size_t i = 0; i < tetMesh.tetSurfaceTriIds.size() / 3; ++i) {
		XMINT3 tri(
			tetMesh.tetSurfaceTriIds[3 * i + 0],
			tetMesh.tetSurfaceTriIds[3 * i + 1],
			tetMesh.tetSurfaceTriIds[3 * i + 2]
		);
		surfaceTriIds.push_back(tri);
	}

	vel.resize(3 * numVerts, XMFLOAT3(0,0,0));
}

void SoftBody::InitPhysics()
{
	// Initialize inverse mass of vertices and the rest volume of tetrahedra
    std::fill(invMass.begin(), invMass.end(), 0.0f);
    std::fill(restVol.begin(), restVol.end(), 0.0f);

	// For each tetrahedron:
	// compute its rest volume and the inverse mass contribution to its four vertices.
	// Assume uniform density = 1.0 so that mass = volume (from mass = density * volume).
	// Each vertex gets 1/4 of the tetrahedron's mass so that mass is evenly distributed.
	// A vertex can belong to multiple tetrahedra, so we accumulate the inverse mass contributions.
	// Note: if a solid is composed of two tetrahedra, its volume and mass will be bigger than that
	// of a single tetrahedron. Consequently, the mass of a vertex shared by both tetrahedra
	// will be the sum of the mass contributions from each tetrahedron.
    for (int i = 0; i < numTets; ++i)
    {
        float vol = GetTetVolume(i);
        restVol[i] = vol;
        float pInvMass = (vol > 0.0f) ? (1.0f / (vol / 4.0f)) : 0.0f;
        const XMINT4& tet = tetIds[i];
        invMass[tet.x] += pInvMass;
        invMass[tet.y] += pInvMass;
        invMass[tet.z] += pInvMass;
        invMass[tet.w] += pInvMass;
    }

	// Compute the rest length of each edge in the mesh
    for (size_t i = 0; i < edgeLengths.size(); ++i)
    {
        const XMINT2& id = edgeIds[i];
        XMVECTOR va = XMLoadFloat3(&pos[id.x]);
        XMVECTOR vb = XMLoadFloat3(&pos[id.y]);
        edgeLengths[i] = XMVectorGetX(XMVector3Length(XMVectorSubtract(va, vb)));
    }
}

void SoftBody::Translate(const XMFLOAT3& delta)
{
    XMVECTOR vdelta = XMLoadFloat3(&delta);
    for (int i = 0; i < numVerts; ++i)
    {
        XMVECTOR vpos = XMLoadFloat3(&pos[i]);
        vpos = XMVectorAdd(vpos, vdelta);
        XMStoreFloat3(&pos[i], vpos);

        XMVECTOR vprev = XMLoadFloat3(&prevPos[i]);
        vprev = XMVectorAdd(vprev, vdelta);
        XMStoreFloat3(&prevPos[i], vprev);
    }
}

void SoftBody::Scale(const XMFLOAT3& scale)
{
	XMMATRIX scaleMat = XMMatrixScaling(scale.x, scale.y, scale.z);

	for (int i = 0; i < numVerts; ++i)
	{
		// Scale pos
		XMVECTOR p = XMVectorSet(pos[i].x, pos[i].y, pos[i].z, 0.0f);
		p = XMVector3Transform(p, scaleMat);
		XMStoreFloat3(&pos[i], p);

		// Scale prevPos
		XMVECTOR pp = XMVectorSet(prevPos[i].x, prevPos[i].y, prevPos[i].z, 0.0f);
		pp = XMVector3Transform(pp, scaleMat);
		XMStoreFloat3(&prevPos[i], pp);
	}
}

void SoftBody::RotateRollPitchYaw(const XMFLOAT3& angle)
{
	XMMATRIX rot = XMMatrixRotationRollPitchYaw(angle.y, angle.z, angle.x);

	for (int i = 0; i < numVerts; ++i)
	{
		// Rotate pos
		XMVECTOR p = XMVectorSet(pos[i].x, pos[i].y, pos[i].z, 0.0f);
		p = XMVector3Transform(p, rot);
		XMStoreFloat3(&pos[i], p);

		// Rotate prevPos
		XMVECTOR pp = XMVectorSet(prevPos[i].x, prevPos[i].y, prevPos[i].z, 0.0f);
		pp = XMVector3Transform(pp, rot);
		XMStoreFloat3(&prevPos[i], pp);
	}
}

void SoftBody::Squash()
{
	for (int i = 0; i < numVerts; i++)
	{
		pos[i].y = 0.5f;
		prevPos[i].y = 0.5f;
	}
	if (!gPhysicsScene.paused)
		gPhysicsScene.Run();
}

float SoftBody::GetTetVolume(int nr)
{
	const XMINT4& t = tetIds[nr];
    XMVECTOR v0 = XMLoadFloat3(&pos[t.x]);
    XMVECTOR v1 = XMLoadFloat3(&pos[t.y]);
    XMVECTOR v2 = XMLoadFloat3(&pos[t.z]);
    XMVECTOR v3 = XMLoadFloat3(&pos[t.w]);

    XMVECTOR d1 = XMVectorSubtract(v1, v0);
    XMVECTOR d2 = XMVectorSubtract(v2, v0);
    XMVECTOR d3 = XMVectorSubtract(v3, v0);

	// From linear algebra:
	// The volume of a parallelepiped specified by three vectors is given by the triple scalar product of the vectors.
	// The scalar triple product is unchanged under a circular shift of its three operands:
	// A * (B x C) = B * (C x A) = C * (A x B).
	// In this case, A = (v1 - v0), B = (v2 - v0) and C = (v3 - v0).
	// Tetrahedron volume is 1/6 of the parallelepiped volume:
	// V = 1/6 * (C * (A x B))
	// This is from the fact a tetrahedron is a pyramid with a triangular base, and the volume of a pyramid is
	// 1/3 * base_area * height.
	// In this case, the base area is 1/2 * |A x B| and the height is the projection of C onto the normal of the base, that is
	// base_area = 1/2 * |A x B|
	// height = C * ((A x B) / |A x B|) = (C * (A x B)) / |A x B|
	// So, putting all together:
	// V = 1/3 * base_area * height = 1/3 * (1/2 * |A x B|) * ((C * (A x B)) / |A x B|)
	//   = 1/6 * |A x B| * (C * (A x B)) / |A x B|
	//   = 1/6 * (C * (A x B))
    XMVECTOR cross = XMVector3Cross(d1, d2);
    float dot = XMVectorGetX(XMVector3Dot(cross, d3));
    return dot / 6.0f;
}

void SoftBody::PreSolve(float dt, const XMFLOAT3& gravity)
{
    XMVECTOR grav = XMLoadFloat3(&gravity);

	// For each vertex:
    for (int i = 0; i < numVerts; i++)
    {
		// Skip fixed vertices
        if (invMass[i] == 0.0f)
            continue;

        // Load vectors
        XMVECTOR vel_i = XMLoadFloat3(&vel[i]);
        XMVECTOR pos_i = XMLoadFloat3(&pos[i]);
        XMVECTOR prev_i = XMLoadFloat3(&prevPos[i]);

        // Explicit integration:
        vel_i = XMVectorAdd(vel_i, XMVectorScale(grav, dt));  // v += g * dt (new velocity)
        prev_i = pos_i;                                       // prev = pos (current position stored as previous)
        pos_i = XMVectorAdd(pos_i, XMVectorScale(vel_i, dt)); // pos += v * dt (new predicted position)

        // Ground collision:
		// If the new position is below y = 0, reset it to the previous position and clamp y = 0.
        if (XMVectorGetY(pos_i) < 0.0f)
		{
            pos_i = prev_i;
            pos_i = XMVectorSetY(pos_i, 0.0f);
        }

        // Store results
        XMStoreFloat3(&vel[i], vel_i);
        XMStoreFloat3(&pos[i], pos_i);
        XMStoreFloat3(&prevPos[i], prev_i);
    }
}

void SoftBody::SolveEdges(float compliance, float dt)
{
    float alpha = compliance / (dt * dt);

	// For each edge:
    for (size_t i = 0; i < edgeLengths.size(); ++i)
    {
		// Compute the inverse mass sum of the two vertices of the edge
        const XMINT2& id = edgeIds[i];
        int id0 = id.x;
        int id1 = id.y;
        float w0 = invMass[id0];
        float w1 = invMass[id1];
        float w = w0 + w1;
        if (w == 0.0f) continue;

        // Edge predicted positions
        XMVECTOR p0 = XMLoadFloat3(&pos[id0]);
        XMVECTOR p1 = XMLoadFloat3(&pos[id1]);

        // Current length
        XMVECTOR diff = XMVectorSubtract(p0, p1);
        float len = XMVectorGetX(XMVector3Length(diff));
        if (len == 0.0f) continue;

        // Normalized direction vector of the edge (from p1 towards p0)
        XMVECTOR dir = XMVectorScale(diff, 1.0f / len);

        // Retrieve the current value of the edge constraint and compute the scaling factor to apply to the positions.
		// For an edge (distance) constraint C, the gradient of C is a unit vector along the edge,
		// so its squared modulo is 1 and the denominator simplifies to w0 + w1 (the sum of inverse masses).
		// That’s why here the gradient does not appear explicitly.
		//
		// Hybrid PBD/XPBD approach:
		// compute the scaling factor s as in PBD, but add the compliance term in the denominator as in XPBD.
		// This incorporates a physically meaningful stiffness into the constraint projection, instead of
		// tuning stiffness empirically via iterations as in plain PBD.
		//
		// If C(p0, p1) = length(p0 - p1) - restLen = |p0 - p1| - restLen,
		// with p0 = (x0, y0, z0) and p1 = (x1, y1, z1).
		// We can define d = (dx, dy, dz) = p0 - p1 = (x0 - x1, y0 - y1, z0 - z1) and |d| = r = sqrt((dx)^2 + (dy)^2 + (dz)^2)
		// so that
		// C(p0, p1) = r - restLen
		//
		// The gradient of C w.r.t. p0 and p1 is given by the partial derivatives:
		// dC/dp0 = (dC/dx0, dC/dy0, dC/dz0)
		// dC/dp1 = (dC/dx1, dC/dy1, dC/dz1)
		//
		// To compute the partial derivatives, we use the chain rule. For example, for dC/dx0:
		// let's start from the following definitions:
		// dx = x0 - x1
		// f = (dx^2 + dy^2 + dz^2)
		// r = sqrt(f)
		// and remember that the derivative of sqrt(x) = x^(1/2) is d(sqrt)/dx = (1/2) * x^(-1/2) = 1/(2 * sqrt(x))
		// and also that the derivative of a constant (like restLen) is zero.
		// So we have:
		// dC/dx0 = dr/dx0 - d(restLen)/dx0 = dr/dx0 - 0
		// To compute dr/dx0, we use the chain rule passing through f:
		// dr/dx0 = dr/df * df/dx0 = (1/(2 * sqrt(f))) * df/dx0 = (1/(2 * r)) * df/dx0
		//        = (1/(2 * r)) * d(dx^2 + dy^2 + dz^2)/dx0 = (1/(2 * r)) * (d(dx^2)/dx0 + d(dy^2)/dx0 + d(dz^2)/dx0)
		//        = (1/(2 * r)) * (2 * dx * d(dx)/dx0 + 0 + 0) = (1/(2 * r)) * (2 * dx * d(x0 - x1)/dx0)
		//        = (1/(2 * r)) * (2 *  dx * dx0/dx0 - 2 * dx * dx1/dx0) = (1/(2 * r)) * (2 * dx * 1 - 2 * dx * 0)
		//        = (1/(2 * r)) * (2 * dx) = dx / r = (x0 - x1) / |p0 - p1|
		//
		// The same procedure can be applied to compute the other partial derivatives. So, we have:
		// dC/dx0 = (x0 - x1) / |p0 - p1|
		// dC/dy0 = (y0 - y1) / |p0 - p1|
		// dC/dz0 = (z0 - z1) / |p0 - p1|
		// dC/dx1 = -(x0 - x1) / |p0 - p1|
		// dC/dy1 = -(y0 - y1) / |p0 - p1|
		// dC/dz1 = -(z0 - z1) / |p0 - p1|
		//
		// So, in vector form, we can write:
		// dC/dp0 = (p0 - p1) / |p0 - p1| = dir,
		// dC/dp1 = -(p0 - p1) / |p0 - p1| = -dir.
		// with dir = p0 - p1 / |p0 - p1| = (x0 - x1 / |p0 - p1|, (y0 - y1) / |p0 - p1|, (z0 - z1) / |p0 - p1|)
		//
		// The formula for the scaling factor s is:
		// s = -C(p0, p1) / (w0 * |dC/dp0|^2 + w1 * |dC/dp1|^2 + alpha)
		// but since |dC/dp0| = |dC/dp1| = |dir| = 1 (unit vectors), the formula simplifies to:
		// s = -C(p0, p1) / (w0 + w1 + alpha)
		// For a step by step derivation of this formula, see:
		// Comparison of XPBD and Projective Dynamics (Dennis Ledwon, Master Thesis, 2024)
		// https://github.com/dzenaut/thesis/blob/main/thesis.pdf
		// sections 3.1, 3.2.3, 3.3 and 4.1.5
		// But remember that in this sample we are using a hybrid PBD/XPBD approach,
		// so the formula for s is inspired by PBD but includes the compliance term alpha in
		// the denominator as in XPBD (no lambda accumulation in involed as in pure XPBD).
        float restLen = edgeLengths[i];
        float C = len - restLen;
        float s = -C / (w + alpha);

        // Apply position updates
		// The formulas for the position updates are:
		// p0 += s * w0 * dC/dp0
		// p1 -= s * w1 * dC/dp1
		// The negative sign in the second equation is because dC/dp1 = -dir.
		// These formulas make sense since we want to move p0 and p1 in opposite directions
		// along the edge direction, in order to satisfy the distance constraint asap.
		// Also, the position updates are weighted by the inverse masses w0 and w1,
		// so that a vertex with higher inverse mass (lower mass) moves more than a vertex
		// with lower inverse mass (higher mass).
        p0 = XMVectorAdd(p0, XMVectorScale(dir, s * w0));
        p1 = XMVectorAdd(p1, XMVectorScale(dir, -s * w1));
        XMStoreFloat3(&pos[id0], p0);
        XMStoreFloat3(&pos[id1], p1);
    }
}

void SoftBody::SolveVolumes(float compliance, float dt)
{
    float alpha = compliance / (dt * dt);

	// For each tetrahedron:
    for (int i = 0; i < numTets; ++i)
    {
		// Get vertex indices
        const XMINT4& tet = tetIds[i];
        int ids[4] = { tet.x, tet.y, tet.z, tet.w };
        float w = 0.0f;
        XMVECTOR grads[4];

        // The constraint for a tetrahedron is:
        //   C(p0, p1, p2, p3) = V(p0, p1, p2, p3) - restVolume
        // where V is the signed volume of the tetrahedron with vertices p0, p1, p2, p3.
        //
        // The signed volume V can be written as:
        //   V = 1/6 * det([p1 - p0, p2 - p0, p3 - p0])
        //     = 1/6 * ((p1 - p0) x (p2 - p0)) * (p3 - p0)
        //
        // The constraint is then:
        //   C(p0, p1, p2, p3) = 1/6 * ((p1 - p0) x (p2 - p0)) · (p3 - p0) - restVolume
        //
        // Remember that the gradient of the dot product a * b w.r.t. the vector a is b, and vice versa.
        // Indeed, if we rewrite it in component form:
        // a * b = ax * bx + ay * by + az * bz
        // Then, the partial derivatives of a * b w.r.t. ax, ay and az are:
        // d(a * b)/dax = d(ax * bx + ay * by + az * bz)/dax = d(ax * bx)/dax + d(ay * by)/dax + d(az * bz)/dax = bx + 0 + 0 = bx
        // d(a * b)/day = d(ax * bx + ay * by + az * bz)/day = d(ax * bx)/day + d(ay * by)/day + d(az * bz)/day = 0 + by + 0 = by
        // d(a * b)/daz = d(ax * bx + ay * by + az * bz)/daz = d(ax * bx)/daz + d(ay * by)/daz + d(az * bz)/daz = 0 + 0 + bz = bz
        //
        // Also remember that a * (b x c) = b * (c x a) = c * (a x b) (scalar triple product).
        //
        // Let’s denote:
        //   a = p1 - p0
        //   b = p2 - p0
        //   c = p3 - p0
        //
        // So that:
        //   V = 1/6 * (a x b) * c
        //
        // The constraint is then:
        //   C(p0, p1, p2, p3) = 1/6 * (a x b) · c - restVolume
        //
        // Before deriving the gradients, let's define some preliminary concepts:
        // - The derivative of the scalar triple product function f(a, b, c) = (a x b) * c, w.r.t. to a is:
        //   df/da = (b x c)
        //   Indeed, if we define a = (ax, ay, az), b = (bx, by, bz), c = (cx, cy, cz) and (a x b) = (ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx),
        //   So, f(a, b, c) = (a x b) * c = (ay * bz - az * by) * cx + (az * bx - ax * bz) * cy + (ax * by - ay * bx) * cz
        //   The derivative of f w.r.t. ax, ay and az are:
        //   df/dax = d[(ay * bz - az * by) * cx + (az * bx - ax * bz) * cy + (ax * by - ay * bx) * cz]/dax
        //          = 0 + (-bz) * cy + by * cz = (by * cz - bz * cy)
        //   df/day = d[(ay * bz - az * by) * cx + (az * bx - ax * bz) * cy + (ax * by - ay * bx) * cz]/day
        //          = bz * cx + 0 + (-bx) * cz = (bz * cx - bx * cz)
        //   df/daz = d[(ay * bz - az * by) * cx + (az * bx - ax * bz) * cy + (ax * by - ay * bx) * cz]/daz
        //          = (-by) * cx + bx * cy + 0 = (bx * cy - by * cx)
        //   So, we can write the gradient of f w.r.t. a as:
        //   df/da = (by * cz - bz * cy, bz * cx - bx * cz, bx * cy - by * cx) = (b x c)
        //   Note that the vector (by * cz - bz * cy, bz * cx - bx * cz, bx * cy - by * cx) is indeed the expansion of the cross product b x c.
        //   Similarly, the derivatives w.r.t. b and c are:
        //   df/db = (c x a)
        //   df/dc = (a x b)
        //
        // - da/dp0 = -I, db/dp0 = -I, dc/dp0 = -I
        //   da/dp1 = I,  db/dp2 = I,  dc/dp3 = I
        //   the other combinations are zero matrices.
        //   Indeed, we know that a = p1 - p0, p0 = (x0, y0, z0) and p1 = (x1, y1, z1).
        //   So, a = (ax, ay, az) = (x1 - x0, y1 - y0, z1 - z0).
        //   The derivatives of the components of a w.r.t. the components of p1 are:
        //   dax/dx1 = d(x1 - x0)/dx1 = 1
        //   day/dy1 = d(y1 - y0)/dy1 = 1
        //   daz/dz1 = d(z1 - z0)/dz1 = 1
        //   and zero for the other combinations. So that, in matrix form:
        //   da/dp1 = [dax/dx1  dax/dy1  dax/dz1, day/dx1  day/dy1  day/dz1, daz/dx1  daz/dy1  daz/dz1] = [1 0 0, 0 1 0, 0 0 1] = I (identity matrix).
        //   Indeed, remeber that the derivative of a vector u w.r.t. another vector v is the matrix:
        //   du/dv = [dux/dvx  dux/dvy  dux/dvz, duy/dvx  duy/dvy  duy/dvz, duz/dvx  duz/dvy  duz/dvz]
        //   The derivatives of the components of a w.r.t. the components of p0 are:
        //   dax/dx0 = d(x1 - x0)/dx0 = -1
        //   day/dy0 = d(y1 - y0)/dy0 = -1
        //   daz/dz0 = d(z1 - z0)/dz0 = -1
        //   and zero for the other combinations. So that, in matrix form:
        //   da/dp0 = [dax/dx0  dax/dy0  dax/dz0, day/dx0  day/dy0  day/dz0, daz/dx0  daz/dy0  daz/dz0] = [-1 0 0, 0 -1 0, 0 0 -1] = -I.
        //   Similar considerations hold for b and c.
        //
        // So, we have that the gradient of C w.r.t. a is (restVolume is constant, so its derivative is zero):
        //   dC/da = 1/6 * (b x c)
        //
        // Since C depends on a and a depends on p1, the chain rule gives:
        //  dC/dp1 = dC/da * da/dp1 = 1/6 * (b x c) * I = 1/6 * (b x c)
        // (we are deriving w.r.t. p1, so p0 in a is treated as constant;
        // similarly, b and c in C are treated as constants with respect to p1;
        // the only variables are a for C and p1 for a, so the chain rule applied works as illustrated above).
        //
        // Similar considerations hold for other derivatives. So the full set of gradients is:
        //  dC/dp1 = 1/6 * (b x c)
        //  dC/dp2 = 1/6 * (c x a)
        //  dC/dp3 = 1/6 * (a x b)
        //
        // As for dC/dp0, the derivation is a bit involved, so let's go step by step:
		// From vectorial calculus, we know that if a function f depends on a vector v = (a, b, c), we know that
		//  df = nabla_f · dv = (df/da, df/db, df/dc) · (da, db, dc) = df/da * da + df/db * db + df/dc * dc
		// And that if a, b and c depend, in turn, on another variable p0, then:
		//  da = da/dp0 * dp0
		//  db = db/dp0 * dp0
		//  dc = dc/dp0 * dp0
		// Substituting in the previous equation gives:
		//  df = nabla_f · dv = (df/da, df/db, df/dc) · (da/dp0 * dp0, db/dp0 * dp0, dc/dp0 * dp0)
		//                    = df/da * da/dp0 * dp0 + df/db * db/dp0 * dp0 + df/dc * dc/dp0 * dp0
		// At this point, remembering that a,b and c depend on p0, we can write the gradient of C w.r.t. p0 as:
		//  dC/dp0 = dC/da * da/dp0 + dC/db * db/dp0 + dC/dc * dc/dp0
		// But we know that:
		//  dC/da = 1/6 * (b x c)
		//  dC/db = 1/6 * (c x a)
		//  dC/dc = 1/6 * (a x b)
		//  and also that:
		//  da/dp0 = -I
		//  db/dp0 = -I
		//  dc/dp0 = -I
        // So, substituting all this in the equation for dC/dp0 gives:
		//  dC/dp0 = 1/6 * (b x c) * -I + 1/6 * (c x a) * -I + 1/6 * (a x b) * -I
		//         = -1/6 * ( (b x c) + (c x a) + (a x b) )
		// Expanding the cross products and summing them up gives:
		//  b x c = ((p2 - p0) x (p3 - p0))
		//        = p2 x p3 - p2 x p0 - p0 x p3 + p0 x p0
		//        = p2 x p3 - p2 x p0 - p0 x p3 + 0
		//        = p2 x p3 - p2 x p0 - p0 x p3
		//  c x a = ((p3 - p0) x (p1 - p0))
		//        = p3 x p1 - p3 x p0 - p0 x p1 + p0 x p0
		//        = p3 x p1 - p3 x p0 - p0 x p1 + 0
		//        = p3 x p1 - p3 x p0 - p0 x p1
		//  a x b = ((p1 - p0) x (p2 - p0))
		//        = p1 x p2 - p1 x p0 - p0 x p2 + p0 x p0
		//        = p1 x p2 - p1 x p0 - p0 x p2 + 0
		//        = p1 x p2 - p1 x p0 - p0 x p2
		//  (b x c) + (c x a) + (a x b) = (p2 x p3 - p2 x p0 - p0 x p3) + (p3 x p1 - p3 x p0 - p0 x p1) + (p1 x p2 - p1 x p0 - p0 x p2)
		//                              = p2 x p3 + p3 x p1 + p1 x p2
		// Substituting this result back in the equation for dC/dp0 gives:
		//  dC/dp0 = -1/6 * (p2 x p3 + p3 x p1 + p1 x p2)
		//         = 1/6 * -(p2 x p3 + p3 x p1 + p1 x p2)
		//         = 1/6 * (p3 x p2 - p1 x p3 - p2 x p1)
		//         = 1/6 * (p3 - p1) x (p2 - p1)
		// Indeed, expanding (p3 - p1) x (p2 - p1) gives:
		//  (p3 - p1) x (p2 - p1) = p3 x p2 - p3 x p1 - p1 x p2 + p1 x p1
		//                        = p3 x p2 - p3 x p1 - p1 x p2 + 0
		//                        = p3 x p2 - p3 x p1 - p1 x p2 (same as above)
		//
		// In the following code, volIdOrder defines the order of the vertices to compute the gradients
		// as explained above. So, for example, for j = 0, volIdOrder returns (1,3,2), which means that dC/dp0
		// will be computed as explained:
		//  dC/dp0 = 1/6 * (p3 - p1) x (p2 - p1)
		// The same logic applies to gradient computation w.r.t. p1, p2 and p3.
        for (int j = 0; j < 4; ++j)
        {
            int id0 = ids[volIdOrder[j][0]];
            int id1 = ids[volIdOrder[j][1]];
            int id2 = ids[volIdOrder[j][2]];

            XMVECTOR v0 = XMLoadFloat3(&pos[id0]);
            XMVECTOR v1 = XMLoadFloat3(&pos[id1]);
            XMVECTOR v2 = XMLoadFloat3(&pos[id2]);

            XMVECTOR d1 = XMVectorSubtract(v1, v0);
            XMVECTOR d2 = XMVectorSubtract(v2, v0);
            grads[j] = XMVectorScale(XMVector3Cross(d1, d2), 1.0f / 6.0f);

			// Note that here w includes the squared length of the gradient vector
            w += invMass[ids[j]] * XMVectorGetX(XMVector3LengthSq(grads[j]));
        }

        if (w == 0.0f) continue;

        // Compute constraint value and scaling
		// The formula for the scaling factor s is:
		// s = -C(p0, p1, p2, p3) / (w0 * |dC/dp0|^2 + w1 * |dC/dp1|^2 + w2 * |dC/dp2|^2 + w3 * |dC/dp3|^2 + alpha)
		// but as noted above, here w already includes the squared lengths of the gradients multiplied,
		// so the formula simplifies to:
		// s = -C(p0, p1, p2, p3) / (w0 + w1 + w2 + w3 + alpha)
		// For a step by step derivation of this formula, see:
		// Comparison of XPBD and Projective Dynamics (Dennis Ledwon, Master Thesis, 2024)
		// https://github.com/dzenaut/thesis/blob/main/thesis.pdf
		// sections 3.1, 3.2.3, 3.3 and 4.1.5
		// But remember that in this sample we are using a hybrid PBD/XPBD approach,
		// so the formula for s is inspired by PBD but includes the compliance term alpha in
		// the denominator as in XPBD (no lambda accumulation in involed as in pure XPBD). factor
        float vol = GetTetVolume(i);
        float restVol_i = restVol[i];
        float C = vol - restVol_i;
        float s = -C / (w + alpha);

        // Apply position updates for each vertex
        // The formula for the position update is:
		//  p_j += s * w_j * dC/dp_j  for j = 0,1,2,3
		// The gradient of a function indicates the direction of the steepest ascent of that function.
		// In this case, dC/dp_j indicates the direction in which, if we move the vertex p_j,
		// the constraint C will increase the most rapidly.
		// So, using the gradient to update the position of the vertex makes sense,
		// because it allows us to adjust the vertex position so that the constraint C is satisfied asap.
		// Also, the position updates are weighted by the inverse masses, so that a vertex with higher
		// inverse mass (lower mass) moves more than a vertex with lower inverse mass (higher mass).
        for (int j = 0; j < 4; ++j)
        {
            int id = ids[j];
            XMVECTOR p = XMLoadFloat3(&pos[id]);
            p = XMVectorAdd(p, XMVectorScale(grads[j], s * invMass[id]));
            XMStoreFloat3(&pos[id], p);
        }
    }
}

void SoftBody::Solve(float dt)
{
    SolveEdges(edgeCompliance, dt);
    SolveVolumes(volCompliance, dt);
}

void SoftBody::PostSolve(float dt)
{
    for (int i = 0; i < numVerts; i++)
    {
		// skip fixed vertices
        if (invMass[i] == 0.0f)
            continue;

		// Compute velocity for the next step:
		// vel = (pos - prevPos) / dt
        XMVECTOR p = XMLoadFloat3(&pos[i]);
        XMVECTOR pp = XMLoadFloat3(&prevPos[i]);
        XMVECTOR v = XMVectorScale(XMVectorSubtract(p, pp), 1.0f / dt);

        XMStoreFloat3(&vel[i], v);
    }
}

void SoftBody::StartGrab(const wi::scene::PickResult pick)
{
    // Find the closest vertex to pick.position
    XMVECTOR pickPos = XMVectorSet(pick.position.x, pick.position.y, pick.position.z, 0.0f);
    float minDistSq = std::numeric_limits<float>::max();
    grabId = -1;

    for (int i = 0; i < numVerts; i++)
    {
        XMVECTOR p = XMLoadFloat3(&pos[i]);
        float d2 = XMVectorGetX(XMVector3LengthSq(XMVectorSubtract(p, pickPos)));
        if (d2 < minDistSq) {
            minDistSq = d2;
            grabId = i;
        }
    }

    // If a valid vertex is found, fix it at pick.position
    if (grabId >= 0)
    {
        grabInvMass = invMass[grabId];
        invMass[grabId] = 0.0f;
        XMStoreFloat3(&pos[grabId], pickPos);
    }
}

void SoftBody::MoveGrabbed(const std::vector<float>& pos, const std::vector<float>& vel)
{
    if (grabId >= 0)
    {
        // Copy pos to this->pos[grabId]
        XMVECTOR p = XMVectorSet(pos[0], pos[1], pos[2], 0.0f);
        XMStoreFloat3(&this->pos[grabId], p);
    }
}

void SoftBody::EndGrab(const std::vector<float>& pos, std::vector<float>& vel)
{
    if (grabId >= 0)
    {
        // Restore inverse mass
        invMass[grabId] = grabInvMass;

        // Copy vel to this->vel[grabId]
        XMVECTOR v = XMVectorSet(vel[0], vel[1], vel[2], 0.0f);
        XMStoreFloat3(&this->vel[grabId], v);
    }
    grabId = -1;
}
