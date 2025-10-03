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
	// compute its rest volume and the inverse mass contribution to its four vertices
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

    XMVECTOR cross = XMVector3Cross(d1, d2);
    float dot = XMVectorGetX(XMVector3Dot(cross, d3));
    return dot / 6.0f;
}

void SoftBody::PreSolve(float dt, const XMFLOAT3& gravity)
{
    XMVECTOR grav = XMLoadFloat3(&gravity);

    for (int i = 0; i < numVerts; i++)
    {
        if (invMass[i] == 0.0f)
            continue;

        // Load vectors
        XMVECTOR vel_i = XMLoadFloat3(&vel[i]);
        XMVECTOR pos_i = XMLoadFloat3(&pos[i]);
        XMVECTOR prev_i = XMLoadFloat3(&prevPos[i]);

        // Explicit integration: v += g*dt, prev = pos, pos += v*dt
        vel_i = XMVectorAdd(vel_i, XMVectorScale(grav, dt));
        prev_i = pos_i;
        pos_i = XMVectorAdd(pos_i, XMVectorScale(vel_i, dt));

        // Ground collision (y < 0)
        if (XMVectorGetY(pos_i) < 0.0f) {
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
		// Compute the inverse mass sum
        const XMINT2& id = edgeIds[i];
        int id0 = id.x;
        int id1 = id.y;
        float w0 = invMass[id0];
        float w1 = invMass[id1];
        float w = w0 + w1;
        if (w == 0.0f) continue;

        // Edge positions
        XMVECTOR p0 = XMLoadFloat3(&pos[id0]);
        XMVECTOR p1 = XMLoadFloat3(&pos[id1]);

        // Current length
        XMVECTOR diff = XMVectorSubtract(p0, p1);
        float len = XMVectorGetX(XMVector3Length(diff));
        if (len == 0.0f) continue;

        // Direzione normalizzata
        XMVECTOR dir = XMVectorScale(diff, 1.0f / len);

        // Compute the constraint and its scaling factor
        float restLen = edgeLengths[i];
        float C = len - restLen;
        float s = -C / (w + alpha);

        // Apply position updates
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

        // Compute gradients for each vertex of the tetrahedron
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

            w += invMass[ids[j]] * XMVectorGetX(XMVector3LengthSq(grads[j]));
        }

        if (w == 0.0f) continue;

        // Compute constraint value and scaling factor
        float vol = GetTetVolume(i);
        float restVol_i = restVol[i];
        float C = vol - restVol_i;
        float s = -C / (w + alpha);

        // Apply position updates for each vertex
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
        if (invMass[i] == 0.0f)
            continue;

        // Compute velocity: vel = (pos - prevPos) / dt
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
