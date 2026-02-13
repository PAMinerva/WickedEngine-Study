#include <Utility/DirectXMath/DirectXMath.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <wiScene.h>
#include "cloth.h"
#include "meshes.h"
#include "simulation_utils.h"

VisMesh::VisMesh(const RawMesh& rawMesh)
	: indices(rawMesh.faceTriIds)
	, numVerts(rawMesh.verts.size() / 3)
	, numTri(rawMesh.faceTriIds.size() / 3)
	, numIndices(rawMesh.faceTriIds.size())
{
	numVerts = static_cast<int>(rawMesh.verts.size() / 3);
	numTri = static_cast<int>(rawMesh.faceTriIds.size() / 3);
	verts.reserve(numVerts);
	for (int i = 0; i < numVerts; ++i) {
		XMFLOAT3 v(
			rawMesh.verts[3 * i + 0],
			rawMesh.verts[3 * i + 1],
			rawMesh.verts[3 * i + 2]
		);
		verts.push_back(v);
	}
	triIds.reserve(numTri);
	for (int i = 0; i < numTri; ++i) {
		XMINT3 t(
			rawMesh.faceTriIds[3 * i + 0],
			rawMesh.faceTriIds[3 * i + 1],
			rawMesh.faceTriIds[3 * i + 2]
		);
		triIds.push_back(t);
	}
}

WireMesh::WireMesh(
    const RawMesh& rawMesh,
    float stretchingCompliance_,
    float bendingCompliance_
)
    : numVerts(static_cast<int>(rawMesh.verts.size() / 3))
	, numTris(static_cast<int>(rawMesh.faceTriIds.size() / 3))
    , invMass(numVerts, 0.0f)
	, faceTriIds(rawMesh.faceTriIds)
    , stretchingCompliance(stretchingCompliance_)
    , bendingCompliance(bendingCompliance_)
{
    pos.reserve(numVerts);
    prevPos.reserve(numVerts);
	restPos.reserve(numVerts);
    for (int i = 0; i < numVerts; ++i) {
        XMFLOAT3 v(
            rawMesh.verts[3 * i + 0],
            rawMesh.verts[3 * i + 1],
            rawMesh.verts[3 * i + 2]
        );
        pos.push_back(v);
        prevPos.push_back(v);
		restPos.push_back(v);
    }

	std::vector<float> neighbors = FindTriNeighbors(rawMesh.faceTriIds);
	std::vector<uint32_t> triPairIds;

	// For each triangle...
	for (uint32_t i = 0; i < numTris; i++) {

		// For each edge of the triangle (3 edges per triangle)...
		for (uint32_t j = 0; j < 3; j++) {

			// Get the vertex ids of the current edge (id0, id1)
			uint32_t id0 = rawMesh.faceTriIds[3 * i + j];
			uint32_t id1 = rawMesh.faceTriIds[3 * i + (j + 1) % 3];

			// Get the neighbor edgeNr value from neighbors.
			// This value is the index of the triangle that shares the current edge (i, j) with triangle i, or -1 if the edge is a boundary edge.
			int n = neighbors[3 * i + j];

			// Take each edge only once by storing their vertex ids in edgeIds.
			// If n < 0, it means that the edge is a boundary edge and it belongs to only one triangle,
			// so we can take it without any condition.
			// If n >= 0, it means that the edge is shared by two triangles (i and n), so we take it only if id0 < id1 to avoid duplicates.
			if (n < 0 || id0 < id1) {
				edgeIds.push_back(id0);
				edgeIds.push_back(id1);
			}

			// If n >= 0, it means that the edge is shared by two triangles (i and n),
			// so we can create a bending constraint between the two triangles.
			if (n >= 0) {
				// n is the index of the neighbor triangle that shares the current edge (id0, id1) with triangle i,
				// so we can use it to get the vertex ids of the opposite vertices of the neighbor triangle (id2 and id3).
				// It stores (id0, id1, id2, id3) so that each four elements represent the indices of the two triangles sharing
				// the same edge (id0, id1) and their opposite vertices (id2 and id3).
				// Remember that in FindTriNeighbors, edgeNr = triId * 3 + localEdgeId
				// so we can get the neighbor triangle index (ni) and the local edge index (nj) by dividing and taking the modulo of n by 3, respectively.
				uint32_t ni = n / 3;
				uint32_t nj = n % 3;
				uint32_t id2 = rawMesh.faceTriIds[3 * i + (j + 2) % 3];   // Opposite vertex of triangle i
				uint32_t id3 = rawMesh.faceTriIds[3 * ni + (nj + 2) % 3]; // Opposite vertex of triangle ni (neighbor triangle)
				// Store the vertex ids of the shared edge (id0, id1) and the opposite vertices (id2 and id3) of
				// the two triangles (i and ni) in triPairIds.
				triPairIds.push_back(id0);
				triPairIds.push_back(id1);
				triPairIds.push_back(id2);
				triPairIds.push_back(id3);
			}
		}
	}

	// For each unique edge, store the vertex ids of its endpoints (id0, id1) in stretchingIds.
	// This will be used to apply stretching (distance) constraints to edges.
	stretchingIds.reserve(edgeIds.size() / 2);
	for (size_t i = 0; i < edgeIds.size() / 2; ++i) {
		XMINT2 s(
			edgeIds[2 * i + 0],
			edgeIds[2 * i + 1]
		);
		stretchingIds.push_back(s);
	}

	// For each pair of triangles sharing an edge, store the vertex ids of the shared edge (id0, id1)
	// and the opposite vertices (id2 and id3) of the two triangles in bendingIds.
	// This will be used to apply bending (distance) constraints to pairs of triangles sharing an edge.
	bendingIds.reserve(triPairIds.size() / 4);
	for (int i = 0; i < triPairIds.size() / 4; ++i) {
		XMINT4 t(
			triPairIds[4 * i + 0],
			triPairIds[4 * i + 1],
			triPairIds[4 * i + 2],
			triPairIds[4 * i + 3]
		);
		bendingIds.push_back(t);
	}

    stretchingLengths.resize(stretchingIds.size());
    bendingLengths.resize(bendingIds.size());

	vel.resize(3 * numVerts, XMFLOAT3(0,0,0));
}

void WireMesh::InitPhysics()
{
	// Initialize inverse mass of vertices
    std::fill(invMass.begin(), invMass.end(), 0.0f);

   uint32_t numTris = faceTriIds.size() / 3;

    // For each triangle, compute its area and distribute the inverse mass to its vertices.
	// Remember that the cross product of two edges of a triangle gives a vector whose length is twice the area of the triangle.
    for (uint32_t i = 0; i < numTris; ++i)
    {
        uint32_t id0 = faceTriIds[3 * i + 0];
        uint32_t id1 = faceTriIds[3 * i + 1];
        uint32_t id2 = faceTriIds[3 * i + 2];

        XMVECTOR v0 = XMLoadFloat3(&pos[id0]);
        XMVECTOR v1 = XMLoadFloat3(&pos[id1]);
        XMVECTOR v2 = XMLoadFloat3(&pos[id2]);
        XMVECTOR e0 = XMVectorSubtract(v1, v0);
        XMVECTOR e1 = XMVectorSubtract(v2, v0);
        XMVECTOR c = XMVector3Cross(e0, e1);
        float A = 0.5f * XMVectorGetX(XMVector3Length(c));
        float pInvMass = (A > 0.0f) ? (1.0f / (A * 3.0f)) : 0.0f; // 1 / (A * 3) == (1 / A) / 3 -> area-based mass distribution, each vertex gets 1/3 of the triangle's mass
        invMass[id0] += pInvMass;
        invMass[id1] += pInvMass;
        invMass[id2] += pInvMass;
    }

	// Compute the rest length of each edge in the mesh
    for (size_t i = 0; i < stretchingLengths.size(); ++i)
    {
        uint32_t id0 = stretchingIds[i].x;
        uint32_t id1 = stretchingIds[i].y;
        XMVECTOR va = XMLoadFloat3(&pos[id0]);
        XMVECTOR vb = XMLoadFloat3(&pos[id1]);
        stretchingLengths[i] = std::sqrt(XMVectorGetX(XMVector3LengthSq(XMVectorSubtract(va, vb))));
    }

	// Compute the rest length between opposite vertices of pairs of triangles sharing an edge
	// Remember that vertex ids of pairs of triangles sharing an edge are stored in bendingIds as (id0, id1, id2, id3)
	// where (id0, id1) are the vertex ids of the shared edge and (id2, id3) are the vertex ids of the opposite vertices
	// of the two triangles.
    for (size_t i = 0; i < bendingLengths.size(); ++i)
    {
        uint32_t id0 = bendingIds[i].z;
        uint32_t id1 = bendingIds[i].w;
        XMVECTOR va = XMLoadFloat3(&pos[id0]);
        XMVECTOR vb = XMLoadFloat3(&pos[id1]);
        bendingLengths[i] = std::sqrt(XMVectorGetX(XMVector3LengthSq(XMVectorSubtract(va, vb))));
    }

	// Set fixed vertices (invMass = 0) for the top-left and top-right corners of the mesh to keep it hanging.
    float minX = std::numeric_limits<float>::max();
    float maxX = -std::numeric_limits<float>::max();
    float maxY = -std::numeric_limits<float>::max();

    for (int i = 0; i < numVerts; ++i)
    {
        minX = std::min(minX, pos[i].x);
        maxX = std::max(maxX, pos[i].x);
        maxY = std::max(maxY, pos[i].y);
    }
    float eps = 0.0001f;

    for (int i = 0; i < numVerts; ++i)
    {
        float x = pos[i].x;
        float y = pos[i].y;
        if ((y > maxY - eps) && (x < minX + eps || x > maxX - eps))
            invMass[i] = 0.0f;
    }
}

void WireMesh::PreSolve(float dt, const XMFLOAT3& gravity)
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

void WireMesh::Solve(float dt)
{
    SolveStretching(dt);
    SolveBending(dt);
}

void WireMesh::PostSolve(float dt)
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

void WireMesh::SolveStretching(float dt)
{
    float alpha = stretchingCompliance / (dt * dt);

	// For each unique edge:
    for (size_t i = 0; i < stretchingLengths.size(); ++i)
    {
		// Get the vertex ids of the edge's endpoints (id0, id1)
        const XMINT2& id = stretchingIds[i];
        int id0 = id.x;
        int id1 = id.y;

		// Compute the total inverse mass of the edge's endpoints
        float w0 = invMass[id0];
        float w1 = invMass[id1];
        float w = w0 + w1;
        if (w == 0.0f)
            continue;

        // Current length of the edge
        XMVECTOR p0 = XMLoadFloat3(&pos[id0]);
        XMVECTOR p1 = XMLoadFloat3(&pos[id1]);
        XMVECTOR diff = XMVectorSubtract(p0, p1);
        float len = XMVectorGetX(XMVector3Length(diff));
        if (len == 0.0f)
            continue;

        // Normalize the difference vector to get the direction of the constraint force
        XMVECTOR dir = XMVectorScale(diff, 1.0f / len);

		// Retrieve the rest length of the edge and compute the constraint violation (C) as the
		// difference between the current length and the rest length.
		// Compute the scaling factor (s) for the constraint correction to apply to the edge's endpoints:
		// See 10-SoftBodies sample for further explanation.
        float restLen = stretchingLengths[i];
        float C = len - restLen;
        float s = -C / (w + alpha);

        // Apply the position correction to the edge's endpoints:
        p0 = XMVectorAdd(p0, XMVectorScale(dir, s * w0));
        p1 = XMVectorAdd(p1, XMVectorScale(dir, -s * w1));

		// Store the corrected positions back to pos.
        XMStoreFloat3(&pos[id0], p0);
        XMStoreFloat3(&pos[id1], p1);
    }
}

void WireMesh::SolveBending(float dt)
{
    float alpha = bendingCompliance / (dt * dt);

    for (size_t i = 0; i < bendingLengths.size(); ++i)
    {
        const XMINT4& id = bendingIds[i];
        int id0 = id.z;  // terzo vertice del primo triangolo
        int id1 = id.w;  // terzo vertice del secondo triangolo

        float w0 = invMass[id0];
        float w1 = invMass[id1];
        float w = w0 + w1;
        if (w == 0.0f)
            continue;

        XMVECTOR p0 = XMLoadFloat3(&pos[id0]);
        XMVECTOR p1 = XMLoadFloat3(&pos[id1]);
        XMVECTOR diff = XMVectorSubtract(p0, p1);
        float len = XMVectorGetX(XMVector3Length(diff));
        if (len == 0.0f)
            continue;

        XMVECTOR dir = XMVectorScale(diff, 1.0f / len);

        float restLen = bendingLengths[i];
        float C = len - restLen;
        float s = -C / (w + alpha);

        p0 = XMVectorAdd(p0, XMVectorScale(dir, s * w0));
        p1 = XMVectorAdd(p1, XMVectorScale(dir, -s * w1));

        XMStoreFloat3(&pos[id0], p0);
        XMStoreFloat3(&pos[id1], p1);
    }
}

std::vector<float> WireMesh::FindTriNeighbors(const std::vector<uint32_t>& triIds)
{
    std::vector<Edge> edges;
	uint32_t numTris = triIds.size() / 3;
	//
    // create edge list
	//
	// For each triangle, create 3 edges (id0, id1, edgeNr) where id0 and id1 are the vertex indices of a
	// triangle edge and edgeNr is a unique number for each edge in the mesh
	// (e.g., edgeNr = 3 * triangleIndex + localEdgeIndex; where localEdgeIndex is 0, 1, or 2 for the three edges of a triangle).
	// It orders the vertex indices in each edge such that id0 < id1 to ensure that each common edge from two different triangles will have the same (id0, id1) pair.
    for (uint32_t i = 0; i < numTris; i++) {
        for (uint32_t j = 0; j < 3; j++) {
            uint32_t id0 = triIds[3 * i + j];
            uint32_t id1 = triIds[3 * i + (j + 1) % 3];
            edges.push_back({
                std::min(id0, id1),
                std::max(id0, id1),
                3 * i + j
            });
        }
    }

    // sort so that common edges are next to each other
    std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
        return (a.id0 < b.id0) || (a.id0 == b.id0 && a.id1 < b.id1);
    });

    // find matching edges
	// Initialize a vector of size 3 * numTris with -1.0f, where each element:
	// - corresponds to an edge of a triangle(3 edges per triangle, hence 3 * numTris)
	// - will store the edgeNr of the neighboring triangle.
	// If an edge has no neighbor (i.e., it's a boundary edge), it will remain -1.0f.
    std::vector<float> neighbors(3 * numTris, -1.0f);

	// When two consecutive edges in the sorted edge list share the same (id0, id1),
	// it indicates that they are shared between two triangles. In this case, we update
	// the neighbors vector by storing each edge's edgeNr in the corresponding position
	// of the other, marking them as neighbors.
    size_t nr = 0;
    while (nr < edges.size()) {
        Edge e0 = edges[nr];
        nr++;
        if (nr < edges.size()) {
            Edge e1 = edges[nr];
            if (e0.id0 == e1.id0 && e0.id1 == e1.id1) {
                neighbors[e0.edgeNr] = static_cast<float>(e1.edgeNr);
                neighbors[e1.edgeNr] = static_cast<float>(e0.edgeNr);
            }
            nr++;
        }
    }

    return neighbors;
}

void WireMesh::Translate(const XMFLOAT3& delta)
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

void WireMesh::Scale(const XMFLOAT3& scale)
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

void WireMesh::RotateRollPitchYaw(const XMFLOAT3& angle)
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

void WireMesh::StartGrab(const wi::scene::PickResult pick)
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

void WireMesh::MoveGrabbed(const std::vector<float>& pos, const std::vector<float>& vel)
{
    if (grabId >= 0)
    {
        // Copy pos to this->pos[grabId]
        XMVECTOR p = XMVectorSet(pos[0], pos[1], pos[2], 0.0f);
        XMStoreFloat3(&this->pos[grabId], p);
    }
}

void WireMesh::EndGrab(const std::vector<float>& pos, std::vector<float>& vel)
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
