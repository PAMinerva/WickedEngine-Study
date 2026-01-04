#pragma once
#include <memory>
#include <vector>
#include <array>
#include <cstddef>
#include <wiScene.h>
#include "meshes.h"

struct VisualMesh
{
    VisualMesh(const VisMesh &visMesh);

    std::vector<XMFLOAT3> verts;
	std::vector<XMINT3> triIds;
    int numVerts = 0;
    int numTris = 0;
};

class TetraMesh {
public:
    TetraMesh(
		const TetMesh &tetMesh,
        float edgeCompliance = 100.0f,
        float volCompliance = 0.0f
    );

    void InitPhysics();
	void PreSolve(float dt, const XMFLOAT3& gravity);
    void Solve(float dt);
    void PostSolve(float dt);

	float GetTetVolume(int nr);
	void SolveEdges(float compliance, float dt);
	void SolveVolumes(float compliance, float dt);

    void Squash();
	void Translate(const XMFLOAT3& delta);
	void RotateRollPitchYaw(const XMFLOAT3& angle);
	void Scale(const XMFLOAT3& scale);

	void StartGrab(const wi::scene::PickResult pick);
    void MoveGrabbed(const std::vector<float>& pos, const std::vector<float>& vel);
    void EndGrab(const std::vector<float>& pos, std::vector<float>& vel);

    int numVerts = 0;
    int numTets = 0;

	// Tetrahedron index stored as first component,
	// followed by the first three barycentric coordinates
	// (the fourth can be derived as 1 - (b0 + b1 + b2))
	std::vector<XMFLOAT4> skinningInfo;

    std::vector<XMFLOAT3> pos;
    std::vector<XMFLOAT3> prevPos;
    std::vector<XMFLOAT3> vel;

    std::vector<XMINT4> tetIds;
    std::vector<XMINT2> edgeIds;
    std::vector<XMINT3> surfaceTriIds;

    std::vector<float> restVol;
    std::vector<float> edgeLengths;
    std::vector<float> invMass;

    float edgeCompliance = 100.0f;
    float volCompliance = 0.0f;

    std::vector<float> grads;

    float grabInvMass = 0.0f;
	int grabId = -1;

    std::array<std::array<int,3>,4> volIdOrder = {{
        {1,3,2}, {0,2,3}, {0,3,1}, {0,1,2}
    }};
};

struct SoftBodySkinningParams
{
    XMFLOAT3 scale = {1.0f, 1.0f, 1.0f};
    XMFLOAT3 rotate = {0.0f, 0.0f, 0.0f};
    XMFLOAT3 translate = {0.0f, 0.0f, 0.0f};
    float edgeCompliance = 100.0f;
    float volCompliance = 0.0f;
};

// Each SoftBodySkinning instance contains a TetraMesh (for physics simulation)
// and a VisualMesh (for rendering), along with their associated entity IDs.
struct SoftBodySkinning {
	std::unique_ptr<TetraMesh> tetMesh;
	std::unique_ptr<VisualMesh> visMesh;
	uint64_t tetEntity;
	uint64_t visEntity;
};


