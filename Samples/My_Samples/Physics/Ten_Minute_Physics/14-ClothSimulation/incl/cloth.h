#pragma once
#include <memory>
#include <vector>
#include <array>
#include <cstddef>
#include <wiScene.h>
#include "meshes.h"

struct VisMesh
{
    VisMesh(const RawMesh &rawMesh);

    std::vector<XMFLOAT3> verts;
	std::vector<XMINT3> triIds;
	std::vector<uint32_t> indices;
    uint32_t numVerts = 0;
    uint32_t numTri = 0;
	uint32_t numIndices = 0;
};

struct Edge {
    int id0;
    int id1;
    int edgeNr;
};

class WireMesh {

public:
    WireMesh(
		const RawMesh &rawMesh,
        float stretchingCompliance = 0.0f,
        float bendingCompliance = 1.0f
    );

    void InitPhysics();
	void PreSolve(float dt, const XMFLOAT3& gravity);
    void Solve(float dt);
    void PostSolve(float dt);

	std::vector<float> FindTriNeighbors(const std::vector<uint32_t>& triIds);

	void SolveStretching(float dt);
	void SolveBending(float dt);

	void Translate(const XMFLOAT3& delta);
	void RotateRollPitchYaw(const XMFLOAT3& angle);
	void Scale(const XMFLOAT3& scale);

	void StartGrab(const wi::scene::PickResult pick);
    void MoveGrabbed(const std::vector<float>& pos, const std::vector<float>& vel);
    void EndGrab(const std::vector<float>& pos, std::vector<float>& vel);

    int numVerts = 0;
    int numTris = 0;

    std::vector<XMFLOAT3> pos;
    std::vector<XMFLOAT3> prevPos;
	std::vector<XMFLOAT3> restPos;
    std::vector<XMFLOAT3> vel;

    std::vector<float> invMass;

    std::vector<XMINT4> bendingIds;
    std::vector<XMINT2> stretchingIds;
    std::vector<float> stretchingLengths;
	std::vector<float> bendingLengths;

    std::vector<uint32_t> faceTriIds;
	std::vector<uint16_t> edgeIds;

    float stretchingCompliance = 0.0f;
    float bendingCompliance = 1.0f;

    float grabInvMass = 0.0f;
	int grabId = -1;
};

struct SimulationParams
{
    XMFLOAT3 scale = {1.0f, 1.0f, 1.0f};
    XMFLOAT3 rotate = {0.0f, 0.0f, 0.0f};
    XMFLOAT3 translate = {0.0f, 0.0f, 0.0f};
    float stretchingCompliance = 0.0f;
    float bendingCompliance = 1.0f;
};

struct SimulationObject {
	std::unique_ptr<WireMesh> wireMesh;
	std::unique_ptr<VisMesh> visMesh;
	uint64_t wireEntity;
	uint64_t visEntity;
};


