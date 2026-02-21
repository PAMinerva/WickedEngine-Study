#pragma once
#include <memory>
#include <vector>
#include <cstdint>
#include <wiScene.h>
#include "hashing.h"

struct SimulationParams
{
    int numX = 30;
    int numY = 200;
    float spacing = 0.01f;
    float thickness = 0.01f;
    float stretchingCompliance = 0.0f;
    float shearCompliance = 0.0001f;
    float bendingCompliance = 1.0f;
    XMFLOAT3 translate = {0.0f, 0.0f, 0.0f};
};

class ClothMesh
{
public:
    ClothMesh(const SimulationParams& params);

    // Main simulation method (replaces PreSolve/Solve/PostSolve)
    void Simulate(float frameDt, int numSubSteps, const XMFLOAT3& gravity);

    void SolveConstraints(float dt);
    void SolveGroundCollisions();
    void SolveCollisions(float dt);

    // Update bending compliance at runtime
    void UpdateBendingCompliance(float newCompliance);

    void StartGrab(const wi::scene::PickResult& pick);
    void MoveGrabbed(const std::vector<float>& pos, const std::vector<float>& vel);
    void EndGrab(const std::vector<float>& pos, std::vector<float>& vel);

    // Particle data
    int numParticles = 0;
    int numX = 0;
    int numY = 0;
    float thickness = 0.01f;
    bool handleCollisions = true;

    std::vector<XMFLOAT3> pos;
    std::vector<XMFLOAT3> prevPos;
    std::vector<XMFLOAT3> restPos;
    std::vector<XMFLOAT3> vel;
    std::vector<float> invMass;

    // Constraint data
    int numConstraints = 0;
	int numBendingConstraints = 0;
    std::vector<XMINT2> constraintIds;      // Pairs of particle indices
    std::vector<float> restLens;            // Rest lengths
    std::vector<float> compliances;         // Per-constraint compliance

    // Compliance values (stored for runtime updates)
    float stretchingCompliance = 0.0f;
    float shearCompliance = 0.0001f;
	float bendingCompliance = 1.0f;

    // Triangle indices for rendering
    std::vector<uint32_t> triIds;
    std::vector<uint32_t> edgeIds;          // For wireframe visualization
    int numTris = 0;

    // Spatial hashing for self-collision
    std::unique_ptr<SpatialHashing::Hash> hash;

    // Temporary vectors for calculations (avoid allocations)
    std::vector<XMFLOAT3> vecs;

    // Grabbing
    float grabInvMass = 0.0f;
    int grabId = -1;
};

// Each SimulationObject contains a Cloth and associated entity IDs
struct SimulationObject
{
    std::unique_ptr<ClothMesh> cloth;
    uint64_t frontEntity;   // Front face mesh (red)
    uint64_t backEntity;    // Back face mesh (yellow)
	uint64_t wireEntity;
};
