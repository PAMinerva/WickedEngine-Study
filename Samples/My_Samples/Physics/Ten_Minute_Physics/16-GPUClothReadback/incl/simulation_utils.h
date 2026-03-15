#pragma once

#include "cloth.h"
#include <memory>
#include <wiScene.h>
#include <wiGUI.h>
#include <wiGraphics.h>

namespace simulation
{
    void init_simulation(uint64_t wireShaderID);
    void simulate(wi::graphics::CommandList cmd, float frameDt);

    void update_mesh(const ClothMesh& cloth, wi::scene::MeshComponent& mesh, wi::graphics::CommandList cmd);

	void create_wire_mesh(std::unique_ptr<SimulationObject>& instance, const SimulationParams& params);

    std::unique_ptr<SimulationObject> create_simulation_object(const SimulationParams& params);
    void remove_simulation_object(SimulationObject& obj);

    struct PhysicsScene
    {
        XMFLOAT3 gravity = {0.0f, -10.0f, 0.0f};
		float dt = 1.0f / 60.0f;
        int numSubsteps = 30;
        bool paused = true;
        int solveType = 0;  // 0=coloring, 1=jacobi
        std::vector<std::unique_ptr<SimulationObject>> objects;

        // Sphere collision parameters
        float sphereCenter[3] = {0.0f, 1.5f, 0.0f};
        float sphereRadius = 0.5f;

        inline void Run() { paused = !paused; }
    };
}

extern simulation::PhysicsScene gPhysicsScene;
