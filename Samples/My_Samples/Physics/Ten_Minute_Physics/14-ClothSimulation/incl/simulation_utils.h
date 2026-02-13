#pragma once

#include "meshes.h"
#include "cloth.h"
#include <memory>
#include <wiScene.h>
#include <wiGUI.h>

namespace simulation
{
	void init_simulation(uint64_t visShaderID, uint64_t wireShaderID, const std::string &modelPath);
    void simulate(float dt);

    void update_tetMesh(const WireMesh& tetMesh, wi::scene::MeshComponent& mesh, bool updateGPUBuffer);
    void update_visMesh(const SimulationObject& softbody, wi::scene::MeshComponent& mesh, bool updateGPUBuffer);

    std::unique_ptr<SimulationObject> create_simulation_object(const SimulationParams& params, const std::string& modelPath);
    void remove_simulation_object(SimulationObject& sbi);
	void new_simulation_object(const std::string& modelPath);

	struct PhysicsScene
	{
		XMFLOAT3 gravity = {0.0f, -10.0f, 0.0f};
		int numSubsteps = 10;
		float dt = 1.0f / 60.0f;
		bool paused = true;
		std::vector<std::unique_ptr<SimulationObject>> objects;

		inline void Run() { paused = !paused; }
	};

}

extern simulation::PhysicsScene gPhysicsScene;
