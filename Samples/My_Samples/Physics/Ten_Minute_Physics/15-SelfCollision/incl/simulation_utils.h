#pragma once

#include "cloth.h"
#include <memory>
#include <wiScene.h>
#include <wiGUI.h>

namespace simulation
{
	void init_simulation(uint64_t frontShaderID, uint64_t backShaderID, uint64_t wireShaderID);
    void simulate(float dt);

    void update_mesh(const ClothMesh& tetMesh, wi::scene::MeshComponent& mesh, bool updateGPUBuffer);

    std::unique_ptr<SimulationObject> create_simulation_object(const SimulationParams& params);
    void remove_simulation_object(SimulationObject& sbi);
	void new_simulation_object();

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
