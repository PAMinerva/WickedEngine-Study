#pragma once

#include "meshes.h"
#include "softBody.h"
#include <memory>
#include <wiScene.h>
#include <wiGUI.h>

namespace simulation
{
	void init_simulation(/*uint64_t visShaderID,*/ uint64_t wireShaderID, const std::string& modelPath);
    void simulate(float dt);

    void update_tetMesh(const TetraMesh& tetMesh, wi::scene::MeshComponent& mesh, bool updateGPUBuffer);
    void update_visMesh(const SoftBodySkinning& softbody, wi::scene::MeshComponent& mesh, bool updateGPUBuffer);

    std::unique_ptr<SoftBodySkinning> create_softbody_skinning(const SoftBodySkinningParams& params, const std::string& modelPath);
    void remove_softbody_skinning(SoftBodySkinning& sbi);
	void new_softbody_skinning(const std::string& modelPath);

	struct PhysicsScene
	{
		XMFLOAT3 gravity = {0.0f, -10.0f, 0.0f};
		int numSubsteps = 10;
		float dt = 1.0f / 60.0f;
		bool paused = true;
		std::vector<std::unique_ptr<SoftBodySkinning>> objects;

		inline void Run() { paused = !paused; }

		void Squash()
		{
			for (auto& instance : objects)
				instance->tetMesh->Squash();
		}
	};

}

extern simulation::PhysicsScene gPhysicsScene;
