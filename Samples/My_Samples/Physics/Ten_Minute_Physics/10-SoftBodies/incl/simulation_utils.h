#pragma once

#include "softBody.h"
#include <memory>
#include <wiScene.h>
#include <wiGUI.h>

namespace simulation
{
	void init_physics();
    void simulate(float dt);

    void update_mesh(const SoftBody& softbody, wi::scene::MeshComponent& mesh, bool updateGPUBuffer);

    std::unique_ptr<SoftBodyInstance> create_softbody_instance(const SoftBodyParams& params);
    void remove_softbody_instance(SoftBodyInstance& sbi);
	void new_softbody_instance();

	struct PhysicsScene
	{
		XMFLOAT3 gravity = {0.0f, -10.0f, 0.0f};
		int numSubsteps = 10;
		float dt = 1.0f / 60.0f;
		bool paused = true;
		std::vector<std::unique_ptr<SoftBodyInstance>> objects;

		inline void Run() { paused = !paused; }

		void Squash()
		{
			for (auto& instance : objects)
				instance->softBody->Squash();
		}
	};

}

extern simulation::PhysicsScene gPhysicsScene;
