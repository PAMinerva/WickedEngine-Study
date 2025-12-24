#pragma once

#include "balls.h"
#include <memory>
#include <wiScene.h>
#include <wiGUI.h>

namespace simulation
{
    void simulate(float dt);

	struct PhysicsScene
	{
		XMFLOAT3 gravity = {0.0f, 0.0f, 0.0f};
		// int numSubsteps = 10;
		// float dt = 1.0f / 60.0f;
		bool paused = true;
		std::unique_ptr<Balls> balls;

		inline void Run() { paused = !paused; }
	};

}

extern simulation::PhysicsScene gPhysicsScene;
