#include <memory>
#include <random>
#include "simulation_utils.h"
#include "balls.h"

simulation::PhysicsScene gPhysicsScene;

namespace simulation
{
	void simulate(float dt)
	{
		if (gPhysicsScene.paused)
			return;

		gPhysicsScene.balls->Simulate(dt, gPhysicsScene.gravity);
	}
}
