#pragma once

#include <memory>
#include <vector>
#include <array>
#include <cstddef>
#include <wiECS.h>
#include <wiScene.h>
#include "hashing.h"

class Balls {
	public:
		Balls(float radius,
			  const std::vector<XMFLOAT3>& positions,
			  const std::vector<XMFLOAT3>& velocities,
			  const std::array<XMFLOAT3, 2>& worldBounds);

		void CreateBalls();
		void CreateMaterials();

		void UpdateSpherePosition(int32_t i);
		void SetSphereMaterial(size_t sphereIndex, bool inCollision);

		void Simulate(float dt, const XMFLOAT3& gravity);

		float radius;
		int numBalls;
		bool showCollisions = false;
		wi:: ecs::Entity white_mat_entity = wi::ecs::INVALID_ENTITY;
		wi::ecs::Entity mesh_entity = wi::ecs::INVALID_ENTITY;
		std::vector<XMFLOAT3> pos;
		std::vector<XMFLOAT3> prevPos;
		std::vector<XMFLOAT3> vel;
		std::vector<uint32_t> sphereEntities;
		std::array<XMFLOAT3, 2> worldBounds;

		std::unique_ptr<SpatialHashing::Hash> hash;
};
