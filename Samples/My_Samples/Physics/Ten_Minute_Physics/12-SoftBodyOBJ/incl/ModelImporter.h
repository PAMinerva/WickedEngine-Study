#pragma once

std::pair<wi::ecs::Entity, wi::ecs::Entity> ImportModel_OBJ(const std::string& fileName, wi::scene::Scene& scene, bool transform_to_LH = false);
