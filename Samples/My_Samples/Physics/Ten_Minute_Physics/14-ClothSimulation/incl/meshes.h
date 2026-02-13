#pragma once
#include <cstdint>
#include <string>
#include <vector>

struct RawMesh
{
    std::string name;
    std::vector<float> verts;
    std::vector<uint32_t> faceTriIds;
};

extern std::string gMeshPath;

RawMesh& getMeshData();
