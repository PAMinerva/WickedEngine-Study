#pragma once
#include <string>
#include <vector>

struct TetMesh
{
    std::string name;
    std::vector<float> verts;
    std::vector<int> tetIds;
    std::vector<int> tetEdgeIds;
};

struct VisMesh
{
    std::string name;
    std::vector<float> verts;
    std::vector<int> triIds;
};


extern std::string gTetMeshPath;

TetMesh& getTetMesh();
// VisMesh& getDragonVisMesh();
