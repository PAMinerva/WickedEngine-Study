#include <fstream>
#include <iterator>
#include <sstream>
#include "meshes.h"

std::string gMeshPath;

RawMesh& getMeshData()
{
	static RawMesh mesh;
	static bool initialized = false;
	if (!initialized) {
		mesh.name = "TetMesh";
		std::ifstream file(gMeshPath);
		std::string token;

		// Read until "verts"
		while (file >> token && token != "verts") {}

		// Read verts until "faceTriIds"
		while (file >> token && token != "faceTriIds") {
			std::istringstream iss(token);
			float value;
			if (iss >> value && iss.eof()) {
				mesh.verts.push_back(value);
			}
		}

		// Read faceTriIds until EOF
		while (file >> token) {
			std::istringstream iss(token);
			int idx;
			if (iss >> idx && iss.eof()) {
				mesh.faceTriIds.push_back(idx);
			}
		}

		initialized = true;
	}
	return mesh;
}
