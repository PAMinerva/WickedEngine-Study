#include <fstream>
#include <iterator>
#include <sstream>
#include "meshes.h"

VisMesh& getDragonVisMesh()
{
    static VisMesh mesh;
    static bool initialized = false;
    if (!initialized) {
        mesh.name = "dragonVisMesh";
        std::ifstream file("models/dragon_vis.txt");
        std::string token;

        // Read until "verts"
        while (file >> token && token != "verts") {}

        // Read verts until "indices"
        while (file >> token && token != "indices") {
            std::istringstream iss(token);
            float value;
            if (iss >> value && iss.eof()) {
                mesh.verts.push_back(value);
            }
        }

        // Read indices until EOF
        while (file >> token) {
            std::istringstream iss(token);
            int idx;
            if (iss >> idx && iss.eof()) {
                mesh.triIds.push_back(idx);
            }
        }

        initialized = true;
    }
    return mesh;
}

TetMesh& getDragonTetMesh()
{
	static TetMesh mesh;
	static bool initialized = false;
	if (!initialized) {
		mesh.name = "dragonTetMesh";
		std::ifstream file("models/dragon_tet.txt");
		std::string token;

		// Read until "verts"
		while (file >> token && token != "verts") {}

		// Read verts until "tetIds"
		while (file >> token && token != "tetIds") {
			std::istringstream iss(token);
			float value;
			if (iss >> value && iss.eof()) {
				mesh.verts.push_back(value);
			}
		}

		// Read tetIds until "tetEdgeIds"
		while (file >> token && token != "tetEdgeIds") {
			std::istringstream iss(token);
			int idx;
			if (iss >> idx && iss.eof()) {
				mesh.tetIds.push_back(idx);
			}
		}

		// Read tetEdgeIds until EOF
		while (file >> token) {
			std::istringstream iss(token);
			int idx;
			if (iss >> idx && iss.eof()) {
				mesh.tetEdgeIds.push_back(idx);
			}
		}

		initialized = true;
	}
	return mesh;
}
