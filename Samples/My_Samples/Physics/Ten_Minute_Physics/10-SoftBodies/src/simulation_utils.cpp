#include <memory>
#include <random>
#include "simulation_utils.h"
#include "softBody.h"

static wi::ecs::Entity mat_entity = wi::ecs::INVALID_ENTITY;
static uint64_t shader_id = 0;

simulation::PhysicsScene gPhysicsScene;
wi::gui::Label label_tets;

namespace simulation
{
	void init_physics(uint64_t shaderID)
	{
		shader_id = shaderID;

		SoftBodyParams params{};
		params.rotate = {0.0f, 0.0f, 3.14f}; // rotate 180 degrees around Z axis clockwise (180 x (pi/180) = pi = 3.14 rad)
		params.edgeCompliance = 100.0f;
		params.volCompliance = 0.0f;

		auto instance = create_softbody_instance(params);
		int numTets = instance->softBody->numTets;
		gPhysicsScene.objects.push_back(std::move(instance));

		std::string text = std::to_string(numTets) + " tetrahedrons";
		label_tets.SetText(text);
	}

	void simulate(float dt)
	{
		if (gPhysicsScene.paused)
			return;

		// float subdt = dt / float(gPhysicsScene.numSubsteps); // Variable timestep not recommended
		float subdt = gPhysicsScene.dt / float(gPhysicsScene.numSubsteps);

		for (int step = 0; step < gPhysicsScene.numSubsteps; ++step)
		{
			// 1. PreSolve: explicit integration
			for (auto &instance : gPhysicsScene.objects)
				instance->softBody->PreSolve(subdt, gPhysicsScene.gravity);

			// 2. Solve: project constraints (edge, volume, etc.)
			// Projecting a set of points according to a constraint means moving (correcting) the points such that they satisfy the constraint
			for (auto &instance : gPhysicsScene.objects)
				instance->softBody->Solve(subdt);

			// 3. PostSolve: update velocity
			for (auto &instance : gPhysicsScene.objects)
				instance->softBody->PostSolve(subdt);
		}
	}

	void new_softbody_instance()
	{
		float a = 3.0f;

		// Random number generator setup
		static std::random_device rd;
		static std::mt19937 gen(rd());

		std::uniform_real_distribution<float> dist_x(-a, a);
		std::uniform_real_distribution<float> dist_y(0, a);
		std::uniform_real_distribution<float> dist_z(0, a);

		SoftBodyParams params{};
		params.translate = {dist_x(gen), dist_y(gen), dist_z(gen)};
		params.edgeCompliance = 100.0f;
		params.volCompliance = 0.0f;

		auto instance = create_softbody_instance(params);
		gPhysicsScene.objects.push_back(std::move(instance));

		size_t numTets = 0;
		for (auto &instance : gPhysicsScene.objects)
		{
			numTets += instance->softBody->numTets;
		}
		std::string text = std::to_string(numTets) + " tetrahedrons";
		label_tets.SetText(text);
	}

	void update_mesh(const SoftBody &softBody,
					wi::scene::MeshComponent &mesh,
					bool updateGPUBuffer)
	{
		// Update the vertex positions of the mesh, triangle by triangle
		size_t numTris = bunnyMesh.tetSurfaceTriIds.size() / 3;
		for (size_t t = 0; t < numTris; ++t) // for each surface triangle...
		{
			// Indices of the vertices for the current triangle
			uint32_t i0 = bunnyMesh.tetSurfaceTriIds[t * 3 + 0];
			uint32_t i1 = bunnyMesh.tetSurfaceTriIds[t * 3 + 1];
			uint32_t i2 = bunnyMesh.tetSurfaceTriIds[t * 3 + 2];

			// New positions from physics simulation
			XMFLOAT3 p0 = softBody.pos[i0];
			XMFLOAT3 p1 = softBody.pos[i1];
			XMFLOAT3 p2 = softBody.pos[i2];

			// Compute the new triangle normal
			XMVECTOR v0 = XMLoadFloat3(&p0);
			XMVECTOR v1 = XMLoadFloat3(&p1);
			XMVECTOR v2 = XMLoadFloat3(&p2);
			XMVECTOR n = XMVector3Cross(v1 - v0, v2 - v0);
			n = XMVector3Normalize(n);

			XMFLOAT3 fn;
			XMStoreFloat3(&fn, n);

			// Update vertex positions
			size_t baseIndex = t * 3;
			mesh.vertex_positions[baseIndex + 0] = p0;
			mesh.vertex_positions[baseIndex + 1] = p1;
			mesh.vertex_positions[baseIndex + 2] = p2;

			// Update vertex normals
			// In flat shading, the three vertices of each triangle share the same normal
			mesh.vertex_normals[baseIndex + 0] = fn;
			mesh.vertex_normals[baseIndex + 1] = fn;
			mesh.vertex_normals[baseIndex + 2] = fn;
		}

		// Update AABB of the mesh
		XMFLOAT3 _min = XMFLOAT3(std::numeric_limits<float>::max(),
								std::numeric_limits<float>::max(),
								std::numeric_limits<float>::max());
		XMFLOAT3 _max = XMFLOAT3(std::numeric_limits<float>::lowest(),
								std::numeric_limits<float>::lowest(),
								std::numeric_limits<float>::lowest());
		for (size_t i = 0; i < mesh.vertex_positions.size(); ++i)
		{
			const XMFLOAT3 &pos = mesh.vertex_positions[i];
			_min = wi::math::Min(_min, pos);
			_max = wi::math::Max(_max, pos);
		}
		mesh.aabb = wi::primitive::AABB(_min, _max);

		if (!updateGPUBuffer)
			return;

		// Update the GPU buffer where are stored the vertex positions and normals of the mesh
		auto device = wi::graphics::GetDevice();
		auto cmd = device->BeginCommandList();
		if (mesh.position_format == wi::scene::MeshComponent::Vertex_POS16::FORMAT)
		{
			std::vector<wi::scene::MeshComponent::Vertex_POS16> pos16;
			pos16.reserve(mesh.vertex_positions.size());
			for (size_t i = 0; i < mesh.vertex_positions.size(); ++i)
			{
				wi::scene::MeshComponent::Vertex_POS16 vert;
				const XMFLOAT3 &pos = mesh.vertex_positions[i];
				uint8_t wind = 0xFF;
				vert.FromFULL(mesh.aabb, pos, wind);
				pos16.push_back(vert);
			}

			device->UpdateBuffer(&mesh.generalBuffer,
								pos16.data(),
								cmd,
								mesh.vb_pos_wind.size, // pos16.size() * sizeof(wi::scene::MeshComponent::Vertex_POS16),
								mesh.vb_pos_wind.offset);
		}
		else if (mesh.position_format == wi::scene::MeshComponent::Vertex_POS32::FORMAT)
		{
			std::vector<wi::scene::MeshComponent::Vertex_POS32> pos32;
			pos32.reserve(mesh.vertex_positions.size());
			for (size_t i = 0; i < mesh.vertex_positions.size(); ++i)
			{
				wi::scene::MeshComponent::Vertex_POS32 vert;
				const XMFLOAT3 &pos = mesh.vertex_positions[i];
				uint8_t wind = 0xFF;
				vert.FromFULL(pos);
				pos32.push_back(vert);
			}

			device->UpdateBuffer(&mesh.generalBuffer,
								pos32.data(),
								cmd,
								mesh.vb_pos_wind.size, // pos32.size() * sizeof(wi::scene::MeshComponent::Vertex_POS32),
								mesh.vb_pos_wind.offset);
		}
		else
		{
			assert("Unsupported vertex format!" && false);
		}

		std::vector<wi::scene::MeshComponent::Vertex_NOR> nor32;
		nor32.reserve(mesh.vertex_normals.size());
		for (size_t i = 0; i < mesh.vertex_normals.size(); ++i)
		{
			wi::scene::MeshComponent::Vertex_NOR vert;
			vert.FromFULL(mesh.vertex_normals[i]);
			nor32.push_back(vert);
		}

		device->UpdateBuffer(&mesh.generalBuffer,
							nor32.data(),
							cmd,
							mesh.vb_nor.size,
							mesh.vb_nor.offset);
	}

	std::unique_ptr<SoftBodyInstance> create_softbody_instance(const SoftBodyParams &params)
	{
		using namespace wi::scene;
		using namespace wi::ecs;

		// Create a material (shared by all soft body instances) that returns red color when enlighted
		if (mat_entity == INVALID_ENTITY)
		{
			mat_entity = CreateEntity();
			GetScene().materials.Create(mat_entity);
			MaterialComponent &material = *GetScene().materials.GetComponent(mat_entity);
			material.shaderType = MaterialComponent::SHADERTYPE_PBR;
			material.baseColor = XMFLOAT4(0.941f, 0.125f, 0.0f, 1.0f);
			// material.SetDoubleSided(true);
			material.SetCustomShaderID(shader_id);
		}

		// Create a new entity.
		// This entity (integer) will be used to create and access all components related to this soft body instance.
		Entity entity = CreateEntity();

		GetScene().layers.Create(entity);
		GetScene().transforms.Create(entity); // Transform component is required for applying translations, rotations, scalings, etc. to the object

		// A mesh component stores the mesh data (vertices, indices, subsets, etc.).
		// An object component represents a renderable object in the scene and holds the mesh component entity that describes its geometry.
        // Note that in this sample we compute the physics simulation on soft body mesh data and use it to update data in the
        // mesh component. So we cannot use multiple instance object sharing the same mesh component (we need a mesh per object).
		ObjectComponent &object = GetScene().objects.Create(entity);
		MeshComponent &mesh = GetScene().meshes.Create(entity);
		object.meshID = entity;

		// A subset specifies a range of indices to describe part or whole
		// of the mesh to be rendered with a specific material.
		mesh.subsets.push_back(wi::scene::MeshComponent::MeshSubset());
		mesh.subsets.back().materialID = mat_entity;
		mesh.subsets.back().indexOffset = 0;
		mesh.subsets.back().indexCount = (uint32_t)bunnyMesh.tetSurfaceTriIds.size(); // tetSurfaceTriIds contains indices of surface triangles only

		// Compute the number of surface triangles:
		// For rendering purposes, we are only interested in surface triangles of the tetrahedral mesh.
		size_t numTris = bunnyMesh.tetSurfaceTriIds.size() / 3;
		mesh.vertex_positions.reserve(numTris * 3);
		mesh.vertex_normals.reserve(numTris * 3);
		mesh.indices.reserve(numTris * 3);

		// For each surface triangle of the tetrahedral mesh
		for (size_t t = 0; t < numTris; ++t)
		{
			// Indices of the vertices for this triangle in clocwise order
			uint32_t i0 = bunnyMesh.tetSurfaceTriIds[t * 3 + 0];
			uint32_t i1 = bunnyMesh.tetSurfaceTriIds[t * 3 + 1];
			uint32_t i2 = bunnyMesh.tetSurfaceTriIds[t * 3 + 2];

			// Positions of the three vertices for this triangle
			XMFLOAT3 p0(bunnyMesh.verts[i0 * 3 + 0],
						bunnyMesh.verts[i0 * 3 + 1],
						bunnyMesh.verts[i0 * 3 + 2]);

			XMFLOAT3 p1(bunnyMesh.verts[i1 * 3 + 0],
						bunnyMesh.verts[i1 * 3 + 1],
						bunnyMesh.verts[i1 * 3 + 2]);

			XMFLOAT3 p2(bunnyMesh.verts[i2 * 3 + 0],
						bunnyMesh.verts[i2 * 3 + 1],
						bunnyMesh.verts[i2 * 3 + 2]);

			// Compute the new triangle normal
			XMVECTOR v0 = XMLoadFloat3(&p0);
			XMVECTOR v1 = XMLoadFloat3(&p1);
			XMVECTOR v2 = XMLoadFloat3(&p2);
			// We are in a left-handed coordinate system, so the vector resulting from the cross product
			// "sees" the first vector rotated clockwise towards the second vector.
			// In this case, the vertices of each triangle are given in clockwise order, so the first
			// and second vectors will be, respectively, (v1 - v0) and (v2 - v0).
			XMVECTOR n = XMVector3Cross(v1 - v0, v2 - v0);
			n = XMVector3Normalize(n);

			XMFLOAT3 fn;
			XMStoreFloat3(&fn, n);

			// Add the three vertices of this triangle to the vertex_position array of the mesh component.
			// Remember that we are only storing vertices of surface triangles.
			uint32_t baseIndex = (uint32_t)mesh.vertex_positions.size();
			mesh.vertex_positions.push_back(p0);
			mesh.vertex_positions.push_back(p1);
			mesh.vertex_positions.push_back(p2);

			// In flat shading, the three vertices of each triangle share the same normal
			mesh.vertex_normals.push_back(fn);
			mesh.vertex_normals.push_back(fn);
			mesh.vertex_normals.push_back(fn);

			// Add the three vertex indices of this triangle to the index array of the mesh component.
			// Note that these indices are relative to the vertex_position array.
			mesh.indices.push_back(baseIndex + 0);
			mesh.indices.push_back(baseIndex + 1);
			mesh.indices.push_back(baseIndex + 2);
		}

		// Create a new SoftBodyInstance object that will contain the SoftBody object and the entity
		// we will use to access the mesh and other components associated with this soft body object.
		// SoftBody will store and manage all physics-related data (some taken from bunnyMesh) and operations
		// (vertex positions will include all vertices of the tetrahedral mesh, not only surface vertices).
		auto instance = std::make_unique<SoftBodyInstance>();
		instance->softBody = std::make_unique<SoftBody>(bunnyMesh, params.edgeCompliance, params.volCompliance);
		instance->entity = entity;

		// Apply initial transformations to the SoftBody object according to the given parameters
		instance->softBody->Scale(params.scale);
		instance->softBody->RotateRollPitchYaw(params.rotate);
		instance->softBody->Translate(params.translate);

		// Initialize the physics-related data structures of the SoftBody object
		instance->softBody->InitPhysics();

		// Update mesh data according to the current state of the SoftBody object.
		simulation::update_mesh(*instance->softBody, mesh, false);

		// Upload the mesh data to a buffer created as a sub-resource of a general GPU buffer
		// containing all mesh data of the scene, if available.
		// Otherwise, create a new GPU buffer for this mesh only.
		mesh.CreateRenderData();

		return instance;
	}

	void remove_softbody_instance(SoftBodyInstance& sbi)
	{
		using namespace wi::scene;
		using namespace wi::ecs;

		if (mat_entity != INVALID_ENTITY)
		{	
			GetScene().materials.Remove(mat_entity);
			mat_entity = INVALID_ENTITY;
		}

		GetScene().layers.Remove(sbi.entity);
		GetScene().transforms.Remove(sbi.entity);
		GetScene().objects.Remove(sbi.entity);
		GetScene().meshes.Remove(sbi.entity);
	}
}
