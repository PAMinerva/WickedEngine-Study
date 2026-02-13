#include "simulation_utils.h"
#include "meshes.h"
#include "cloth.h"
#include <memory>
#include <random>
#include <string>
#include <wiScene_Components.h>

static wi::ecs::Entity visual_mat_entity = wi::ecs::INVALID_ENTITY;
static wi::ecs::Entity wire_mat_entity = wi::ecs::INVALID_ENTITY;
static uint64_t vis_shader_id = 0;
static uint64_t wire_shader_id = 0;

simulation::PhysicsScene gPhysicsScene;
// wi::gui::Label label_tets;
wi::gui::Label label_tris;
wi::gui::Label label_verts;

namespace simulation
{
void init_simulation(uint64_t visShaderID, uint64_t wireShaderID, const std::string &modelPath)
{
    vis_shader_id = visShaderID;
    wire_shader_id = wireShaderID;

    SimulationParams params{};

    auto instance = create_simulation_object(params, modelPath);
    int numTris = instance->wireMesh->numTris;
    int numVerts = instance->visMesh->numVerts;

    gPhysicsScene.objects.push_back(std::move(instance));

    std::string text = std::to_string(numTris) + " triangles";
    label_tris.SetText(text);

    text = std::to_string(numVerts) + " vertices";
    label_verts.SetText(text);
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
            instance->wireMesh->PreSolve(subdt, gPhysicsScene.gravity);

        // 2. Solve: project constraints (edge, volume, etc.)
        // Projecting a set of points according to a constraint means moving
        // (correcting) the points such that they satisfy the constraint
        for (auto &instance : gPhysicsScene.objects)
            instance->wireMesh->Solve(subdt);

        // 3. PostSolve: update velocity
        for (auto &instance : gPhysicsScene.objects)
            instance->wireMesh->PostSolve(subdt);
    }
}

void new_simulation_object(const std::string &modelPath)
{
    float a = 3.0f;

    // Random number generator setup
    static std::random_device rd;
    static std::mt19937 gen(rd());

    std::uniform_real_distribution<float> dist_x(-a, a);
    std::uniform_real_distribution<float> dist_y(0, a);
    std::uniform_real_distribution<float> dist_z(0, a);

    SimulationParams params{};
    params.translate = {dist_x(gen), dist_y(gen), dist_z(gen)};

    auto instance = create_simulation_object(params, modelPath);
    gPhysicsScene.objects.push_back(std::move(instance));

    size_t numTri = 0;
	size_t numVerts = 0;
	for (auto &instance : gPhysicsScene.objects)
	{
		numVerts += instance->visMesh->numVerts;
		numTri += instance->visMesh->numTri;
	}

    std::string text = std::to_string(numTri) + " triangles";
    label_tris.SetText(text);

	text = std::to_string(numVerts) + " vertices";
	label_verts.SetText(text);
}

void update_visMesh(const SimulationObject &softbody,
                    wi::scene::MeshComponent &mesh, bool updateGPUBuffer)
{
	size_t numVerts = softbody.visMesh->numVerts;
    size_t numTri = softbody.visMesh->numTri;

	mesh.vertex_positions.resize(numVerts);

	// Copy updated visual mesh positions into the mesh component
    for (size_t i = 0; i < numVerts; ++i)
    {
		// first update the vertex positions of the visual mesh with the current positions of
		// the wire mesh (which are updated by physics simulation)
		softbody.visMesh->verts[i] = softbody.wireMesh->pos[i];
        mesh.vertex_positions[i] = softbody.visMesh->verts[i];
    }

	// Recompute normals
    mesh.vertex_normals.resize(mesh.vertex_positions.size(), XMFLOAT3(0, 0, 0));
    for (size_t t = 0; t < numTri; ++t)
    {
		const XMINT3& tri = softbody.visMesh->triIds[t];
		uint32_t i0 = tri.x;
		uint32_t i1 = tri.y;
		uint32_t i2 = tri.z;

        const XMFLOAT3 &p0 = mesh.vertex_positions[i0];
        const XMFLOAT3 &p1 = mesh.vertex_positions[i1];
        const XMFLOAT3 &p2 = mesh.vertex_positions[i2];

        XMVECTOR v0 = XMLoadFloat3(&p0);
        XMVECTOR v1 = XMLoadFloat3(&p1);
        XMVECTOR v2 = XMLoadFloat3(&p2);

        XMVECTOR n = XMVector3Cross(v2 - v0, v1 - v0);
        n = XMVector3Normalize(n);

        XMFLOAT3 fn;
        XMStoreFloat3(&fn, n);

        mesh.vertex_normals[i0].x += fn.x;
        mesh.vertex_normals[i0].y += fn.y;
        mesh.vertex_normals[i0].z += fn.z;

        mesh.vertex_normals[i1].x += fn.x;
        mesh.vertex_normals[i1].y += fn.y;
        mesh.vertex_normals[i1].z += fn.z;

        mesh.vertex_normals[i2].x += fn.x;
        mesh.vertex_normals[i2].y += fn.y;
        mesh.vertex_normals[i2].z += fn.z;
    }

	// Normalize normals
    for (size_t i = 0; i < mesh.vertex_normals.size(); ++i)
    {
        XMVECTOR n = XMLoadFloat3(&mesh.vertex_normals[i]);
        n = XMVector3Normalize(n);
        XMStoreFloat3(&mesh.vertex_normals[i], n);
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

    // Update the GPU buffer where are stored the vertex positions and normals
    // of the mesh
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

        device->UpdateBuffer(
            &mesh.generalBuffer, pos16.data(), cmd,
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

        device->UpdateBuffer(
            &mesh.generalBuffer, pos32.data(), cmd,
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

    device->UpdateBuffer(&mesh.generalBuffer, nor32.data(), cmd,
                         mesh.vb_nor.size, mesh.vb_nor.offset);
}

void update_tetMesh(const WireMesh& wireMesh, wi::scene::MeshComponent &mesh,
                    bool updateGPUBuffer)
{
    size_t numVerts = wireMesh.pos.size();
    mesh.vertex_positions.resize(numVerts);

	// Copy tet mesh positions into the mesh component
    for (size_t i = 0; i < numVerts; ++i)
    {
        mesh.vertex_positions[i] = wireMesh.pos[i];
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

    // Update the GPU buffer where are stored the vertex positions and normals
    // of the mesh
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

        device->UpdateBuffer(
            &mesh.generalBuffer, pos16.data(), cmd,
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

        device->UpdateBuffer(
            &mesh.generalBuffer, pos32.data(), cmd,
            mesh.vb_pos_wind.size, // pos32.size() * sizeof(wi::scene::MeshComponent::Vertex_POS32),
            mesh.vb_pos_wind.offset);
    }
    else
    {
        assert("Unsupported vertex format!" && false);
    }
}

void create_wireMesh(std::unique_ptr<SimulationObject> &instance, const SimulationParams &params)
{
    using namespace wi::scene;
    using namespace wi::ecs;

	// Create a material for the triangular mesh (i.e., rendered in wireframe mode) with the custom shader created for this sample.
    if (wire_mat_entity == INVALID_ENTITY)
    {
        wire_mat_entity = CreateEntity();
        GetScene().materials.Create(wire_mat_entity);
        MaterialComponent &material = *GetScene().materials.GetComponent(wire_mat_entity);
        material.shaderType = MaterialComponent::SHADERTYPE_UNLIT;
        material.baseColor = XMFLOAT4(1.0f, 1.0f, 1.0f, 1.0f);
        material.SetCustomShaderID(wire_shader_id);
    }

    // Create a new entity for the triangular mesh
    // This entity (integer) will be used to create and access all components
    // related to the wire mesh.
    Entity wireEntity = CreateEntity();

    GetScene().layers.Create(wireEntity);     // Layer component is required if you want to filter entities (see doc.)
    GetScene().transforms.Create(wireEntity); // Transform component is required for applying translations, rotations, scalings, etc. to the object

    // A mesh component stores the mesh data (vertices, indices, subsets, etc.).
    // An object component represents a renderable object in the scene and holds
    // the mesh component entity that describes its geometry. Note that in this
    // sample we compute the physics simulation on the TetraMesh member of the
	// SoftBodySkinning object (which is not the mesh component), and use it to
	// update data in the tetrahedral mesh component.
	// So we cannot use multiple instance object compponents sharing the same
	// mesh component (we need a mesh per object).
    ObjectComponent &tetObject = GetScene().objects.Create(wireEntity);
    MeshComponent &wireMesh = GetScene().meshes.Create(wireEntity);
    tetObject.meshID = wireEntity;

    // instace that will hold the triangular object and the entity we will use to
    // access the mesh and other components associated with this triangular object.
    // WireMesh will store and manage all physics-related data and operations.
    instance->wireMesh = std::make_unique<WireMesh>(getMeshData(), params.stretchingCompliance, params.bendingCompliance);
	instance->wireEntity = wireEntity;

    // A subset specifies a range of indices to describe part or whole
    // of the mesh to be rendered with a specific material.
	// Note that we will render all the triangles of the triangular mesh in
	// wireframe mode so we set the indexCount to the number of triangle indices.
    wireMesh.subsets.push_back(wi::scene::MeshComponent::MeshSubset());
    wireMesh.subsets.back().materialID = wire_mat_entity;
    wireMesh.subsets.back().indexOffset = 0;
    wireMesh.subsets.back().indexCount = instance->wireMesh->faceTriIds.size();

	uint32_t numverts = static_cast<uint32_t>(instance->wireMesh->numVerts);
	uint32_t numIndices = static_cast<uint32_t>(instance->wireMesh->faceTriIds.size());

	// Load the wire mesh data into the mesh component
	wireMesh.vertex_positions.resize(numverts);
	for (size_t i = 0; i < numverts; ++i)
	{
		wireMesh.vertex_positions[i] = instance->wireMesh->pos[i];
	}

	wireMesh.indices.resize(numIndices);
	for (size_t i = 0; i < numIndices; ++i)
	{
		wireMesh.indices[i] = instance->wireMesh->faceTriIds[i];
	}

    // Apply initial transformations to the triangular object according to the given parameters
	// We don't use the transform component corresponding to wireEntity because the physics simulation
	// is computed on the vertex positions stored in the triangular object.
    instance->wireMesh->Scale(params.scale);
    instance->wireMesh->RotateRollPitchYaw(params.rotate);
    instance->wireMesh->Translate(params.translate);

    // Initialize the physics-related data structures of the triangular object
    instance->wireMesh->InitPhysics();

    // Initially hide the triangular mesh
    tetObject.SetRenderable(false);

    // Update mesh data according to the current state of the triangular object.
    update_tetMesh(*instance->wireMesh, wireMesh, false);

    // Upload the mesh data to a buffer created as a sub-resource of a general
    // GPU buffer containing all mesh data of the scene, if available.
    // Otherwise, create a new GPU buffer for this mesh only.
    wireMesh.CreateRenderData();
}

void create_visMesh(std::unique_ptr<SimulationObject> &instance, const SimulationParams &params)
{
    using namespace wi::scene;
    using namespace wi::ecs;

    // Create a material for the visual mesh with the custom shader created for this sample.
    if (visual_mat_entity == INVALID_ENTITY)
    {
        visual_mat_entity = CreateEntity();
        GetScene().materials.Create(visual_mat_entity);
        MaterialComponent &material = *GetScene().materials.GetComponent(visual_mat_entity);
        material.shaderType = MaterialComponent::SHADERTYPE_PBR;
        material.baseColor = XMFLOAT4(0.9686f, 0.1412f, 0.1137f, 1.0f);
        material.roughness = 0.8f;
        // material.SetDoubleSided(true);
        material.SetCustomShaderID(vis_shader_id);
    }

    // Create a new entity for the visual mesh
    // This entity (integer) will be used to create and access all components
    // related to the visual mesh.
    Entity visEntity = CreateEntity();
    GetScene().layers.Create(visEntity);
    GetScene().transforms.Create(visEntity);

    ObjectComponent &visObject = GetScene().objects.Create(visEntity);
    MeshComponent &visMesh = GetScene().meshes.Create(visEntity);
    visObject.meshID = visEntity;

    // instace that will hold the VisualMesh object and the entity we will use to
    // access the mesh and other components associated with this VisualMesh object.
    // VisualMesh will store and manage all visual mesh-related data.
    instance->visMesh = std::make_unique<VisMesh>(getMeshData());
    instance->visEntity = visEntity;

	// Note that here we set the indexCount to the number of triangle indices
    visMesh.subsets.push_back(wi::scene::MeshComponent::MeshSubset());
    visMesh.subsets.back().materialID = visual_mat_entity;
    visMesh.subsets.back().indexOffset = 0;
    visMesh.subsets.back().indexCount = instance->visMesh->numIndices;

    // Apply the same transformations to visual mesh vertices
	// Note that we don't use the transform component of the visEntity because
	// the VisualMesh data is stored and managed inside the SoftBodySkinning object,
	// which is used to update the visual mesh vertex positions after skinning.
    XMMATRIX scaleMat = XMMatrixScaling(params.scale.x, params.scale.y, params.scale.z);
    XMMATRIX rotMat = XMMatrixRotationRollPitchYaw(params.rotate.y, params.rotate.z, params.rotate.x);
    XMVECTOR transl = XMLoadFloat3(&params.translate);

    for (size_t i = 0; i < instance->visMesh->verts.size(); ++i)
    {
        XMVECTOR v = XMLoadFloat3(&instance->visMesh->verts[i]);
        v = XMVector3Transform(v, scaleMat);
        v = XMVector3Transform(v, rotMat);
        v = XMVectorAdd(v, transl);
        XMStoreFloat3(&instance->visMesh->verts[i], v);
    }

	// Load the visual mesh data into the mesh component
	size_t numVisVerts = instance->visMesh->numVerts;
	size_t numTri = instance->visMesh->numTri;
	size_t numIndices = instance->visMesh->numIndices;
    visMesh.vertex_positions.resize(numVisVerts);
    visMesh.vertex_normals.resize(numVisVerts);
    visMesh.indices.resize(numIndices);

    for (size_t i = 0; i < numVisVerts; ++i)
    {
		visMesh.vertex_positions[i] = instance->visMesh->verts[i];
    }

	for (size_t i = 0; i < numIndices; ++i)
	{
		visMesh.indices[i] = instance->visMesh->indices[i];
	}

	// Compute initial normals for the visual mesh
    for (size_t t = 0; t < numTri; ++t)
    {
		uint32_t i0 = visMesh.indices[t * 3 + 0];
		uint32_t i1 = visMesh.indices[t * 3 + 1];
		uint32_t i2 = visMesh.indices[t * 3 + 2];

        const XMFLOAT3 &p0 = visMesh.vertex_positions[i0];
        const XMFLOAT3 &p1 = visMesh.vertex_positions[i1];
        const XMFLOAT3 &p2 = visMesh.vertex_positions[i2];

        XMVECTOR v0 = XMLoadFloat3(&p0);
        XMVECTOR v1 = XMLoadFloat3(&p1);
        XMVECTOR v2 = XMLoadFloat3(&p2);

		// We are in a left-handed coordinate system, so the vector resulting from the cross product
		// "sees" the first operand vector rotated clockwise towards the second operand vector.
		// In this case, the vertices of each triangle are given in counter-clockwise order, so the first
		// and second vectors will be, respectively, (v2 - v0) and (v1 - v0).
        XMVECTOR n = XMVector3Cross(v2 - v0, v1 - v0);
        n = XMVector3Normalize(n);

        XMFLOAT3 fn;
        XMStoreFloat3(&fn, n);

        visMesh.vertex_normals[i0].x += fn.x;
        visMesh.vertex_normals[i0].y += fn.y;
        visMesh.vertex_normals[i0].z += fn.z;

        visMesh.vertex_normals[i1].x += fn.x;
        visMesh.vertex_normals[i1].y += fn.y;
        visMesh.vertex_normals[i1].z += fn.z;

        visMesh.vertex_normals[i2].x += fn.x;
        visMesh.vertex_normals[i2].y += fn.y;
        visMesh.vertex_normals[i2].z += fn.z;
    }

    // Normalize all normals
    for (size_t i = 0; i < visMesh.vertex_normals.size(); ++i)
    {
        XMVECTOR n = XMLoadFloat3(&visMesh.vertex_normals[i]);
        n = XMVector3Normalize(n);
        XMStoreFloat3(&visMesh.vertex_normals[i], n);
    }

    visMesh.CreateRenderData();
}

std::unique_ptr<SimulationObject> create_simulation_object(const SimulationParams &params, const std::string &modelPath)
{
    using namespace wi::scene;
    using namespace wi::ecs;

    // Create a Simulation object
    auto instance = std::make_unique<SimulationObject>();
	gMeshPath = modelPath;

    create_wireMesh(instance, params);
    create_visMesh(instance, params);

    return instance;
}

void remove_simulation_object(SimulationObject &sbi)
{
    using namespace wi::scene;
    using namespace wi::ecs;

    if (visual_mat_entity != INVALID_ENTITY)
    {
        GetScene().materials.Remove(visual_mat_entity);
        visual_mat_entity = INVALID_ENTITY;
    }

	if (wire_mat_entity != INVALID_ENTITY)
	{
		GetScene().materials.Remove(wire_mat_entity);
		wire_mat_entity = INVALID_ENTITY;
	}

    GetScene().layers.Remove(sbi.wireEntity);
    GetScene().transforms.Remove(sbi.wireEntity);
    GetScene().objects.Remove(sbi.wireEntity);
    GetScene().meshes.Remove(sbi.wireEntity);

	GetScene().layers.Remove(sbi.visEntity);
	GetScene().transforms.Remove(sbi.visEntity);
	GetScene().objects.Remove(sbi.visEntity);
	GetScene().meshes.Remove(sbi.visEntity);
}
} // namespace simulation
