#include "simulation_utils.h"
#include "hashing.h"
#include "meshes.h"
#include "softBody.h"
#include <memory>
#include <random>
#include <string>
#include <wiScene_Components.h>

static wi::ecs::Entity visual_mat_entity = wi::ecs::INVALID_ENTITY;
static wi::ecs::Entity wire_mat_entity = wi::ecs::INVALID_ENTITY;
static uint64_t vis_shader_id = 0;
static uint64_t wire_shader_id = 0;

simulation::PhysicsScene gPhysicsScene;
wi::gui::Label label_tets;
wi::gui::Label label_tris;
wi::gui::Label label_verts;

namespace simulation
{
void init_physics(uint64_t visShaderID, uint64_t wireShaderID)
{
    vis_shader_id = visShaderID;
    wire_shader_id = wireShaderID;

    SoftBodySkinningParams params{};
    params.rotate = {0.0f, 0.0f, 3.14f}; // rotate 180 degrees around Z axis clockwise (180 x (pi/180) = pi = 3.14 rad)
    params.edgeCompliance = 0.0f;
    params.volCompliance = 0.0f;

    auto instance = create_softbody_skinning(params);
    int numTets = instance->tetMesh->numTets;
    int numTris = instance->visMesh->numTris;
    int numVerts = instance->visMesh->numVerts;

    gPhysicsScene.objects.push_back(std::move(instance));

    std::string text = std::to_string(numTets) + " tetrahedrons";
    label_tets.SetText(text);

    text = std::to_string(numTris) + " triangles";
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
            instance->tetMesh->PreSolve(subdt, gPhysicsScene.gravity);

        // 2. Solve: project constraints (edge, volume, etc.)
        // Projecting a set of points according to a constraint means moving
        // (correcting) the points such that they satisfy the constraint
        for (auto &instance : gPhysicsScene.objects)
            instance->tetMesh->Solve(subdt);

        // 3. PostSolve: update velocity
        for (auto &instance : gPhysicsScene.objects)
            instance->tetMesh->PostSolve(subdt);
    }
}

void new_softbody_skinning()
{
    float a = 3.0f;

    // Random number generator setup
    static std::random_device rd;
    static std::mt19937 gen(rd());

    std::uniform_real_distribution<float> dist_x(-a, a);
    std::uniform_real_distribution<float> dist_y(0, a);
    std::uniform_real_distribution<float> dist_z(0, a);

    SoftBodySkinningParams params{};
    params.translate = {dist_x(gen), dist_y(gen), dist_z(gen)};
    params.edgeCompliance = 0.0f;
    params.volCompliance = 0.0f;

    auto instance = create_softbody_skinning(params);
    gPhysicsScene.objects.push_back(std::move(instance));

    size_t numTets = 0;
    for (auto &instance : gPhysicsScene.objects)
    {
        numTets += instance->tetMesh->numTets;
    }

	size_t numTris = 0;
	for (auto &instance : gPhysicsScene.objects)
	{
		numTris += instance->visMesh->numTris;
	}

	size_t numVerts = 0;
	for (auto &instance : gPhysicsScene.objects)
	{
		numVerts += instance->visMesh->numVerts;
	}

    std::string text = std::to_string(numTets) + " tetrahedrons";
    label_tets.SetText(text);

	text = std::to_string(numTris) + " triangles";
	label_tris.SetText(text);

	text = std::to_string(numVerts) + " vertices";
	label_verts.SetText(text);
}

void update_visMesh(const SoftBodySkinning &softbody,
                    wi::scene::MeshComponent &mesh, bool updateGPUBuffer)
{
    size_t numVerts = softbody.visMesh->numVerts;
    size_t numTris = getDragonVisMesh().triIds.size() / 3;
    std::vector<XMFLOAT3> &positions = softbody.visMesh->verts;
    const auto &skinningInfo = softbody.tetMesh->skinningInfo;
    const auto &tetIds = softbody.tetMesh->tetIds;
    const auto &pos = softbody.tetMesh->pos;

	// For each visual vertex...
    for (size_t i = 0; i < numVerts; i++)
    {
		// Get the nearest tetrahedron and barycentric coordinates, if any
        const XMFLOAT4 &skin = skinningInfo[i];
        int tetNr = static_cast<int>(skin.x);
        if (tetNr < 0)
        {
            continue;
        }
        float b0 = skin.y;
        float b1 = skin.z;
        float b2 = skin.w;
        float b3 = 1.0f - b0 - b1 - b2;

		// Get the tetrahedron vertex indices...
        const XMINT4 &tet = tetIds[tetNr];
        int id0 = tet.x;
        int id1 = tet.y;
        int id2 = tet.z;
        int id3 = tet.w;

		// ...and their current positions
        XMVECTOR v0 = XMLoadFloat3(&pos[id0]);
        XMVECTOR v1 = XMLoadFloat3(&pos[id1]);
        XMVECTOR v2 = XMLoadFloat3(&pos[id2]);
        XMVECTOR v3 = XMLoadFloat3(&pos[id3]);

		// Compute the new visual vertex position using barycentric interpolation
        XMVECTOR result = XMVectorZero();
        result = XMVectorAdd(result, XMVectorScale(v0, b0));
        result = XMVectorAdd(result, XMVectorScale(v1, b1));
        result = XMVectorAdd(result, XMVectorScale(v2, b2));
        result = XMVectorAdd(result, XMVectorScale(v3, b3));

        XMStoreFloat3(&positions[i], result);
    }

	// Copy updated visual mesh positions into the mesh component
    for (size_t i = 0; i < numVerts; ++i)
    {
        mesh.vertex_positions[i] = positions[i];
    }

	// Recompute normals
    mesh.vertex_normals.resize(mesh.vertex_positions.size(), XMFLOAT3(0, 0, 0));
    for (size_t t = 0; t < numTris; ++t)
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

        XMVECTOR n = XMVector3Cross(v1 - v0, v2 - v0);
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

void update_tetMesh(const TetraMesh& tetMesh, wi::scene::MeshComponent &mesh,
                    bool updateGPUBuffer)
{
    size_t numVerts = tetMesh.pos.size();
    mesh.vertex_positions.resize(numVerts);

	// Copy tet mesh positions into the mesh component
    for (size_t i = 0; i < numVerts; ++i)
    {
        mesh.vertex_positions[i] = tetMesh.pos[i];
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

void create_TetraDragon_Meshes()
{
    using namespace wi::scene;
    using namespace wi::ecs;
}

void computeSkinningInfo(SoftBodySkinning &softBodyInstance)
{
    // Initialize the spatial hashing structures for the visual mesh vertices
    auto hash = new SpatialHashing::Hash(0.05, softBodyInstance.visMesh->numVerts);
    hash->Create(softBodyInstance.visMesh->verts);

    std::vector<XMFLOAT4> &skinningInfo = softBodyInstance.tetMesh->skinningInfo;
    const std::vector<XMFLOAT3> &visVerts = softBodyInstance.visMesh->verts;
    const std::vector<XMFLOAT3> &pos = softBodyInstance.tetMesh->pos;
    const std::vector<XMINT4> &tetIds = softBodyInstance.tetMesh->tetIds;
    int numVisVerts = softBodyInstance.visMesh->numVerts;
    int numTets = softBodyInstance.tetMesh->numTets;

    std::fill(skinningInfo.begin(), skinningInfo.end(), XMFLOAT4(-1, -1, -1, -1));

	// The i-th element will contain the negative of the smallest barycentric coordinate
	// for the i-th visual mesh vertex across all tetrahedra.
	// When this value is zero or negative, it means that the vertex is inside at least one tetrahedron.
    std::vector<float> minDist(numVisVerts, std::numeric_limits<float>::max());
    float border = 0.05f;

	// For each tetrahedron...
    for (int i = 0; i < numTets; ++i)
    {
		// Compute tetrahedron center
        XMVECTOR tetCenter = XMVectorZero();
        tetCenter += XMLoadFloat3(&pos[tetIds[i].x]);
        tetCenter += XMLoadFloat3(&pos[tetIds[i].y]);
        tetCenter += XMLoadFloat3(&pos[tetIds[i].z]);
        tetCenter += XMLoadFloat3(&pos[tetIds[i].w]);
        tetCenter = tetCenter / 4.0f;

		// Compute maximum radius from center to a tet vertex, plus a border
        float rMax = 0.0f;
        XMVECTOR tetVerts[4] = {XMLoadFloat3(&pos[tetIds[i].x]),
								XMLoadFloat3(&pos[tetIds[i].y]),
								XMLoadFloat3(&pos[tetIds[i].z]),
								XMLoadFloat3(&pos[tetIds[i].w])};
        for (int j = 0; j < 4; ++j)
        {
            float dist = XMVectorGetX(XMVector3Length(tetCenter - tetVerts[j]));
            if (dist > rMax)
                rMax = dist;
        }
        rMax += border;

		// Query the spatial hash for visual mesh vertices within rMax of tetCenter
        XMFLOAT3 tetCenterF3;
        XMStoreFloat3(&tetCenterF3, tetCenter);
        std::vector<XMFLOAT3> queryCenters = {tetCenterF3};
        hash->Query(queryCenters, 0, rMax);

		// Get the list of visual mesh vertex IDs returned by the spatial hash
        const auto &queryIds = hash->GetQueryIds();
        if (queryIds.empty())
            continue;

		// Precompute the inverse matrix for barycentric coordinate calculation:
		//
		// A vertex v can be expressed in barycentric coordinates
		// with respect to the tetrahedron defined by points p0, p1, p2, p3 as:
		//
		//   v = b0*p0 + b1*p1 + b2*p2 + b3*p3      (1)
		//
		// where b0, b1, b2, b3 are the barycentric.
		// Also, if v is inside the tetrahedron, then each barycentric coordinate is >= 0.
		// We know v, p0, p1, p2, p3, and we want to compute b0, b1, b2, b3.
		// We can subtract p3 from both sides of (1). This is equivalent to translating
		// the tetrahedron and the vertex v by the same amount, which does not change
		// the barycentric coordinates:
		//
		//   v - p3 = b0*(p0 - p3) + b1*(p1 - p3) + b2*(p2 - p3)      (2)
		//
		// b3 disappears from the equation so that we have only three unknowns (b0, b1, b2)
		// that can be expressed as a vector b = (b0, b1, b2).
		// We can also create the matrix P whose rows are the vectors (p0 - p3), (p1 - p3), (p2 - p3):
		//
		//   P = [ (p0 - p3)  (p1 - p3)  (p2 - p3) ]^T
		//
		// Then, we can express the previous equation in matrix form:
		//
		//  v - p3 = P * b
		//
		// Solving for b:
		//
		//   b = P^-1 * (v - p3)
		//
		// Now we can derive b3 from (2) as:
		//
		//   v = b0*p0 + b1*p1 + b2*p2 - b0p3 - b1p3 - b2p3 + p3
		//     = b0*p0 + b1*p1 + b2*p2 + (1 - b0 - b1 - b2)*p3
		//
		// That is:
		//
		//  b3 = 1 - b0 - b1 - b2
		//
		// Note that the sum of all barycentric coordinates is always 1:
		//
		//   b0 + b1 + b2 + b3 = 1
		//
		// Anyway, we only need to explicitly find b0, b1, b2 computing the
		// inverse matrix P^-1 once per tetrahedron, and then multiplying it by (v - p3).
        XMVECTOR p0 = tetVerts[0];
        XMVECTOR p1 = tetVerts[1];
        XMVECTOR p2 = tetVerts[2];
        XMVECTOR p3 = tetVerts[3];

        XMVECTOR m0 = p0 - p3;
        XMVECTOR m1 = p1 - p3;
        XMVECTOR m2 = p2 - p3;

        XMMATRIX mat(XMVectorGetX(m0), XMVectorGetY(m0), XMVectorGetZ(m0), 0.0f,
                     XMVectorGetX(m1), XMVectorGetY(m1), XMVectorGetZ(m1), 0.0f,
                     XMVectorGetX(m2), XMVectorGetY(m2), XMVectorGetZ(m2), 0.0f,
                     0.0f, 0.0f, 0.0f, 1.0f);

        XMMATRIX invMat = XMMatrixInverse(nullptr, mat);

        for (size_t idx = 0; idx < queryIds.size(); ++idx)
        {
            int id = queryIds[idx];
            if (minDist[id] <= 0.0f) // We already found a tetrahedron containing this vertex
                continue;            // so continue to the next vertex

			// Get the queried visual mesh vertex position
			// and check if it is actually within rMax from tetCenter
            XMVECTOR v = XMLoadFloat3(&visVerts[id]);
            float dist2 = XMVectorGetX(XMVector3LengthSq(tetCenter - v));
            if (dist2 > rMax * rMax)
                continue;

			// Compute barycentric coordinates of v with respect to the current tetrahedron
            XMVECTOR diff = v - p3;
            XMVECTOR bary = XMVector3Transform(diff, invMat);

            float b0 = XMVectorGetX(bary);
            float b1 = XMVectorGetY(bary);
            float b2 = XMVectorGetZ(bary);
            float b3 = 1.0f - b0 - b1 - b2;

			// Find the minimum barycentric coordinate
			// and negate it
            float minBary = std::min({b0, b1, b2, b3});
            float dist = -minBary;

			// If this is the smallest negative barycentric coordinate...
            if (dist < minDist[id])
            {
				// Store it in minDist
				// and store the skinning info for this visual mesh vertex.
				// Note that we use the same integer ID to refer to the visual
				// mesh vertex, its skinning info, and its minDist value.
				// Also note that we store the tet index as a float in skinningInfo.x
				// and the three barycentric coordinates in skinningInfo.yzw
				// (the fourth barycentric coordinate can be computed from the other three).
                minDist[id] = dist;
                skinningInfo[id] = XMFLOAT4((float)i, b0, b1, b2);
            }
        }
    }
}

std::unique_ptr<SoftBodySkinning> create_softbody_skinning(const SoftBodySkinningParams &params)
{
    using namespace wi::scene;
    using namespace wi::ecs;

    // Create a material (shared by all visual mesh instances) that returns an orangish color when enlighted
    if (visual_mat_entity == INVALID_ENTITY)
    {
        visual_mat_entity = CreateEntity();
        GetScene().materials.Create(visual_mat_entity);
        MaterialComponent &material = *GetScene().materials.GetComponent(visual_mat_entity);
        material.shaderType = MaterialComponent::SHADERTYPE_PBR;
        material.baseColor = XMFLOAT4(0.9686f, 0.3412f, 0.1137f, 1.0f); // 0xF7571D; F7=247, 57=87, 1D=29; 247/255=0.9686, 87/255=0.3412, 29/255=0.1137
        // material.SetDoubleSided(true);
        material.SetCustomShaderID(vis_shader_id);
    }

	// Material that returns white color for tetrahedral meshes (wireframe mode)
    if (wire_mat_entity == INVALID_ENTITY)
    {
        wire_mat_entity = CreateEntity();
        GetScene().materials.Create(wire_mat_entity);
        MaterialComponent &material = *GetScene().materials.GetComponent(wire_mat_entity);
        material.shaderType = MaterialComponent::SHADERTYPE_UNLIT;
        material.baseColor = XMFLOAT4(1.0f, 1.0f, 1.0f, 1.0f);
        material.SetCustomShaderID(wire_shader_id);
    }

    // Create a new entity for the tetrahedral mesh (low-resolution).
    // This entity (integer) will be used to create and access all components
    // related to the tetrahedral mesh.
    Entity tetEntity = CreateEntity();

    GetScene().layers.Create(tetEntity);     // Layer component is required if you want to filter entities (see doc.)
    GetScene().transforms.Create(tetEntity); // Transform component is required for applying translations, rotations, scalings, etc. to the object

    // A mesh component stores the mesh data (vertices, indices, subsets, etc.).
    // An object component represents a renderable object in the scene and holds
    // the mesh component entity that describes its geometry. Note that in this
    // sample we compute the physics simulation on the TetraMesh member of the
	// SoftBodySkinning object (which is not the mesh component), and use it to
	// update data in the tetrahedral mesh component.
	// So we cannot use multiple instance object compponents sharing the same
	// mesh component (we need a mesh per object).
    ObjectComponent &tetObject = GetScene().objects.Create(tetEntity);
    MeshComponent &tetMesh = GetScene().meshes.Create(tetEntity);
    tetObject.meshID = tetEntity;

    // A subset specifies a range of indices to describe part or whole
    // of the mesh to be rendered with a specific material.
	// Note that we will render all the edges of the tetrahedral mesh in
	// wireframe mode so we set the indexCount to the number of edge indices.
    tetMesh.subsets.push_back(wi::scene::MeshComponent::MeshSubset());
    tetMesh.subsets.back().materialID = wire_mat_entity;
    tetMesh.subsets.back().indexOffset = 0;
    tetMesh.subsets.back().indexCount = (uint32_t)getDragonTetMesh().tetEdgeIds.size();

    // Create a SoftBodySkinning instace that will contain the TetraMesh
    // object and the entity we will use to access the mesh and other components
    // associated with this TetraMesh object. TetraMesh will store and manage all
    // physics-related data and operations.
    auto instance = std::make_unique<SoftBodySkinning>();
    instance->tetMesh = std::make_unique<TetraMesh>(getDragonTetMesh(), params.edgeCompliance, params.volCompliance);
    instance->tetEntity = tetEntity;

	// Load the tetrahedral mesh data into the mesh component
	size_t numVerts = instance->tetMesh->numVerts;
	size_t numSurfTriInds = getDragonTetMesh().tetEdgeIds.size();
    tetMesh.vertex_positions.resize(numVerts);
    tetMesh.indices.resize(numSurfTriInds);

    for (size_t i = 0; i < numVerts; ++i)
    {
		tetMesh.vertex_positions[i] = instance->tetMesh->pos[i];
    }

    for (size_t i = 0; i < numSurfTriInds; ++i)
    {
        tetMesh.indices[i] = (uint32_t)getDragonTetMesh().tetEdgeIds[i];
    }

    // Apply initial transformations to the TetraMesh object according to the given parameters
	// We don't use the transform component of the tetEntity because the physics simulation
	// is computed on the vertex positions stored in the TetraMesh object.
    instance->tetMesh->Scale(params.scale);
    instance->tetMesh->RotateRollPitchYaw(params.rotate);
    instance->tetMesh->Translate(params.translate);

    // Initialize the physics-related data structures of the TetraMesh object
    instance->tetMesh->InitPhysics();

    // Initially hide the tetrahedal mesh
    tetObject.SetRenderable(false);

    // Update mesh data according to the current state of the TetraMesh object.
    update_tetMesh(*instance->tetMesh, tetMesh, false);

    // Upload the mesh data to a buffer created as a sub-resource of a general
    // GPU buffer containing all mesh data of the scene, if available.
    // Otherwise, create a new GPU buffer for this mesh only.
    tetMesh.CreateRenderData();

    // Create a new entity for the VISUAL mesh (high-resolution)
    // This entity (integer) will be used to create and access all components
    // related to the visual mesh.
    Entity visEntity = CreateEntity();
    GetScene().layers.Create(visEntity);
    GetScene().transforms.Create(visEntity);

    ObjectComponent &visObject = GetScene().objects.Create(visEntity);
    MeshComponent &visMesh = GetScene().meshes.Create(visEntity);
    visObject.meshID = visEntity;

	// Note that here we set the indexCount to the number of triangle indices
    visMesh.subsets.push_back(wi::scene::MeshComponent::MeshSubset());
    visMesh.subsets.back().materialID = visual_mat_entity;
    visMesh.subsets.back().indexOffset = 0;
    visMesh.subsets.back().indexCount = (uint32_t)getDragonVisMesh().triIds.size();

    // Create the VISUAL mesh (high-resolution)
    instance->visMesh = std::make_unique<VisualMesh>(getDragonVisMesh());
    instance->visEntity = visEntity;

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

	// Compute the barycentric coordinates of each visual mesh vertex
	// with respect to the nearest tetrahedron of the tetrahedral mesh
	// This will be used for skinning the visual mesh during the simulation
	// (i.e., updating visual mesh vertex positions based on the current
	// positions of the tetrahedral mesh vertices).
    computeSkinningInfo(*instance);

	// Load the visual mesh data into the mesh component
	size_t numVisVerts = instance->visMesh->numVerts;
	size_t numTris = instance->visMesh->numTris;
    visMesh.vertex_positions.resize(numVisVerts);
    visMesh.vertex_normals.resize(numVisVerts);
    visMesh.indices.resize(numTris * 3);

    for (size_t i = 0; i < numVisVerts; ++i)
    {
		visMesh.vertex_positions[i] = instance->visMesh->verts[i];
    }

	for (size_t t = 0; t < numTris; ++t)
	{
		const XMINT3& tri = instance->visMesh->triIds[t];
		visMesh.indices[t * 3 + 0] = tri.x;
		visMesh.indices[t * 3 + 1] = tri.y;
		visMesh.indices[t * 3 + 2] = tri.z;
	}

	// Compute initial normals for the visual mesh
    for (size_t t = 0; t < numTris; ++t)
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

        XMVECTOR n = XMVector3Cross(v1 - v0, v2 - v0);
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

    return instance;
}

void remove_softbody_skinning(SoftBodySkinning &sbi)
{
    using namespace wi::scene;
    using namespace wi::ecs;

    if (visual_mat_entity != INVALID_ENTITY)
    {
        GetScene().materials.Remove(visual_mat_entity);
        visual_mat_entity = INVALID_ENTITY;
    }

    GetScene().layers.Remove(sbi.tetEntity);
    GetScene().transforms.Remove(sbi.tetEntity);
    GetScene().objects.Remove(sbi.tetEntity);
    GetScene().meshes.Remove(sbi.tetEntity);
}
} // namespace simulation
