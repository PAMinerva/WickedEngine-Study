#include "simulation_utils.h"
#include "cloth.h"
#include <algorithm>
#include <memory>
#include <string>
#include <cfloat>
#include <wiScene_Components.h>
#include <wiGraphics.h>
#include <wiMath.h>
#include <wiPrimitive.h>

static uint64_t wire_shader_id = 0;

simulation::PhysicsScene gPhysicsScene;
wi::gui::Label label_tris;
wi::gui::Label label_verts;

namespace simulation
{

void init_simulation(uint64_t wireShaderID)
{
    wire_shader_id = wireShaderID;

    SimulationParams params{};

    auto instance = create_simulation_object(params);

    int numTris = instance->cloth->numTris;
    int numVerts = instance->cloth->numParticles;

    instance->cloth->InitGPUBuffers();

    gPhysicsScene.objects.push_back(std::move(instance));

    label_tris.SetText(std::to_string(numTris) + " triangles");
    label_verts.SetText(std::to_string(numVerts) + " vertices");
}

void simulate(wi::graphics::CommandList cmd, float frameDt)
{
    if (gPhysicsScene.paused)
        return;

    for (auto& object : gPhysicsScene.objects)
    {
        object->cloth->gravity = gPhysicsScene.gravity;
        object->cloth->numSubSteps = gPhysicsScene.numSubsteps;
		object->cloth->dt = gPhysicsScene.dt;
        // object->cloth->dt = frameDt;

        object->cloth->params.sphereCenter[0] = gPhysicsScene.sphereCenter[0];
        object->cloth->params.sphereCenter[1] = gPhysicsScene.sphereCenter[1];
        object->cloth->params.sphereCenter[2] = gPhysicsScene.sphereCenter[2];
        object->cloth->params.sphereRadius = gPhysicsScene.sphereRadius;

		// Simulate cloth on GPU and update GPU buffers with new positions and normals
        object->cloth->SimulateGPU(frameDt, cmd, gPhysicsScene.solveType);
        object->cloth->UpdateMeshNormalsGPU(cmd);
    }
}

void create_wire_mesh(std::unique_ptr<SimulationObject>& instance, const SimulationParams& params)
{
    using namespace wi::scene;
    using namespace wi::ecs;

    Scene& scene = GetScene();
    ClothMesh& cloth = *instance->cloth;

    Entity entity = CreateEntity();
    instance->wireEntity = entity;

    TransformComponent& transform = scene.transforms.Create(entity);
    transform.UpdateTransform();

    MeshComponent& mesh = scene.meshes.Create(entity);
    mesh.subsets.resize(1);
    mesh.subsets[0].indexCount = static_cast<uint32_t>(cloth.triIds.size());
    mesh.subsets[0].indexOffset = 0;
    mesh.subsets[0].materialID = CreateEntity();

    MaterialComponent& material = scene.materials.Create(mesh.subsets[0].materialID);
    material.baseColor = XMFLOAT4(1.0f, 0.0f, 0.0f, 1.0f);
    material.shaderType = MaterialComponent::SHADERTYPE_PBR;
    material.customShaderID = wire_shader_id;

    mesh.vertex_positions.resize(cloth.numParticles);
    for (int i = 0; i < cloth.numParticles; i++)
    {
        mesh.vertex_positions[i] = XMFLOAT3(cloth.cpuPos[i].x, cloth.cpuPos[i].y, cloth.cpuPos[i].z);
    }

    mesh.vertex_normals.resize(cloth.numParticles);
    for (int i = 0; i < cloth.numParticles; i++)
    {
        mesh.vertex_normals[i] = cloth.cpuNormals[i];
    }

    mesh.indices.resize(cloth.triIds.size());
    for (size_t i = 0; i < cloth.triIds.size(); i++)
    {
        mesh.indices[i] = cloth.triIds[i];
    }

    // Disable quantization for vertex positions, because (at least in my testing),
    // during collision, 16-bit quantization creates too much visual artifacts.
    mesh.SetQuantizedPositionsDisabled(true);
    mesh.CreateRenderData();

    ObjectComponent& obj = scene.objects.Create(entity);
    obj.meshID = entity;
    obj.SetRenderable(true);
    obj.SetCastShadow(false);

    LayerComponent& layer = scene.layers.Create(entity);
    layer.layerMask = 1;
}

void update_mesh(const ClothMesh& cloth, wi::scene::MeshComponent& mesh, wi::graphics::CommandList cmd)
{
    wi::graphics::GraphicsDevice* device = wi::graphics::GetDevice();

    XMFLOAT3 _min = {  FLT_MAX,  FLT_MAX,  FLT_MAX };
    XMFLOAT3 _max = { -FLT_MAX, -FLT_MAX, -FLT_MAX };

    for (int i = 0; i < cloth.numParticles; i++)
    {
        XMFLOAT3 p(cloth.cpuPos[i].x, cloth.cpuPos[i].y, cloth.cpuPos[i].z);
        mesh.vertex_positions[i] = p;
        _min = wi::math::Min(_min, p);
        _max = wi::math::Max(_max, p);
    }

    mesh.aabb = wi::primitive::AABB(_min, _max);

    if (mesh.vertex_normals.size() == static_cast<size_t>(cloth.numParticles))
    {
        for (int i = 0; i < cloth.numParticles; i++)
        {
            mesh.vertex_normals[i] = cloth.cpuNormals[i];
        }
    }

    // Update positions in generalBuffer (POS32 — full float, no quantization)
    {
        wi::vector<wi::scene::MeshComponent::Vertex_POS32> verts(cloth.numParticles);
        for (int i = 0; i < cloth.numParticles; i++)
            verts[i].FromFULL(mesh.vertex_positions[i]);

        device->UpdateBuffer(
            &mesh.generalBuffer,
            verts.data(),
            cmd,
            mesh.vb_pos_wind.size,
            mesh.vb_pos_wind.offset
        );
    }

    // Update normals in generalBuffer
    {
        wi::vector<wi::scene::MeshComponent::Vertex_NOR> nors(cloth.numParticles);
        for (int i = 0; i < cloth.numParticles; i++)
            nors[i].FromFULL(mesh.vertex_normals[i]);

        device->UpdateBuffer(
            &mesh.generalBuffer,
            nors.data(),
            cmd,
            mesh.vb_nor.size,
            mesh.vb_nor.offset
        );
    }
}

std::unique_ptr<SimulationObject> create_simulation_object(const SimulationParams& params)
{
    auto instance = std::make_unique<SimulationObject>();
    instance->cloth = std::make_unique<ClothMesh>(params);
    create_wire_mesh(instance, params);
    return instance;
}

void remove_simulation_object(SimulationObject& obj)
{
    using namespace wi::scene;
    using namespace wi::ecs;

    Scene& scene = GetScene();

    if (obj.wireEntity != INVALID_ENTITY)
    {
        auto mesh = scene.meshes.GetComponent(obj.wireEntity);
        if (mesh)
        {
            if (mesh->subsets.size() > 0 && mesh->subsets[0].materialID != INVALID_ENTITY)
            {
                scene.materials.Remove(mesh->subsets[0].materialID);
                scene.Entity_Remove(mesh->subsets[0].materialID);
            }
        }

        scene.meshes.Remove(obj.wireEntity);
        scene.objects.Remove(obj.wireEntity);
        scene.transforms.Remove(obj.wireEntity);
        scene.layers.Remove(obj.wireEntity);
        scene.Entity_Remove(obj.wireEntity);
    }

    obj.cloth.reset();
}

} // namespace simulation
