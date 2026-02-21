#include "simulation_utils.h"
#include "cloth.h"
#include <memory>
#include <random>
#include <string>
#include <wiScene_Components.h>

static wi::ecs::Entity front_mat_entity = wi::ecs::INVALID_ENTITY;
static wi::ecs::Entity back_mat_entity = wi::ecs::INVALID_ENTITY;
static uint64_t front_shader_id = 0;
static uint64_t back_shader_id = 0;
static uint64_t wire_shader_id = 0;

simulation::PhysicsScene gPhysicsScene;
// wi::gui::Label label_tets;
wi::gui::Label label_tris;
wi::gui::Label label_verts;

namespace simulation
{
void init_simulation(uint64_t frontShaderID, uint64_t backShaderID, uint64_t wireShaderID)
{
    front_shader_id = frontShaderID;
    back_shader_id = backShaderID;
	wire_shader_id = wireShaderID;

    SimulationParams params{};
    // Set the parameters of the simulation according to your preferences.
    // params.numX = 50;
    // params.numY = 50;
    // params.spacing = 0.02f;
    // params.thickness = 0.01f;
    // params.bendingCompliance = 1.0f;

    auto instance = create_simulation_object(params);

    int numTris = instance->cloth->numTris;
    int numVerts = instance->cloth->numParticles;

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

    for (auto& object : gPhysicsScene.objects)
    {
        object->cloth->Simulate(gPhysicsScene.dt, gPhysicsScene.numSubsteps, gPhysicsScene.gravity);
    }
}

void create_front_mesh(std::unique_ptr<SimulationObject>& instance, const SimulationParams& params)
{
    using namespace wi::scene;
    using namespace wi::ecs;

    Scene& scene = GetScene();
    ClothMesh& cloth = *instance->cloth;

    // Create entity
    Entity entity = CreateEntity();
    instance->frontEntity = entity;

    // Create transform
    TransformComponent& transform = scene.transforms.Create(entity);
    transform.Translate(params.translate);
    transform.UpdateTransform();

    // Create mesh
    MeshComponent& mesh = scene.meshes.Create(entity);
    mesh.subsets.resize(1);
    mesh.subsets[0].indexCount = static_cast<uint32_t>(cloth.triIds.size());
    mesh.subsets[0].indexOffset = 0;
    mesh.subsets[0].materialID = CreateEntity();

    // Create material (red for front face)
    MaterialComponent& material = scene.materials.Create(mesh.subsets[0].materialID);
    material.baseColor = XMFLOAT4(1.0f, 0.0f, 0.0f, 1.0f);  // Red
    material.roughness = 0.8f;
    material.SetCastShadow(true);

    // Set custom shader
    material.customShaderID = front_shader_id;

    // Setup vertex positions
    mesh.vertex_positions.resize(cloth.numParticles);
    for (int i = 0; i < cloth.numParticles; i++)
    {
        mesh.vertex_positions[i] = cloth.pos[i];
    }

    // Setup indices
    mesh.indices.resize(cloth.triIds.size());
    for (size_t i = 0; i < cloth.triIds.size(); i++)
    {
        mesh.indices[i] = cloth.triIds[i];
    }

    // Compute normals
    mesh.ComputeNormals(MeshComponent::COMPUTE_NORMALS_SMOOTH_FAST);

    // Create GPU buffers
    mesh.CreateRenderData();

    // Create object component (makes it renderable)
    ObjectComponent& obj = scene.objects.Create(entity);
    obj.meshID = entity;
    obj.SetRenderable(true);
    obj.SetCastShadow(true);

    // Add to layer
    LayerComponent& layer = scene.layers.Create(entity);
    layer.layerMask = 1;
}

void create_back_mesh(std::unique_ptr<SimulationObject>& instance, const SimulationParams& params)
{
    using namespace wi::scene;
    using namespace wi::ecs;

    Scene& scene = GetScene();
    ClothMesh& cloth = *instance->cloth;

    // Create entity
    Entity entity = CreateEntity();
    instance->backEntity = entity;

    // Create transform
    TransformComponent& transform = scene.transforms.Create(entity);
    transform.Translate(params.translate);
    transform.UpdateTransform();

    // Create mesh
    MeshComponent& mesh = scene.meshes.Create(entity);
    mesh.subsets.resize(1);
    mesh.subsets[0].indexCount = static_cast<uint32_t>(cloth.triIds.size());
    mesh.subsets[0].indexOffset = 0;
    mesh.subsets[0].materialID = CreateEntity();

    // Create material (yellow-ish for back face)
    MaterialComponent& material = scene.materials.Create(mesh.subsets[0].materialID);
    material.baseColor = XMFLOAT4(1.0f, 0.5f, 0.0f, 1.0f);  // Yellow-ish
    material.roughness = 0.8f;
    material.SetCastShadow(true);

    // Set custom shader
    material.customShaderID = back_shader_id;

    // Setup vertex positions (same as front)
    mesh.vertex_positions.resize(cloth.numParticles);
    for (int i = 0; i < cloth.numParticles; i++)
    {
        mesh.vertex_positions[i] = cloth.pos[i];
    }

    // Setup indices (same as front)
    mesh.indices.resize(cloth.triIds.size());
    for (size_t i = 0; i < cloth.triIds.size(); i++)
    {
        mesh.indices[i] = cloth.triIds[i];
    }

    // Compute normals
    mesh.ComputeNormals(MeshComponent::COMPUTE_NORMALS_SMOOTH_FAST);

    // Create GPU buffers
    mesh.CreateRenderData();

    // Create object component
    ObjectComponent& obj = scene.objects.Create(entity);
    obj.meshID = entity;
    obj.SetRenderable(true);
    obj.SetCastShadow(true);  // Back face doesn't cast shadow

    // Add to layer
    LayerComponent& layer = scene.layers.Create(entity);
    layer.layerMask = 1;
}

void create_wire_mesh(std::unique_ptr<SimulationObject>& instance, const SimulationParams& params)
{
    using namespace wi::scene;
    using namespace wi::ecs;

    Scene& scene = GetScene();
    ClothMesh& cloth = *instance->cloth;

    // Create entity
    Entity entity = CreateEntity();
    instance->wireEntity = entity;

    // Create transform
    TransformComponent& transform = scene.transforms.Create(entity);
    transform.Translate(params.translate);
    transform.UpdateTransform();

    // Create mesh
    MeshComponent& mesh = scene.meshes.Create(entity);
    mesh.subsets.resize(1);
    mesh.subsets[0].indexCount = static_cast<uint32_t>(cloth.triIds.size());
    mesh.subsets[0].indexOffset = 0;
    mesh.subsets[0].materialID = CreateEntity();

    // Create material (white wireframe)
    MaterialComponent& material = scene.materials.Create(mesh.subsets[0].materialID);
    material.baseColor = XMFLOAT4(1.0f, 1.0f, 1.0f, 1.0f);  // White
    material.shaderType = MaterialComponent::SHADERTYPE_UNLIT;

    // Set custom shader
    material.customShaderID = wire_shader_id;

    // Setup vertex positions
    mesh.vertex_positions.resize(cloth.numParticles);
    for (int i = 0; i < cloth.numParticles; i++)
    {
        mesh.vertex_positions[i] = cloth.pos[i];
    }

    // Setup indices
    mesh.indices.resize(cloth.triIds.size());
    for (size_t i = 0; i < cloth.triIds.size(); i++)
    {
        mesh.indices[i] = cloth.triIds[i];
    }

    // Create GPU buffers
    mesh.CreateRenderData();

    // Create object component
    ObjectComponent& obj = scene.objects.Create(entity);
    obj.meshID = entity;
    obj.SetRenderable(false);  // Hidden by default
    obj.SetCastShadow(false);

    // Add to layer
    LayerComponent& layer = scene.layers.Create(entity);
    layer.layerMask = 1;
}

void update_mesh(const ClothMesh& cloth, wi::scene::MeshComponent& mesh, bool updateGPUBuffer)
{
    // Update vertex positions from cloth simulation
    for (int i = 0; i < cloth.numParticles; i++)
    {
        mesh.vertex_positions[i] = cloth.pos[i];
    }

    // Recompute normals
    mesh.ComputeNormals(wi::scene::MeshComponent::COMPUTE_NORMALS_SMOOTH_FAST);

    if (updateGPUBuffer)
    {
        mesh.CreateRenderData();
    }
}

std::unique_ptr<SimulationObject> create_simulation_object(const SimulationParams &params)
{
    using namespace wi::scene;
    using namespace wi::ecs;

    // Create a SimulationObject
    auto instance = std::make_unique<SimulationObject>();

    // Create the Cloth (generates grid, constraints, and triangle indices)
    instance->cloth = std::make_unique<ClothMesh>(params);

    // Create the visual meshes
    create_front_mesh(instance, params);
    create_back_mesh(instance, params);
	create_wire_mesh(instance, params);

    return instance;
}

void remove_simulation_object(SimulationObject& obj)
{
    using namespace wi::scene;
    using namespace wi::ecs;

    Scene& scene = GetScene();

    // Remove front entity and its components
    if (obj.frontEntity != INVALID_ENTITY)
    {
        auto mesh = scene.meshes.GetComponent(obj.frontEntity);
        if (mesh)
        {
            // Remove material
            if (mesh->subsets.size() > 0 && mesh->subsets[0].materialID != INVALID_ENTITY)
            {
                scene.materials.Remove(mesh->subsets[0].materialID);
                scene.Entity_Remove(mesh->subsets[0].materialID);
            }
        }

        scene.meshes.Remove(obj.frontEntity);
        scene.objects.Remove(obj.frontEntity);
        scene.transforms.Remove(obj.frontEntity);
        scene.layers.Remove(obj.frontEntity);
        scene.Entity_Remove(obj.frontEntity);
    }

    // Remove back entity and its components
    if (obj.backEntity != INVALID_ENTITY)
    {
        auto mesh = scene.meshes.GetComponent(obj.backEntity);
        if (mesh)
        {
            if (mesh->subsets.size() > 0 && mesh->subsets[0].materialID != INVALID_ENTITY)
            {
                scene.materials.Remove(mesh->subsets[0].materialID);
                scene.Entity_Remove(mesh->subsets[0].materialID);
            }
        }

        scene.meshes.Remove(obj.backEntity);
        scene.objects.Remove(obj.backEntity);
        scene.transforms.Remove(obj.backEntity);
        scene.layers.Remove(obj.backEntity);
        scene.Entity_Remove(obj.backEntity);
    }

    // Remove wire entity and its components
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

    // Reset the cloth pointer
    obj.cloth.reset();
}
} // namespace simulation
