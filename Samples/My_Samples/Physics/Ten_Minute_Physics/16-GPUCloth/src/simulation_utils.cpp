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

	// Setup renderer routing: so_pos/so_nor point to our render buffers
	{
		auto wireMesh = wi::scene::GetScene().meshes.GetComponent(instance->wireEntity);
		if (wireMesh)
			instance->cloth->SetupRendererRouting(*wireMesh);
	}

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

		// Simulate cloth on GPU and update GPU buffers with new positions and normals.
        object->cloth->SimulateGPU(frameDt, cmd, gPhysicsScene.solveType);
        object->cloth->UpdateMeshNormalsGPU(cmd);

		// Update renderer output buffers described by so_pos and so_nor views in the mesh component
		// with the new positions and normals from the simulation.
		// Note that the renderer automatically checks if the mesh component has so_pos/so_nor set up,
		// and if so, it uses the underlying buffers of those views (renderPosBuffer and renderNorBuffer in our case)
		// as vertex buffers instead of the original vertex buffers we created for the wireframe mesh by
		// calling CreateRenderData(). This way we can feed the simulated vertex data directly to the renderer
		// without any CPU readback or copying.
		object->cloth->UpdateStreamoutGPU(cmd);
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

	// Conservative AABB — cloth can move anywhere within this large volume.
	// Without position readback we can't compute exact AABB, so use generous bounds.
	//
	// TODO: For tighter AABB without full readback, implement a parallel reduction
	// compute shader that computes min/max of all positions in posBuffer and writes
	// the result to a small 24-byte buffer (6 floats). Then readback only those
	// 24 bytes and set mesh.aabb from them — orders of magnitude cheaper than
	// reading back all particle positions.
	float halfExtent = cloth.params.spacing * std::max(cloth.numX, cloth.numZ) * 0.5f + 2.0f;
	mesh.aabb = wi::primitive::AABB(
		XMFLOAT3(-halfExtent, -1.0f, -halfExtent),
		XMFLOAT3( halfExtent, cloth.params.clothY + 2.0f, halfExtent));

    ObjectComponent& obj = scene.objects.Create(entity);
    obj.meshID = entity;
    obj.SetRenderable(true);
    obj.SetCastShadow(false);

    LayerComponent& layer = scene.layers.Create(entity);
    layer.layerMask = 1;
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
