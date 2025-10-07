#include "sample.h"
#include <wiECS.h>
#include <wiScene_Components.h>

static wi::ecs::Entity entity;

void SampleRenderPath::Load()
{
	using namespace wi::scene;
	using namespace wi::ecs;

	// In an Entity-Component-System (ECS) architecture, entities and components are fundamental concepts:
	// 
	// Entity:
	// An entity is a general-purpose object that represents a unique identifier. It does not contain any data or behavior by itself.
	// Instead, it acts as an identifier for components. Entities are often implemented as simple integer IDs.
	// In the context of a game, an entity could represent a player, an enemy, a projectile, or any other object in the game world.
	// 
	// Component:
	// A component is a modular piece of data that can be attached to an entity. Each component typically contains specific data
	// related to a particular aspect of the entity, such as its position, velocity, health, or appearance. Usually, components
	// do not contain behavior (methods) and are implemented as plain data structures (structs or classes with only data members).
	// By attaching different combinations of components to entities, you can create complex behaviors and characteristics.
	// However, components can contain methods that set or initialize their data, but these methods should not contain complex logic.
	// 
	// System:
	// A system is responsible for processing entities that have a specific set of components. Systems contain the logic and behavior
	// that operate on the data contained in components. For example, a physics system might update the positions of all entities
	// that have both position and velocity components, while a rendering system might draw all entities that have a mesh component.
	// 
	// In summary, ECS is a design pattern that promotes separation of concerns by dividing data (components) and behavior (systems),
	// with entities acting as unique identifiers that group related components together.
	Scene& scene = GetScene();

	// Create a material for the mesh:
	// We specify that the material won't be affected by lighting (SHADERTYPE_UNLIT),
	// and that it will use vertex colors (SetUseVertexColors(true)) to color the mesh.
	Entity mat_entity = CreateEntity();
	scene.materials.Create(mat_entity);
	MaterialComponent &material = *scene.materials.GetComponent(mat_entity);
	material.shaderType = wi::scene::MaterialComponent::SHADERTYPE_UNLIT;
	material.SetUseVertexColors(true);

	entity = wi::ecs::CreateEntity();
	scene.layers.Create(entity);
	scene.transforms.Create(entity);
	ObjectComponent &object = scene.objects.Create(entity);

    // We store the mesh entity in meshID for later access.
    // Without this, there would be no way to determine which mesh an object is using,
    // as the association between objects and meshes is only maintained by
    // ComponentManager<ObjectComponent>, which is not directly accessible
    // from an ObjectComponent instance.
	auto &mesh = scene.meshes.Create(entity);
	object.meshID = entity;

	mesh.subsets.push_back(wi::scene::MeshComponent::MeshSubset());
	mesh.subsets.back().materialID = mat_entity;
	mesh.indices.resize(3);
	mesh.subsets.back().indexOffset = 0;
	mesh.subsets.back().indexCount = 3;
	mesh.indices[0] = 0;
	mesh.indices[1] = 2;
	mesh.indices[2] = 1;
	mesh.vertex_positions.resize(3);
	mesh.vertex_positions[0] = XMFLOAT3(-1.0f, -0.5f, 0.0f);
	mesh.vertex_positions[1] = XMFLOAT3(0.0f, 1.0f, 0.0f);
	mesh.vertex_positions[2] = XMFLOAT3(1.0f, -0.5f, 0.0f);
	mesh.vertex_colors.resize(3);
	mesh.vertex_colors[0] = wi::math::CompressColor(XMFLOAT4(1.0f, 0.0f, 0.0f, 1.0f));  // rgba Red;
	mesh.vertex_colors[1] = wi::math::CompressColor(XMFLOAT4(0.0f, 1.0f, 0.0f, 1.0f));  // rgba Green;
	mesh.vertex_colors[2] = wi::math::CompressColor(XMFLOAT4(0.0f, 0.0f, 1.0f, 1.0f));  // rgba Blue;

	// Create vertex and index buffers as aliasing resources in a single general buffer accessible by the GPU.
	// Also create an SRV for each buffer type (indices, positions, colors, etc) in the general buffer
	// and put them into a heap (bindless if slots are available in the shader-visible descriptor heap
	// associated with the command list, otherwise put them into a local descriptor heap accessible
	// through SingleDescriptor::allocationhandler; we can do it via the SingleDescriptor object
	// created in CreateRenderData (specifically, in the CreateSubresource method called at the end of
	// CreateRenderData) and stored as an element of the subresources_srv vector of the 
	// Resource_DX12 represeting the general buffer resource).
	mesh.CreateRenderData();

	auto meshtrans = scene.transforms.GetComponent(entity);
	meshtrans->Translate(XMFLOAT3(0.0f, -0.2f, 2.0f));
	meshtrans->Scale(XMFLOAT3(1, 1, 1));

	// Here we call the base class Load method in case it has any additional setup to do.
	RenderPath3D::Load();
}

void SampleRenderPath::Update(float dt)
{
    using namespace wi::scene;
    static float speed = 1.5f;
    static float direction = 1.0f;
    static float delta_x = 0.0f;

    TransformComponent* transform = GetScene().transforms.GetComponent(entity);

    if (transform != nullptr)
    {
        float pos_x = transform->GetPosition().x;
        delta_x = speed * direction * dt;

        if (pos_x > 1.25f)
        {
            direction = -1.0f;
        }
        else if (pos_x < -1.25f)
        {
            direction = 1.0f;
        }

        transform->Translate(XMFLOAT3(delta_x, 0, 0));
    }

	// Here we call the base class Update method in case it has any additional setup to do.
    RenderPath3D::Update(dt);
}

