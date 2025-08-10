#include "sample.h"

void SampleRenderPath::Load()
{
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
	wi::scene::Scene& scene = wi::scene::GetScene();

	auto materialEntity = scene.Entity_CreateMaterial("default");
	auto& material = *scene.materials.GetComponent(materialEntity);
	material.shaderType = wi::scene::MaterialComponent::SHADERTYPE_UNLIT;
	material.SetUseVertexColors(true);
	material.SetDoubleSided(true);

	auto triangleMeshEntity = scene.Entity_CreateMesh("triangle");

	auto& objectTriangle = scene.objects.Create(triangleMeshEntity);
	objectTriangle.meshID = triangleMeshEntity;

	auto& meshTriangle = *scene.meshes.GetComponent(triangleMeshEntity);
	meshTriangle.subsets.push_back(wi::scene::MeshComponent::MeshSubset());
	meshTriangle.subsets.back().materialID = materialEntity;
	meshTriangle.indices.resize(3);
	meshTriangle.subsets.back().indexOffset = (uint32_t)0;
	meshTriangle.subsets.back().indexCount = (uint32_t)3;
	meshTriangle.indices[0] = 0;
	meshTriangle.indices[1] = 2;
	meshTriangle.indices[2] = 1;
	meshTriangle.vertex_positions.resize(3);
	meshTriangle.vertex_positions[0] = XMFLOAT3(-1.0f, -0.5f, 0.0f);
	meshTriangle.vertex_positions[1] = XMFLOAT3(1.0f, -0.5f, 0.0f);
	meshTriangle.vertex_positions[2] = XMFLOAT3(0.0f, 1.0f, 0.0f);
	meshTriangle.vertex_colors.resize(3);
	meshTriangle.vertex_colors[0] = wi::math::CompressColor(XMFLOAT4(1.0f, 0.0f, 0.0f, 1.0f));  // rgba Red;
	meshTriangle.vertex_colors[1] = wi::math::CompressColor(XMFLOAT4(0.0f, 1.0f, 0.0f, 1.0f));  // rgba Green;
	meshTriangle.vertex_colors[2] = wi::math::CompressColor(XMFLOAT4(0.0f, 0.0f, 1.0f, 1.0f));  // rgba Blue;

	// Create vertex and index buffers in a single general buffer accessible by the GPU.
	// Also create an SRV for each buffer type (indices, positions, colors, etc) in the general buffer
	// and put them into a heap (bindless if slots are available) for shader access.
	meshTriangle.CreateRenderData();

	scene.transforms.Create(triangleMeshEntity);
	auto meshtrans = scene.transforms.GetComponent(triangleMeshEntity);
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

	auto triangleMeshEntity = wi::scene::GetScene().Entity_FindByName("triangle");
    TransformComponent* transform = GetScene().transforms.GetComponent(triangleMeshEntity);

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

