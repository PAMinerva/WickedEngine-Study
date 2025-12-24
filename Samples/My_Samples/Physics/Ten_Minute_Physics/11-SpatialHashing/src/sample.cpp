#include <cstddef>
#include <random>
#include <wiColor.h>
#include <wiECS.h>
#include <wiGUI.h>
#include <wiMath.h>
#include <wiScene.h>
#include <wiScene_Components.h>
#include "stdafx.h"
#include "balls.h"
#include "simulation_utils.h"

#define CAMERAMOVESPEED 10.0f
static float camera_pos[3] = {2.0f, 1.5f, -3.0f};  // Camera starting position
static float camera_ang[3] = {8.0f, -35.0f, 0.0f}; // Camera starting orientation (8 degrees down, 35 degrees left)
wi::scene::TransformComponent camera_transform;

wi::gui::Label label_tets;

void init_wicked_scene()
{
    wi::scene::Scene &scene = wi::scene::GetScene();

    // Ambient light
    auto &weather = scene.weathers.Create(wi::ecs::CreateEntity());
    weather.ambient = XMFLOAT3(0.313f, 0.313f, 0.313f);

    // SpotLight
    {
        auto lightEntity = scene.Entity_CreateLight("SpotLight", XMFLOAT3(2, 2, -2));
        auto &light = *scene.lights.GetComponent(lightEntity);
        light.SetType(wi::scene::LightComponent::LightType::SPOT);
        light.intensity = 20.0f;
        // light.SetCastShadow(true);
        // light.SetVisualizerEnabled(true);

        wi::scene::TransformComponent &transform = *scene.transforms.GetComponent(lightEntity);
        // -70 degree to radians:
        // -70 × (π / 180) ≈ -1,22173
        // -40 degree to radians:
        // -40 × (π / 180) ≈ -0,69813
        // I will round them to -1.2 and -0.7
        transform.RotateRollPitchYaw(XMFLOAT3(-1.2f, -0.7f, 0.0f));
        transform.UpdateTransform();
    }

    // Directional Light
    {
        auto lightEntity = scene.Entity_CreateLight("DirLight", XMFLOAT3(0, 3, 0));
        auto &light = *scene.lights.GetComponent(lightEntity);
        light.SetType(wi::scene::LightComponent::LightType::DIRECTIONAL);
        light.intensity = 1.0f;
        light.color = XMFLOAT3(0.333f, 0.314f, 0.353f);
        // light.SetCastShadow(true);
        // light.SetVisualizerEnabled(true);
    }

    // Ground
    {
        auto groundEntity = scene.Entity_CreatePlane("ground_entity");
        auto groundComponent = scene.transforms.GetComponent(groundEntity);
        groundComponent->Scale(XMFLOAT3(10, 0, 10));
    }

    // Grid
    {
        wi::renderer::SetToDrawGridHelper(true);
    }

    // Camera
    {
        // wi::scene::TransformComponent transform;
        // transform.RotateRollPitchYaw(XMFLOAT3(0.3f, 0.0f, 0.0f));
        // transform.Translate(XMFLOAT3(0.0f, 2.0f, -4.0f));
        // transform.UpdateTransform();
        //
        // auto &cam = wi::scene::GetCamera();
        // cam.TransformCamera(transform);
    }
}

// Create a grid of spheres with random initial velocities
void init_positions_and_velocities()
{
    float radius = 0.025f;
    float spacing = 3.0f * radius;
    float velRand = 0.2f;

    XMFLOAT3 worldMin(-1.0f, 0.0f, -1.0f);
    XMFLOAT3 worldMax(1.0f, 2.0f, 1.0f);

	std::array<XMFLOAT3, 2> worldBounds = { worldMin, worldMax };
 
    int numX = static_cast<int>((worldMax.x - worldMin.x - 2.0f * spacing) / spacing);
    int numY = static_cast<int>((worldMax.y - worldMin. y - 2.0f * spacing) / spacing);
    int numZ = static_cast<int>((worldMax.z - worldMin.z - 2.0f * spacing) / spacing);
 
    std::vector<XMFLOAT3> positions;
    std::vector<XMFLOAT3> velocities;
 
    std::random_device rd;
    std:: mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(-velRand, velRand);
 
    for (int xi = 0; xi < numX; ++xi)
    {
        for (int yi = 0; yi < numY; ++yi)
        {
            for (int zi = 0; zi < numZ; ++zi)
            {
                XMFLOAT3 pos;
                pos.x = worldMin.x + spacing + xi * spacing;
                pos.y = worldMin.y + spacing + yi * spacing;
                pos. z = worldMin.z + spacing + zi * spacing;
                positions.push_back(pos);

                XMFLOAT3 vel;
                vel.x = dis(gen);
                vel.y = dis(gen);
                vel.z = dis(gen);
                velocities. push_back(vel);
            }
        }
    }
 
	label_tets.SetText(std::to_string(positions.size()) + " spheres");
    gPhysicsScene.balls = std::make_unique<Balls>(radius, positions, velocities, worldBounds);
}

wi::gui::Label label;
wi::gui::Button run;
void init_gui(SampleRenderPath &srp)
{
    wi::gui::GUI &gui = srp.GetGUI();

    label.Create("Label1");
    label.SetText("Ten Minute Physics: 11 - Spatial Hashing");
    label.font.params.h_align = wi::font::WIFALIGN_CENTER;
    label.SetSize(XMFLOAT2(340, 20));
    gui.AddWidget(&label);

	static wi::gui::Label label2;
    static wi::gui::Button restart;
	static wi::gui::CheckBox show_collisions;

	run.Create("runButton");
	run.SetText("Run Simulation");
	run.SetSize(XMFLOAT2(200, 20));
	run.SetPos(XMFLOAT2(90, 140));
	run.OnClick(
		[&](wi::gui::EventArgs args) 
		{
			gPhysicsScene.Run();

			if (gPhysicsScene.paused)
				run.SetText("Run Simulation");
			else
				run.SetText("Stop Simulation");
		});
	gui.AddWidget(&run);

	restart.Create("restartButton");
	restart.SetText("Restart Simulation");
	restart.SetSize(XMFLOAT2(200, 20));
	restart.SetPos(XMFLOAT2(90, 170));
	restart.OnClick(
		[&](wi::gui::EventArgs args)
		{
			if (gPhysicsScene.balls->white_mat_entity != wi::ecs::INVALID_ENTITY)
			{	
				wi::scene::GetScene().materials.Remove(gPhysicsScene.balls->white_mat_entity);
				gPhysicsScene.balls->white_mat_entity = wi::ecs::INVALID_ENTITY;
			}

			wi::scene::MeshComponent* mesh_component = wi::scene::GetScene().meshes.GetComponent(gPhysicsScene.balls->mesh_entity);
			if (mesh_component)
				mesh_component->DeleteRenderData();

			for (auto& objectID : gPhysicsScene.balls->sphereEntities)
			{
				wi::scene::GetScene().objects.Remove(objectID);
				wi::scene::GetScene().transforms.Remove(objectID);
				wi::scene::GetScene().Entity_Remove(objectID);
			}

			gPhysicsScene.balls.reset();
			init_positions_and_velocities();
		});
	gui.AddWidget(&restart);

	show_collisions.Create("showCollidersCheckbox");
	show_collisions.SetText("Show Collisions         ");
	show_collisions.SetSize(XMFLOAT2(20, 20));
	show_collisions.SetPos(XMFLOAT2(270, 200));
	show_collisions.OnClick(
		[&](wi::gui::EventArgs args)
		{
			gPhysicsScene.balls->showCollisions = args.bValue;
		});
	gui.AddWidget(&show_collisions);

    label_tets.Create("Label3");
    label_tets.font.params.h_align = wi::font::WIFALIGN_CENTER;
    label_tets.SetSize(XMFLOAT2(200, 20));
	label_tets.SetPos(XMFLOAT2(90, 230));
	label_tets.SetColor(wi::Color(0, 120, 0));
    gui.AddWidget(&label_tets);

    label2.Create("Label2");
    label2.SetText("WASD - Move Camera\n"
				   "[R]Mouse - Rotate Camera\n"
				   );
    label2.font.params.h_align = wi::font::WIFALIGN_CENTER;
    label2.SetSize(XMFLOAT2(200, 50));
	label2.SetPos(XMFLOAT2(90, 260));
    gui.AddWidget(&label2);
}

void SampleRenderPath::ResizeLayout()
{
    RenderPath3D::ResizeLayout();

    float screenW = GetLogicalWidth();
    float screenH = GetLogicalHeight();
    label.SetPos(XMFLOAT2(screenW / 2.f - label.scale.x / 2.f, screenH * 0.95f));
}

void SampleApp::Initialize()
{
	wi::Application::Initialize();

	infoDisplay.active = true;
	infoDisplay.watermark = true;
	infoDisplay.fpsinfo = true;
	infoDisplay.resolution = true;
	infoDisplay.heap_allocation_counter = true;
	infoDisplay.colorspace = true;
	infoDisplay.device_name = true;

	renderer.init(canvas);
	renderer.Load();

	ActivatePath(&renderer);
}

void SampleRenderPath::Load()
{
    init_gui(*this);
    init_wicked_scene();

	init_positions_and_velocities();

    // Here we call the base class Load method in case it has any additional setup to do.
    RenderPath3D::Load();
}

void SampleRenderPath::FixedUpdate()
{
	// simulation::simulate(0);
	//
	// for (auto &object : gPhysicsScene.objects)
	// {
	// 	auto meshComponent = wi::scene::GetScene().meshes.GetComponent(object->entity);
	// 	simulation::update_mesh(*object->softBody, *meshComponent, true);
	// }
}

void SampleRenderPath::Update(float dt)
{
	// --- PHYSICS SIMULATION ---
	simulation::simulate(dt);

    // --- CAMERA CONTROL ---
    {
        wi::scene::CameraComponent &camera = wi::scene::GetCamera();
        // XMFLOAT4 currentMouse = wi::input::GetPointer();
        static XMFLOAT4 originalMouse = XMFLOAT4(0, 0, 0, 0);
        static bool camControlStart = true;
        if (camControlStart)
        {
            originalMouse = wi::input::GetPointer();
        }

        if (wi::input::Down(wi::input::MOUSE_BUTTON_MIDDLE) ||
            wi::input::Down(wi::input::MOUSE_BUTTON_RIGHT))
        {
            camControlStart = false;
            // Mouse delta
            float xDif = wi::input::GetMouseState().delta_position.x;
            float yDif = wi::input::GetMouseState().delta_position.y;
            // float xDif = currentMouse.x - originalMouse.x;
            // float yDif = currentMouse.y - originalMouse.y;
            xDif = 0.1f * xDif;
            yDif = 0.1f * yDif;
            wi::input::SetPointer(originalMouse);
            wi::input::HidePointer(true);
            camera_ang[0] += yDif;
            camera_ang[1] += xDif;

            if (camera_ang[0] < -89.999f)
                camera_ang[0] = -89.999f;

            if (camera_ang[0] > 89.999f)
                camera_ang[0] = 89.999f;
        }
        else
        {
            camControlStart = true;
            wi::input::HidePointer(false);
        }

        float movespeed = CAMERAMOVESPEED;
        if (wi::input::Down(wi::input::KEYBOARD_BUTTON_LSHIFT) ||
            wi::input::Down(wi::input::KEYBOARD_BUTTON_RSHIFT))
        {
            movespeed *= 3.0f;
        }
        movespeed *= dt;

        if (wi::input::Down((wi::input::BUTTON)'W') ||
            wi::input::Down(wi::input::KEYBOARD_BUTTON_UP))
        {
            camera_pos[0] += movespeed * camera.At.x;
            camera_pos[1] += movespeed * camera.At.y;
            camera_pos[2] += movespeed * camera.At.z;
        }
        if (wi::input::Down((wi::input::BUTTON)'S') ||
            wi::input::Down(wi::input::KEYBOARD_BUTTON_DOWN))
        {
            camera_pos[0] -= movespeed * camera.At.x;
            camera_pos[1] -= movespeed * camera.At.y;
            camera_pos[2] -= movespeed * camera.At.z;
        }

        XMFLOAT3 dir_right;
        XMStoreFloat3(&dir_right, camera.GetRight());
        if (wi::input::Down((wi::input::BUTTON)'D') ||
            wi::input::Down(wi::input::KEYBOARD_BUTTON_RIGHT))
        {
            camera_pos[0] += movespeed * dir_right.x;
            camera_pos[1] += movespeed * dir_right.y;
            camera_pos[2] += movespeed * dir_right.z;
        }
        if (wi::input::Down((wi::input::BUTTON)'A') ||
            wi::input::Down(wi::input::KEYBOARD_BUTTON_LEFT))
        {
            camera_pos[0] -= movespeed * dir_right.x;
            camera_pos[1] -= movespeed * dir_right.y;
            camera_pos[2] -= movespeed * dir_right.z;
        }

        camera_transform.ClearTransform();
        camera_transform.Translate(
            XMFLOAT3(camera_pos[0],
					 camera_pos[1],
					 camera_pos[2]));
		camera_transform.RotateRollPitchYaw(
			XMFLOAT3(wi::math::DegreesToRadians(camera_ang[0]),
			         wi::math::DegreesToRadians(camera_ang[1]),
		             wi::math::DegreesToRadians(camera_ang[2]))
		);
        camera_transform.UpdateTransform();
        camera.TransformCamera(camera_transform);
        camera.UpdateCamera();
    }

    // Here we call the base class Update method in case it has any additional update to do.
    RenderPath3D::Update(dt);
}
