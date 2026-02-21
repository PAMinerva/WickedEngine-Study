#include <cstddef>
#include <wiColor.h>
#include <wiECS.h>
#include <wiEnums.h>
#include <wiGUI.h>
#include <wiMath.h>
#include <wiRenderer.h>
#include "stdafx.h"
#include "grabber.h"
#include "cloth.h"
#include "simulation_utils.h"

#define CAMERAMOVESPEED 10.0f
static float camera_pos[3] = {0.5f, 1.0f, -0.8f};	// Camera starting position
static float camera_ang[3] = {40.0f, -30.0f, 0.0f};	// Camera starting orientation (40 degrees down, 30 degrees to the left)
wi::scene::TransformComponent camera_transform;

static uint64_t front_shader_id = 0;
static uint64_t back_shader_id = 0;
static uint64_t wire_shader_id = 0;

static bool mouse_down = false;
static Grabber grabber;

extern wi::gui::Label label_tris;
extern wi::gui::Label label_verts;

void init_wicked_scene()
{
    wi::scene::Scene &scene = wi::scene::GetScene();

    // Ambient light
    auto &weather = scene.weathers.Create(wi::ecs::CreateEntity());
    weather.ambient = XMFLOAT3(0.313f, 0.313f, 0.313f);

	//PointLight1
	{
		auto lightEntity = scene.Entity_CreateLight("PointLight1", XMFLOAT3(-2, 2, -2));
		auto &light = *scene.lights.GetComponent(lightEntity);
		light.SetType(wi::scene::LightComponent::LightType::POINT);
		light.intensity = 5.0f;
		light.SetCastShadow(true);
		// light.SetVisualizerEnabled(true);
	}

	//PointLight2
	{
		auto lightEntity = scene.Entity_CreateLight("PointLight2", XMFLOAT3(2, 2, -2));
		auto &light = *scene.lights.GetComponent(lightEntity);
		light.SetType(wi::scene::LightComponent::LightType::POINT);
		light.intensity = 5.0f;
		light.SetCastShadow(true);
	}

	//SpotLight
	{
		auto lightEntity = scene.Entity_CreateLight("SpotLight", XMFLOAT3(0, 2, -2));
		auto &light = *scene.lights.GetComponent(lightEntity);
		light.SetType(wi::scene::LightComponent::LightType::SPOT);
        light.intensity = 5.0f;
        light.SetCastShadow(true);
        // light.SetVisualizerEnabled(true);

        wi::scene::TransformComponent &transform = *scene.transforms.GetComponent(lightEntity);
        // -40 degree (counter-clockwise around the Z-axis; Roll) to radians:
        // -40 × (π / 180) ≈ -0,69813
        // I will round it to -0.7
        transform.RotateRollPitchYaw(XMFLOAT3(-0.7f, 0.0f, 0.0f));
        transform.UpdateTransform();

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

void create_custom_shader()
{
	using namespace wi::graphics;
    using namespace wi::renderer;

    GraphicsDevice* device = wi::graphics::GetDevice();

    PipelineStateDesc desc;
    PipelineState pso;

    // ========== FRONT FACE SHADER (red) ==========
    static RasterizerState rs_front;
    rs_front.fill_mode = FillMode::SOLID;
    rs_front.cull_mode = CullMode::BACK;      // Cull back faces
    rs_front.front_counter_clockwise = true;  // CCW = front
    rs_front.depth_clip_enable = true;

    CustomShader frontShader;
    frontShader.name = "ClothFrontFace";
    frontShader.filterMask = wi::enums::FILTER_OPAQUE;

    // Prepass
    desc.vs = GetShader(wi::enums::VSTYPE_OBJECT_PREPASS);
    desc.ps = GetShader(wi::enums::PSTYPE_OBJECT_PREPASS);
    desc.dss = GetDepthStencilState(wi::enums::DSSTYPE_DEFAULT);
    desc.rs = &rs_front;
    device->CreatePipelineState(&desc, &pso);
    frontShader.pso[wi::enums::RENDERPASS_PREPASS] = pso;

    // Main Pass
    desc.vs = GetShader(wi::enums::VSTYPE_OBJECT_COMMON);
    desc.ps = GetShader(wi::enums::PSTYPE_OBJECT_PERMUTATION_BEGIN);
    desc.dss = GetDepthStencilState(wi::enums::DSSTYPE_DEPTHREADEQUAL);
    desc.bs = GetBlendState(wi::enums::BSTYPE_OPAQUE);
    device->CreatePipelineState(&desc, &pso);
    frontShader.pso[wi::enums::RENDERPASS_MAIN] = pso;

    // Shadow Pass
    desc.vs = GetShader(wi::enums::VSTYPE_SHADOW);
    desc.ps = nullptr;
    desc.bs = GetBlendState(wi::enums::BSTYPE_OPAQUE);
    desc.rs = GetRasterizerState(wi::enums::RSTYPE_SHADOW);
    desc.dss = GetDepthStencilState(wi::enums::DSSTYPE_SHADOW);
    device->CreatePipelineState(&desc, &pso);
    frontShader.pso[wi::enums::RENDERPASS_SHADOW] = pso;

    // ========== BACK FACE SHADER (yellow-ish) ==========
    static RasterizerState rs_back;
    rs_back.fill_mode = FillMode::SOLID;
    rs_back.cull_mode = CullMode::FRONT; // Cull front faces
    rs_back.front_counter_clockwise = true;
    rs_back.depth_clip_enable = true;

    CustomShader backShader;
    backShader.name = "ClothBackFace";
    backShader.filterMask = wi::enums::FILTER_OPAQUE;

    // Prepass
    desc = PipelineStateDesc{};
    desc.vs = GetShader(wi::enums::VSTYPE_OBJECT_PREPASS);
    desc.ps = GetShader(wi::enums::PSTYPE_OBJECT_PREPASS);
    desc.dss = GetDepthStencilState(wi::enums::DSSTYPE_DEFAULT);
    desc.rs = &rs_back;
    device->CreatePipelineState(&desc, &pso);
    backShader.pso[wi::enums::RENDERPASS_PREPASS] = pso;

    // Main Pass
    desc.vs = GetShader(wi::enums::VSTYPE_OBJECT_COMMON);
    desc.ps = GetShader(wi::enums::PSTYPE_OBJECT_PERMUTATION_BEGIN);
    desc.dss = GetDepthStencilState(wi::enums::DSSTYPE_DEPTHREADEQUAL);
    desc.bs = GetBlendState(wi::enums::BSTYPE_OPAQUE);
    device->CreatePipelineState(&desc, &pso);
    backShader.pso[wi::enums::RENDERPASS_MAIN] = pso;

    // Shadow Pass
    desc.vs = GetShader(wi::enums::VSTYPE_SHADOW);
    desc.ps = nullptr;
    desc.bs = GetBlendState(wi::enums::BSTYPE_OPAQUE);
    desc.rs = GetRasterizerState(wi::enums::RSTYPE_SHADOW);
    desc.dss = GetDepthStencilState(wi::enums::DSSTYPE_SHADOW);
    device->CreatePipelineState(&desc, &pso);
    backShader.pso[wi::enums::RENDERPASS_SHADOW] = pso;

    // ========== WIREFRAME SHADER (opzionale, per debug) ==========
    static RasterizerState rs_wire;
    rs_wire.fill_mode = FillMode::WIREFRAME;
    rs_wire.cull_mode = CullMode::NONE;
    rs_wire.front_counter_clockwise = true;

    CustomShader wireShader;
    wireShader.name = "ClothWireframe";
    wireShader.filterMask = wi::enums::FILTER_OPAQUE;

    // Prepass
    desc = PipelineStateDesc{};
    desc.vs = GetShader(wi::enums::VSTYPE_OBJECT_PREPASS);
    desc.ps = GetShader(wi::enums::PSTYPE_OBJECT_PREPASS);
    desc.dss = GetDepthStencilState(wi::enums::DSSTYPE_DEFAULT);
    desc.rs = &rs_wire;
    device->CreatePipelineState(&desc, &pso);
    wireShader.pso[wi::enums::RENDERPASS_PREPASS] = pso;

    // Main Pass
    desc.vs = GetShader(wi::enums::VSTYPE_OBJECT_COMMON);
    desc.ps = GetShader(static_cast<wi::enums::SHADERTYPE>(
        wi::enums::PSTYPE_OBJECT_PERMUTATION_BEGIN +
        wi::scene::MaterialComponent::SHADERTYPE_UNLIT));
    desc.dss = GetDepthStencilState(wi::enums::DSSTYPE_DEPTHREADEQUAL);
    desc.bs = GetBlendState(wi::enums::BSTYPE_OPAQUE);
    device->CreatePipelineState(&desc, &pso);
    wireShader.pso[wi::enums::RENDERPASS_MAIN] = pso;

	// Register custom shaders
    front_shader_id = RegisterCustomShader(frontShader);
    back_shader_id = RegisterCustomShader(backShader);
	wire_shader_id = RegisterCustomShader(wireShader);
}

wi::gui::Label label;
wi::gui::Button run;
void init_gui(SampleRenderPath &srp)
{
    wi::gui::GUI &gui = srp.GetGUI();

    label.Create("Label1");
    label.SetText("Ten Minute Physics: 15 - Self Collision");
    label.font.params.h_align = wi::font::WIFALIGN_CENTER;
    label.SetSize(XMFLOAT2(340, 20));
    gui.AddWidget(&label);

	static wi::gui::Label label2;
    static wi::gui::Button restart;
	static wi::gui::CheckBox handle_collisions;
    static wi::gui::Slider compliance;
	static wi::gui::CheckBox show_wireMesh;

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
			for (auto &object : gPhysicsScene.objects)
			{
				auto frontMeshComponent = wi::scene::GetScene().meshes.GetComponent(object->frontEntity);
				frontMeshComponent->DeleteRenderData();

				auto backMeshComponent = wi::scene::GetScene().meshes.GetComponent(object->backEntity);
				backMeshComponent->DeleteRenderData();

				auto wireMeshComponent = wi::scene::GetScene().meshes.GetComponent(object->wireEntity);
				wireMeshComponent->DeleteRenderData();

				simulation::remove_simulation_object(*object);
			}

			gPhysicsScene.objects.clear();
			simulation::init_simulation(front_shader_id, back_shader_id, wire_shader_id);
		});
	gui.AddWidget(&restart);

	handle_collisions.Create("handleCollisionsCheckbox");
	handle_collisions.SetText("Handle Self-Collisions    ");
	handle_collisions.SetSize(XMFLOAT2(20, 20));
	handle_collisions.SetPos(XMFLOAT2(270, 200));
	handle_collisions.SetCheck(true);
	handle_collisions.OnClick(
		[&](wi::gui::EventArgs args)
		{
			for (auto& object : gPhysicsScene.objects)
			{
				object->cloth->handleCollisions = args.bValue;
			}
		});
	gui.AddWidget(&handle_collisions);

	show_wireMesh.Create("showWireMeshCheckbox");
	show_wireMesh.SetText("Show Wireframe Mesh    ");
	show_wireMesh.SetSize(XMFLOAT2(20, 20));
	show_wireMesh.SetPos(XMFLOAT2(270, 230));
	show_wireMesh.OnClick(
		[&](wi::gui::EventArgs args)
		{
			for (size_t i = 0; i < gPhysicsScene.objects.size(); i++)
			{
				// Show/hide wireframe
				auto objWireComponent = wi::scene::GetScene().objects.GetComponent(gPhysicsScene.objects[i]->wireEntity);
				objWireComponent->SetRenderable(args.bValue);

				// Show/hide front and back faces (opposite of wireframe)
				auto objFrontComponent = wi::scene::GetScene().objects.GetComponent(gPhysicsScene.objects[i]->frontEntity);
				objFrontComponent->SetRenderable(!args.bValue);

				auto objBackComponent = wi::scene::GetScene().objects.GetComponent(gPhysicsScene.objects[i]->backEntity);
				objBackComponent->SetRenderable(!args.bValue);
			}
		});
	gui.AddWidget(&show_wireMesh);

	compliance.Create(0, 10, 1, 10, "slider1");
	compliance.SetText("Compliance: ");
	compliance.SetSize(XMFLOAT2(200, 20));
	compliance.SetPos(XMFLOAT2(90, 260));
	compliance.SetValue(1);
	compliance.OnSlide(
		[](wi::gui::EventArgs args)
		{
			for (size_t i = 0; i < gPhysicsScene.objects.size(); i++)
			{
				// Update bending compliance for all bending constraints
				// In 15-selfCollision, bending constraints are the last two types (indices 4 and 5)
				// but since we store compliance per-constraint, we need to update them all
				// or track which constraints are bending constraints

				// For now, update the stored bendingCompliance value
				// The actual per-constraint compliances were set at construction time
				// To properly update, we'd need to regenerate constraints or track bending indices

				// Simple approach: store a bendingCompliance member and use it in SolveConstraints
				gPhysicsScene.objects[i]->cloth->bendingCompliance = args.fValue;
			}
		});
	gui.AddWidget(&compliance);

	// label_tets.Create("Label3");
	// label_tets.font.params.h_align = wi::font::WIFALIGN_CENTER;
	// label_tets.SetSize(XMFLOAT2(200, 20));
	// label_tets.SetPos(XMFLOAT2(90, 290));
	// label_tets.SetColor(wi::Color(0, 120, 0));
	// gui.AddWidget(&label_tets);

    label_tris.Create("Label4");
    label_tris.font.params.h_align = wi::font::WIFALIGN_CENTER;
    label_tris.SetSize(XMFLOAT2(200, 20));
	label_tris.SetPos(XMFLOAT2(90, 290));
	label_tris.SetColor(wi::Color(0, 120, 0));
    gui.AddWidget(&label_tris);

    label_verts.Create("Label5");
    label_verts.font.params.h_align = wi::font::WIFALIGN_CENTER;
    label_verts.SetSize(XMFLOAT2(200, 20));
	label_verts.SetPos(XMFLOAT2(90, 320));
	label_verts.SetColor(wi::Color(0, 120, 0));
    gui.AddWidget(&label_verts);

    label2.Create("Label2");
    label2.SetText("WASD - Move Camera\n"
				   "[R]Mouse - Rotate Camera\n"
				   "[L]Mouse - Pick Object");
    label2.font.params.h_align = wi::font::WIFALIGN_CENTER;
    label2.SetSize(XMFLOAT2(200, 60));
	label2.SetPos(XMFLOAT2(90, 350));
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
	static wi::jobsystem::context init_ctx;
    static bool isCustomShaderReady = false;
    static bool initializationStarted = false;

	// Create custom shader and dependent physics simulation initialization in a background job
	// to avoid stalling the main thread.
	// This is because shader compilation may take some time, especially for the shaders in the
	// [PSTYPE_OBJECT_PERMUTATION_BEGIN, PSTYPE_OBJECT_PERMUTATION_END] range and also in the
	// [PSTYPE_OBJECT_TRANSPARENT_PERMUTATION_BEGIN, PSTYPE_OBJECT_TRANSPARENT_PERMUTATION_END] range.
	// PSOs using these shaders may take significant time to compile in GPU machine code so the PSO creation
	// and relative shader comppilation in GPU machine code continues in background even after engine was
	// initialized to avoid stalls.
    if (!initializationStarted)
    {
        initializationStarted = true;

        wi::jobsystem::Execute(init_ctx,
							   [](wi::jobsystem::JobArgs args)
								{
									create_custom_shader();
									simulation::init_simulation(front_shader_id, back_shader_id, wire_shader_id);
									isCustomShaderReady = true;
								});
    }

    if (!isCustomShaderReady)
		return;

	// --- PHYSICS SIMULATION ---
	simulation::simulate(dt);

	// Update all three meshes for each simulation object
	for (auto& object : gPhysicsScene.objects)
	{
		auto frontMesh = wi::scene::GetScene().meshes.GetComponent(object->frontEntity);
		simulation::update_mesh(*object->cloth, *frontMesh, true);

		auto backMesh = wi::scene::GetScene().meshes.GetComponent(object->backEntity);
		simulation::update_mesh(*object->cloth, *backMesh, true);

		// Wireframe only if visible
		auto wireObj = wi::scene::GetScene().objects.GetComponent(object->wireEntity);
		if (wireObj && wireObj->IsRenderable())
		{
			auto wireMesh = wi::scene::GetScene().meshes.GetComponent(object->wireEntity);
			simulation::update_mesh(*object->cloth, *wireMesh, true);
		}
	}

    // --- GRAB LOGIC ---
    auto pointer = wi::input::GetPointer();
    int mouseX = (int)pointer.x;
    int mouseY = (int)pointer.y;

    // Mouse PRESS: start grab
    if (wi::input::Press(wi::input::MOUSE_BUTTON_LEFT))
    {
        wi::primitive::Ray ray = wi::renderer::GetPickRay(mouseX, mouseY,
                                                          wi::scene::GetCamera().canvas,
                                                          wi::scene::GetCamera());
        wi::scene::PickResult pick = wi::scene::Pick(ray);

        for (auto& object : gPhysicsScene.objects)
        {
            if (pick.entity == object->frontEntity ||
                pick.entity == object->backEntity /*||
                pick.entity == object->wireEntity*/)
            {
                grabber.start(pick, object->cloth.get());
                mouse_down = true;
                gPhysicsScene.paused = false;
                run.SetText("Stop Simulation");
                break;
            }
        }
    }

    // Mouse MOVE: drag object
    if (mouse_down && wi::input::Down(wi::input::MOUSE_BUTTON_LEFT))
    {
        wi::primitive::Ray ray = wi::renderer::GetPickRay(
            mouseX, mouseY, wi::scene::GetCamera().canvas,
            wi::scene::GetCamera());
        wi::scene::PickResult pick = wi::scene::Pick(ray);

        grabber.move(ray, pick);
    }

    // Mouse RELEASE: end grab
    if (mouse_down && wi::input::Release(wi::input::MOUSE_BUTTON_LEFT))
    {
        grabber.end();
        mouse_down = false;
    }

    grabber.increaseTime(dt);

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
