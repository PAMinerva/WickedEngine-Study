#include <wiColor.h>
#include <wiECS.h>
#include <wiEnums.h>
#include <wiGUI.h>
#include <wiMath.h>
#include <wiRenderer.h>
#include <wiGraphics.h>
#include "stdafx.h"
#include "cloth.h"
#include "simulation_utils.h"

// Camera position, orientation (40 degrees down) and speed:
#define CAMERAMOVESPEED 10.0f
static float camera_pos[3] = {0.0f, 4.0f, -4.0f};
static float camera_ang[3] = {40.0f, 0.0f, 0.0f};
wi::scene::TransformComponent camera_transform;

static uint64_t wire_shader_id = 0;
static bool mouse_down = false;

extern wi::gui::Label label_tris;
extern wi::gui::Label label_verts;

void init_wicked_scene()
{
    wi::scene::Scene& scene = wi::scene::GetScene();

    auto& weather = scene.weathers.Create(wi::ecs::CreateEntity());
    weather.ambient = XMFLOAT3(0.313f, 0.313f, 0.313f);

    // PointLight1
    {
        auto lightEntity = scene.Entity_CreateLight("PointLight1", XMFLOAT3(-2, 4, -2));
        auto& light = *scene.lights.GetComponent(lightEntity);
        light.SetType(wi::scene::LightComponent::LightType::POINT);
        light.intensity = 5.0f;
        light.SetCastShadow(true);
    }

    // PointLight2
    {
        auto lightEntity = scene.Entity_CreateLight("PointLight2", XMFLOAT3(2, 4, -2));
        auto& light = *scene.lights.GetComponent(lightEntity);
        light.SetType(wi::scene::LightComponent::LightType::POINT);
        light.intensity = 5.0f;
        light.SetCastShadow(true);
    }

    // Ground
    {
        auto groundEntity = scene.Entity_CreatePlane("ground_entity");
        auto groundComponent = scene.transforms.GetComponent(groundEntity);
        groundComponent->Scale(XMFLOAT3(10, 0, 10));
    }

    // Sphere for collision
    {
        auto sphereEntity = scene.Entity_CreateSphere("collision_sphere", 1.0f);
        auto sphereTransform = scene.transforms.GetComponent(sphereEntity);
        sphereTransform->Translate(XMFLOAT3(0.0f, 1.5f, 0.0f));
        sphereTransform->Scale(XMFLOAT3(0.5f, 0.5f, 0.5f));
        sphereTransform->UpdateTransform();

        auto meshComp = scene.meshes.GetComponent(sphereEntity);
        if (meshComp && meshComp->subsets.size() > 0)
        {
            auto matComp = scene.materials.GetComponent(meshComp->subsets[0].materialID);
            if (matComp)
            {
                matComp->baseColor = XMFLOAT4(0.8f, 0.8f, 0.8f, 1.0f);
            }
        }
    }

    // Grid
    {
        wi::renderer::SetToDrawGridHelper(true);
    }
}

void create_custom_shader()
{
    using namespace wi::graphics;
    using namespace wi::renderer;

    GraphicsDevice* device = wi::graphics::GetDevice();

    PipelineStateDesc desc;
    PipelineState pso;

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
	// MaterialComponent::SHADERTYPE_PBR (see create_wire_mesh) needs
	// PSTYPE_OBJECT_PERMUTATION_BEGIN to work, see LoadShaders in wiRenderer.cpp for details
    desc.vs = GetShader(wi::enums::VSTYPE_OBJECT_COMMON);
    desc.ps = GetShader(wi::enums::PSTYPE_OBJECT_PERMUTATION_BEGIN);
    desc.dss = GetDepthStencilState(wi::enums::DSSTYPE_DEPTHREADEQUAL);
    desc.bs = GetBlendState(wi::enums::BSTYPE_OPAQUE);
    device->CreatePipelineState(&desc, &pso);
    wireShader.pso[wi::enums::RENDERPASS_MAIN] = pso;

    wire_shader_id = RegisterCustomShader(wireShader);
}

wi::gui::Label label;
wi::gui::Button run;

void init_gui(SampleRenderPath& srp)
{
    wi::gui::GUI& gui = srp.GetGUI();

    label.Create("Label1");
    label.SetText("Ten Minute Physics: 16 - GPU Cloth Simulation");
    label.font.params.h_align = wi::font::WIFALIGN_CENTER;
    label.SetSize(XMFLOAT2(340, 20));
    gui.AddWidget(&label);

    static wi::gui::Label label2;
    static wi::gui::Button restart;
    static wi::gui::Button solverToggle;

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
            for (auto& object : gPhysicsScene.objects)
            {
                auto wireMeshComponent = wi::scene::GetScene().meshes.GetComponent(object->wireEntity);
                if (wireMeshComponent)
                    wireMeshComponent->DeleteRenderData();

                simulation::remove_simulation_object(*object);
            }

            gPhysicsScene.objects.clear();
            simulation::init_simulation(wire_shader_id);
        });
    gui.AddWidget(&restart);

    solverToggle.Create("solverToggle");
    solverToggle.SetText("Solver: Coloring");
    solverToggle.SetSize(XMFLOAT2(200, 20));
    solverToggle.SetPos(XMFLOAT2(90, 200));
    solverToggle.OnClick(
        [](wi::gui::EventArgs args)
        {
            gPhysicsScene.solveType = (gPhysicsScene.solveType == 0) ? 1 : 0;
            if (gPhysicsScene.solveType == 0)
                solverToggle.SetText("Solver: Coloring");
            else
                solverToggle.SetText("Solver: Jacobi");
        });
    gui.AddWidget(&solverToggle);

    label_tris.Create("Label4");
    label_tris.font.params.h_align = wi::font::WIFALIGN_CENTER;
    label_tris.SetSize(XMFLOAT2(200, 20));
    label_tris.SetPos(XMFLOAT2(90, 240));
    label_tris.SetColor(wi::Color(0, 120, 0));
    gui.AddWidget(&label_tris);

    label_verts.Create("Label5");
    label_verts.font.params.h_align = wi::font::WIFALIGN_CENTER;
    label_verts.SetSize(XMFLOAT2(200, 20));
    label_verts.SetPos(XMFLOAT2(90, 270));
    label_verts.SetColor(wi::Color(0, 120, 0));
    gui.AddWidget(&label_verts);

    label2.Create("Label2");
    label2.SetText("WASD - Move Camera\n"
                   "[R]Mouse - Rotate Camera\n"
                   "[L]Mouse - Grab Cloth");
    label2.font.params.h_align = wi::font::WIFALIGN_CENTER;
    label2.SetSize(XMFLOAT2(200, 60));
    label2.SetPos(XMFLOAT2(90, 310));
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
    RenderPath3D::Load();
}

void SampleRenderPath::FixedUpdate()
{
    // GPU simulation is done in Update()
}

void SampleRenderPath::Update(float dt)
{
    static wi::jobsystem::context init_ctx;
    static bool isCustomShaderReady = false;
    static bool initializationStarted = false;

    if (!initializationStarted)
    {
        initializationStarted = true;
        wi::jobsystem::Execute(init_ctx,
            [](wi::jobsystem::JobArgs args)
            {
                create_custom_shader();
                simulation::init_simulation(wire_shader_id);
                isCustomShaderReady = true;
            });
    }

    if (!isCustomShaderReady)
        return;

    wi::graphics::GraphicsDevice* device = wi::graphics::GetDevice();
    wi::graphics::CommandList cmd = device->BeginCommandList();

    // --- PHYSICS SIMULATION ---
    simulation::simulate(cmd, dt);

    // --- GRAB LOGIC ---
    auto pointer = wi::input::GetPointer();
    int mouseX = (int)pointer.x;
    int mouseY = (int)pointer.y;

    if (wi::input::Press(wi::input::MOUSE_BUTTON_LEFT))
    {
        wi::primitive::Ray ray = wi::renderer::GetPickRay(mouseX, mouseY,
            wi::scene::GetCamera().canvas, wi::scene::GetCamera());

        XMFLOAT3 rayOrigin(ray.origin.x, ray.origin.y, ray.origin.z);
        XMFLOAT3 rayDir(ray.direction.x, ray.direction.y, ray.direction.z);

        for (auto& object : gPhysicsScene.objects)
        {
            if (object->cloth->StartGrabGPU(rayOrigin, rayDir))
            {
                mouse_down = true;
                gPhysicsScene.paused = false;
                run.SetText("Stop Simulation");
                break;
            }
        }
    }

    if (mouse_down && wi::input::Down(wi::input::MOUSE_BUTTON_LEFT))
    {
        wi::primitive::Ray ray = wi::renderer::GetPickRay(
            mouseX, mouseY, wi::scene::GetCamera().canvas, wi::scene::GetCamera());

        XMFLOAT3 rayOrigin(ray.origin.x, ray.origin.y, ray.origin.z);
        XMFLOAT3 rayDir(ray.direction.x, ray.direction.y, ray.direction.z);

        for (auto& object : gPhysicsScene.objects)
        {
            object->cloth->DragGPU(rayOrigin, rayDir, cmd);
        }
    }

    if (mouse_down && wi::input::Release(wi::input::MOUSE_BUTTON_LEFT))
    {
        for (auto& object : gPhysicsScene.objects)
        {
            object->cloth->EndGrabGPU(cmd);
        }
        mouse_down = false;
    }

    // --- CAMERA CONTROL ---
    {
        wi::scene::CameraComponent& camera = wi::scene::GetCamera();
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
            float xDif = wi::input::GetMouseState().delta_position.x;
            float yDif = wi::input::GetMouseState().delta_position.y;
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
        camera_transform.Translate(XMFLOAT3(camera_pos[0], camera_pos[1], camera_pos[2]));
        camera_transform.RotateRollPitchYaw(
            XMFLOAT3(wi::math::DegreesToRadians(camera_ang[0]),
                     wi::math::DegreesToRadians(camera_ang[1]),
                     wi::math::DegreesToRadians(camera_ang[2])));
        camera_transform.UpdateTransform();
        camera.TransformCamera(camera_transform);
        camera.UpdateCamera();
    }

    RenderPath3D::Update(dt);
}
