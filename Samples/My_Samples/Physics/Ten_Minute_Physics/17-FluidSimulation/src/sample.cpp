#include <wiColor.h>
#include <wiECS.h>
#include <wiEnums.h>
#include <wiGUI.h>
#include <wiMath.h>
#include <wiRenderer.h>
#include <wiGraphics.h>
#include <wiImage.h>
#include <wiSprite.h>
#include <wiPrimitive.h>
#include "stdafx.h"
#include "fluid_sim.h"

// ============================================================================
// Globals
// ============================================================================
static FluidSimulation3D gFluid;

// Camera
#define CAMERAMOVESPEED 6.0f
static float camera_pos[3] = {0.0f, 3.0f, -8.0f};
static float camera_ang[3] = {20.0f, 0.0f, 0.0f};
wi::scene::TransformComponent camera_transform;

// Obstacle/emitter dragging (left-click, constant depth)
enum class DragTarget { NONE, OBSTACLE, SMOKE_SRC };
static DragTarget dragTarget = DragTarget::NONE;
static float dragDepth = 0.0f;
static XMFLOAT3 lastDragCenter = {0, 0, 0};

// GUI widgets
static wi::gui::Label titleLabel;
static wi::gui::Button runButton;
static wi::gui::Button resetButton;
static wi::gui::ComboBox renderModeCombo;
static wi::gui::CheckBox macCormackCheckbox;
static wi::gui::Slider iterSlider;
static wi::gui::Slider sorSlider;
static wi::gui::Slider absorptionSlider;
static wi::gui::Slider buoyancySlider;
static wi::gui::Slider dissipationSlider;
static wi::gui::Slider vorticitySlider;
static wi::gui::Label infoLabel;

// ============================================================================
// Scene Setup
// ============================================================================

void init_wicked_scene()
{
    wi::scene::Scene& scene = wi::scene::GetScene();

    auto& weather = scene.weathers.Create(wi::ecs::CreateEntity());
    weather.ambient = XMFLOAT3(0.08f, 0.08f, 0.12f);

    // Directional light
    {
        auto lightEntity = scene.Entity_CreateLight("DirLight", XMFLOAT3(0, 5, 0));
        auto& light = *scene.lights.GetComponent(lightEntity);
        light.SetType(wi::scene::LightComponent::LightType::DIRECTIONAL);
        light.intensity = 3.0f;
        light.SetCastShadow(true);

        auto* transform = scene.transforms.GetComponent(lightEntity);
        transform->RotateRollPitchYaw(XMFLOAT3(
            wi::math::DegreesToRadians(60.0f),
            wi::math::DegreesToRadians(-30.0f),
            0.0f));
        transform->UpdateTransform();
    }

    // Point light
    {
        auto lightEntity = scene.Entity_CreateLight("PointLight1", XMFLOAT3(2, 4, -2));
        auto& light = *scene.lights.GetComponent(lightEntity);
        light.SetType(wi::scene::LightComponent::LightType::POINT);
        light.intensity = 4.0f;
    }

    // Ground plane
    {
        auto groundEntity = scene.Entity_CreatePlane("ground");
        auto groundTransform = scene.transforms.GetComponent(groundEntity);
        groundTransform->Scale(XMFLOAT3(10, 0, 10));
    }

    // Grid helper
    wi::renderer::SetToDrawGridHelper(true);
}

// ============================================================================
// GUI
// ============================================================================

void init_gui(SampleRenderPath& srp)
{
    wi::gui::GUI& gui = srp.GetGUI();
    float yPos = 130.0f;
    // Reserve some space on the left because several Wicked Engine widgets
    // (for example sliders, checkboxes, combo boxes) render their caption
    // right-aligned against the widget anchor, so the text appears to the
    // left of the control itself.
    float xPos = 170.0f;
    float w = 220.0f;
    float h = 20.0f;
    float spacing = 30.0f;

    titleLabel.Create("TitleLabel");
    titleLabel.SetText("Ten Minute Physics: 17 - Fluid Simulation (3D GPU)");
    titleLabel.font.params.h_align = wi::font::WIFALIGN_CENTER;
    titleLabel.SetSize(XMFLOAT2(400, 20));
    gui.AddWidget(&titleLabel);

    runButton.Create("RunButton");
    runButton.SetText("Run Simulation");
    runButton.SetSize(XMFLOAT2(w, h));
    runButton.SetPos(XMFLOAT2(xPos, yPos));
    runButton.OnClick([](wi::gui::EventArgs args)
    {
        gFluid.paused = !gFluid.paused;
        runButton.SetText(gFluid.paused ? "Run Simulation" : "Pause Simulation");
    });
    gui.AddWidget(&runButton);
    yPos += spacing;

    resetButton.Create("ResetButton");
    resetButton.SetText("Reset");
    resetButton.SetSize(XMFLOAT2(w, h));
    resetButton.SetPos(XMFLOAT2(xPos, yPos));
    resetButton.OnClick([](wi::gui::EventArgs args)
    {
        gFluid.Reset();
        gFluid.paused = true;
        runButton.SetText("Run Simulation");
    });
    gui.AddWidget(&resetButton);
    yPos += spacing;

    renderModeCombo.Create("RenderModeCombo");
    renderModeCombo.SetText("Visualization  ");
    renderModeCombo.SetSize(XMFLOAT2(w, h));
    renderModeCombo.SetPos(XMFLOAT2(xPos, yPos));
    renderModeCombo.AddItem("Smoke");
    renderModeCombo.AddItem("Pressure");
    renderModeCombo.AddItem("Velocity");
    renderModeCombo.SetSelected(0);
    renderModeCombo.OnSelect([](wi::gui::EventArgs args)
    {
        gFluid.renderMode = args.iValue;
    });
    gui.AddWidget(&renderModeCombo);
    yPos += spacing;

    macCormackCheckbox.Create("MacCormackCB");
    macCormackCheckbox.SetText("MacCormack Advection  ");
    macCormackCheckbox.SetSize(XMFLOAT2(h, h));
    macCormackCheckbox.SetPos(XMFLOAT2(xPos, yPos));
    macCormackCheckbox.SetCheck(gFluid.useMacCormack);
    macCormackCheckbox.OnClick([](wi::gui::EventArgs args)
    {
        gFluid.useMacCormack = args.bValue;
    });
    gui.AddWidget(&macCormackCheckbox);
    yPos += spacing;

    iterSlider.Create(5, 100, 40, 95, "Pressure Iters  ");
    iterSlider.SetSize(XMFLOAT2(w, h));
    iterSlider.SetPos(XMFLOAT2(xPos, yPos));
    iterSlider.SetValue((float)gFluid.params.numPressureIters);
    iterSlider.OnSlide([](wi::gui::EventArgs args)
    {
        gFluid.params.numPressureIters = (int)args.fValue;
    });
    gui.AddWidget(&iterSlider);
    yPos += spacing;

    sorSlider.Create(1.0f, 1.99f, 1.9f, 99, "SOR Omega  ");
    sorSlider.SetSize(XMFLOAT2(w, h));
    sorSlider.SetPos(XMFLOAT2(xPos, yPos));
    sorSlider.SetValue(gFluid.params.overRelaxation);
    sorSlider.OnSlide([](wi::gui::EventArgs args)
    {
        gFluid.params.overRelaxation = args.fValue;
    });
    gui.AddWidget(&sorSlider);
    yPos += spacing;

    absorptionSlider.Create(1.0f, 200.0f, 80.0f, 199, "Smoke Density  ");
    absorptionSlider.SetSize(XMFLOAT2(w, h));
    absorptionSlider.SetPos(XMFLOAT2(xPos, yPos));
    absorptionSlider.SetValue(gFluid.smokeAbsorption);
    absorptionSlider.OnSlide([](wi::gui::EventArgs args)
    {
        gFluid.smokeAbsorption = args.fValue;
    });
    gui.AddWidget(&absorptionSlider);
    yPos += spacing;

    buoyancySlider.Create(0.0f, 5.0f, 1.5f, 200, "Buoyancy  ");
    buoyancySlider.SetSize(XMFLOAT2(w, h));
    buoyancySlider.SetPos(XMFLOAT2(xPos, yPos));
    buoyancySlider.SetValue(gFluid.params.buoyancy);
    buoyancySlider.OnSlide([](wi::gui::EventArgs args)
    {
        gFluid.params.buoyancy = args.fValue;
    });
    gui.AddWidget(&buoyancySlider);
    yPos += spacing;

    dissipationSlider.Create(0.90f, 1.0f, 0.995f, 100, "Smoke Retain  ");
    dissipationSlider.SetSize(XMFLOAT2(w, h));
    dissipationSlider.SetPos(XMFLOAT2(xPos, yPos));
    dissipationSlider.SetValue(gFluid.params.dissipation);
    dissipationSlider.OnSlide([](wi::gui::EventArgs args)
    {
        gFluid.params.dissipation = args.fValue;
    });
    gui.AddWidget(&dissipationSlider);
    yPos += spacing;

    vorticitySlider.Create(0.0f, 2.5f, 0.5f, 250, "Vorticity  ");
    vorticitySlider.SetSize(XMFLOAT2(w, h));
    vorticitySlider.SetPos(XMFLOAT2(xPos, yPos));
    vorticitySlider.SetValue(gFluid.params.vorticityStrength);
    vorticitySlider.OnSlide([](wi::gui::EventArgs args)
    {
        gFluid.params.vorticityStrength = args.fValue;
    });
    gui.AddWidget(&vorticitySlider);
    yPos += spacing;

    infoLabel.Create("InfoLabel");
    infoLabel.SetText("WASD - Move Camera\n"
                      "[R]Mouse - Rotate Camera\n"
                      "[L]Mouse - Grab Obstacle/Emitter");
    infoLabel.font.params.h_align = wi::font::WIFALIGN_CENTER;
    infoLabel.SetSize(XMFLOAT2(w, 60));
    infoLabel.SetPos(XMFLOAT2(xPos, yPos));
    gui.AddWidget(&infoLabel);
}

// ============================================================================
// SampleRenderPath
// ============================================================================

void SampleRenderPath::ResizeLayout()
{
    RenderPath3D::ResizeLayout();

    float screenW = GetLogicalWidth();
    float screenH = GetLogicalHeight();
    titleLabel.SetPos(XMFLOAT2(screenW / 2.f - titleLabel.scale.x / 2.f, screenH * 0.95f));
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
    // Deferred initialization: wait for engine shaders to finish compiling
    // before loading our custom shaders.
    // NOTE: initialization runs on the main thread to avoid a race condition
    // with the global shader source path — background threads that call
    // SetShaderSourcePath can corrupt on-demand compilation of engine shaders
    // (e.g., objectPS permutations triggered by DrawBox/DrawSphere).
    static bool fluidReady = false;
    static bool initStarted = false;

    if (!initStarted)
    {
        initStarted = true;
        gFluid.Initialize(gFluid.params);
        gFluid.paused = false;
        runButton.SetText("Pause Simulation");
        fluidReady = true;
    }

    // --- FLUID SIMULATION + RENDERING (only when ready) ---
    if (fluidReady)
    {
        wi::graphics::GraphicsDevice* device = wi::graphics::GetDevice();
        wi::graphics::CommandList cmd = device->BeginCommandList();

        gFluid.Simulate(cmd);

        wi::scene::CameraComponent& camera = wi::scene::GetCamera();
        uint32_t screenW = (uint32_t)camera.canvas.GetPhysicalWidth();
        uint32_t screenH = (uint32_t)camera.canvas.GetPhysicalHeight();
        gFluid.RenderVolume(cmd, camera, screenW, screenH);

        // Draw domain wireframe box
        wi::primitive::AABB domainBox(gFluid.domainMin, gFluid.domainMax);
        wi::renderer::DrawBox(domainBox, XMFLOAT4(0.4f, 0.7f, 1.0f, 1.0f));

        // Draw smoke emitter sphere (warm orange wireframe)
        {
            XMFLOAT3 sc = gFluid.params.smokeSrcCenter;
            float sr = gFluid.params.smokeSrcRadius;
            wi::renderer::DrawSphere(
                wi::primitive::Sphere(sc, sr),
                XMFLOAT4(1.0f, 0.5f, 0.1f, 0.8f));
        }
    }

    // --- DUAL-SPHERE DRAGGING with left mouse button ---
    // Ray-cast against both obstacle and smoke emitter, pick nearest
    if (fluidReady)
    {
        auto pointer = wi::input::GetPointer();
        int mouseX = (int)pointer.x;
        int mouseY = (int)pointer.y;

        // Ray/sphere picking helper.
        //
        // Ray equation:
        //   P(t) = O + tD
        // where:
        //   O = ray origin
        //   D = ray direction
        //   t = distance along the ray
        //
        // Sphere equation:
        //   |P - C|^2 = r^2
        // where:
        //   C = sphere center
        //   r = sphere radius
        //
        // Substituting the ray into the sphere:
        //   |O + tD - C|^2 = r^2
        //
        // Define:
        //   oc = O - C
        // so oc is the vector from the sphere center to the ray origin.
        // Then:
        //   |oc + tD|^2 = r^2
        //
        // Expanding the squared length:
        //   (oc + tD) · (oc + tD) = r^2
        //   oc·oc + 2t(oc·D) + t^2(D·D) = r^2
        //
        // Here we assume D is normalized, so D·D = 1:
        //   t^2 + 2(oc·D)t + (oc·oc - r^2) = 0
        //
        // This is a quadratic equation in t. We store it in compact form as:
        //   b = oc·D
        //   c = oc·oc - r^2
        // giving:
        //   t^2 + 2bt + c = 0
        //
        // Geometric meaning:
        //   oc = vector from sphere center to ray origin
        //   b  = projection of oc onto the ray direction
        //   c  = squared distance from origin to center minus radius^2
        //      -> c > 0 : origin is outside the sphere
        //      -> c = 0 : origin lies on the sphere
        //      -> c < 0 : origin is inside the sphere
        //
        // Solving the quadratic:
        //   t = [-2b +- sqrt((2b)^2 - 4c)] / 2
        //     = -b +- sqrt(b^2 - c)
        //
        // We call:
        //   disc = b^2 - c
        // which is the discriminant.
        //
        // Discriminant meaning:
        //   disc < 0 : no real solution, the ray misses the sphere
        //   disc = 0 : one solution, the ray is tangent to the sphere
        //   disc > 0 : two solutions, the ray enters and exits the sphere
        //
        // The nearer hit is:
        //   tNear = -b - sqrt(disc)
        // and the farther hit is:
        //   tFar  = -b + sqrt(disc)
        //
        // For picking we only want the first intersection in front of the
        // camera, so we return tNear if it is positive. Negative t means the
        // hit lies behind the ray origin, so it is treated as "no hit" and
        // we return -1 in that case.
        auto raySphere = [](const wi::primitive::Ray& ray, const XMFLOAT3& center, float radius) -> float
        {
            float ocx = ray.origin.x - center.x;
            float ocy = ray.origin.y - center.y;
            float ocz = ray.origin.z - center.z;
            float b = ocx * ray.direction.x + ocy * ray.direction.y + ocz * ray.direction.z;
            float c = ocx * ocx + ocy * ocy + ocz * ocz - radius * radius;
            float disc = b * b - c;
            if (disc < 0.0f) return -1.0f;
            float t = -b - std::sqrt(disc);
            return (t > 0.0f) ? t : -1.0f;
        };

        // Start grab:
        // when the left mouse button is pressed and nothing is currently being
        // dragged, build a pick ray from the current mouse position through the
        // camera. Test that ray against both draggable spheres:
        //   1. the solid obstacle sphere
        //   2. the smoke emitter sphere
        //
        // Each raySphere() call returns the first positive hit distance t along
        // the ray, or -1 if that sphere is not hit in front of the camera.
        // If both are hit, we pick the smaller positive t, meaning the object
        // visually closest to the camera along that mouse ray. That becomes the
        // active drag target and we store its hit distance in dragDepth so the
        // object can continue moving on a plane of constant depth while the
        // mouse is dragged.
        if (wi::input::Press(wi::input::MOUSE_BUTTON_LEFT) && dragTarget == DragTarget::NONE)
        {
            wi::primitive::Ray ray = wi::renderer::GetPickRay(mouseX, mouseY,
                wi::scene::GetCamera().canvas, wi::scene::GetCamera());

            float tObs = raySphere(ray, gFluid.params.obsCenter, gFluid.params.obsRadius);
            float tSrc = raySphere(ray, gFluid.params.smokeSrcCenter, gFluid.params.smokeSrcRadius);

            float bestT = -1.0f;
            DragTarget bestTarget = DragTarget::NONE;

            if (tObs > 0.0f) { bestT = tObs; bestTarget = DragTarget::OBSTACLE; }
            if (tSrc > 0.0f && (bestT < 0.0f || tSrc < bestT)) { bestT = tSrc; bestTarget = DragTarget::SMOKE_SRC; }

            if (bestTarget != DragTarget::NONE)
            {
                dragTarget = bestTarget;
                dragDepth = bestT;
                lastDragCenter = (bestTarget == DragTarget::OBSTACLE) ?
                    gFluid.params.obsCenter : gFluid.params.smokeSrcCenter;
            }
        }

        // Drag update:
        // while the left mouse button stays down, rebuild the current pick ray
        // from the latest mouse position and evaluate the point that lies at
        // the previously stored hit distance:
        //   newCenter = ray.origin + dragDepth * ray.direction
        //
        // This does not keep the sphere locked to a world-space plane. Instead,
        // it keeps the dragged object at the same distance from the camera
        // along the current mouse ray, which feels like "grabbing" the sphere
        // and sliding it under the cursor without pushing it farther away or
        // pulling it closer.
        //
        // After that, clamp the center so the whole sphere remains inside the
        // fluid domain box. The radius is used in the clamp so the sphere
        // surface never crosses the domain boundaries.
        if (dragTarget != DragTarget::NONE && wi::input::Down(wi::input::MOUSE_BUTTON_LEFT))
        {
            wi::primitive::Ray ray = wi::renderer::GetPickRay(mouseX, mouseY,
                wi::scene::GetCamera().canvas, wi::scene::GetCamera());

            XMFLOAT3 newCenter;
            newCenter.x = ray.origin.x + dragDepth * ray.direction.x;
            newCenter.y = ray.origin.y + dragDepth * ray.direction.y;
            newCenter.z = ray.origin.z + dragDepth * ray.direction.z;

            float r = (dragTarget == DragTarget::OBSTACLE) ?
                gFluid.params.obsRadius : gFluid.params.smokeSrcRadius;

            // Clamp inside domain
            newCenter.x = std::max(gFluid.domainMin.x + r, std::min(gFluid.domainMax.x - r, newCenter.x));
            newCenter.y = std::max(gFluid.domainMin.y + r, std::min(gFluid.domainMax.y - r, newCenter.y));
            newCenter.z = std::max(gFluid.domainMin.z + r, std::min(gFluid.domainMax.z - r, newCenter.z));

            // Apply the drag result to whichever sphere is currently selected.
            //
            // Obstacle:
            // the obstacle is not just a visual marker; it actively couples to
            // the fluid simulation. For that reason we update both:
            //   1. its new center position
            //   2. an estimated linear velocity
            //
            // The velocity is approximated with a finite difference:
            //   vel = (currentPosition - previousPosition) / dt
            // using the center from the previous drag frame. This lets the
            // solver know how fast the obstacle is moving, so the moving sphere
            // can inject motion into the surrounding smoke.
            // See fluid_simulation.hlsl: in the obstacle setup pass, the
            // compute shader writes cb.obsVelX/Y/Z into the local velocity
            // field on obstacle cells and later uses those values again when
            // enforcing obstacle boundary conditions.
            //
            // Smoke source:
            // the emitter sphere only needs its position updated because it is
            // used as a source region for smoke injection and as a visual glow
            // region in the raymarcher. It does not need a drag velocity here.
            if (dragTarget == DragTarget::OBSTACLE)
            {
                XMFLOAT3 vel;
                vel.x = (newCenter.x - lastDragCenter.x) / std::max(dt, 1e-6f);
                vel.y = (newCenter.y - lastDragCenter.y) / std::max(dt, 1e-6f);
                vel.z = (newCenter.z - lastDragCenter.z) / std::max(dt, 1e-6f);
                gFluid.SetObstacle(newCenter, vel);
            }
            else
            {
                gFluid.SetSmokeSrc(newCenter);
            }

            // Store the new center as lastDragCenter for the next frame, so we can
            // estimate velocity from the center displacement in the next drag update.
            lastDragCenter = newCenter;
        }

        // Release:
        // when the left mouse button is no longer held, terminate the current
        // drag operation and clear the active target so a new pick can start
        // on a later click.
        //
        // For the obstacle we also explicitly reset its velocity to zero. While
        // dragging, the obstacle velocity is estimated every frame from the
        // center displacement and injected into the fluid solver so the moving
        // sphere can push the smoke. Once the drag stops, keeping the previous
        // velocity would make the obstacle appear to continue moving from the
        // simulation's point of view, so we send the same position again with a
        // zero velocity to mark it as stationary.
        // See fluid_simulation.hlsl:
        //   - in the obstacle setup pass, cb.obsVelX/Y/Z is copied into
        //     texU/texV/texW for obstacle cells
        //   - in the solid boundary pass, cb.obsVelX/Y/Z is reused as the
        //     obstacle boundary velocity when enforcing no-penetration
        //
        // Resetting the drag with obsVel = 0 ensures both passes treat the
        // sphere as stopped immediately after mouse release.
        if (dragTarget != DragTarget::NONE && !wi::input::Down(wi::input::MOUSE_BUTTON_LEFT))
        {
            if (dragTarget == DragTarget::OBSTACLE)
                gFluid.SetObstacle(gFluid.params.obsCenter, XMFLOAT3(0, 0, 0));
            dragTarget = DragTarget::NONE;
        }
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

void SampleRenderPath::Compose(wi::graphics::CommandList cmd) const
{
    wi::graphics::GraphicsDevice* device = wi::graphics::GetDevice();

    // Step 1: draw the already rendered 3D scene to the backbuffer.
    // RenderPath3D has already rasterized the normal scene content into its
    // own render target, including the ground plane and the debug helper
    // primitives issued earlier in Update() such as DrawBox() and DrawSphere().
    // Those helpers are regular debug line geometry, not volumetric rendering.
    // Compose() does not rebuild that scene; it simply copies the resulting
    // canvas-sized 2D image.
    wi::image::Params fx;
    fx.blendFlag = wi::enums::BLENDMODE_OPAQUE;
    fx.quality = wi::image::QUALITY_LINEAR;
    fx.enableFullScreen();
    wi::image::Draw(GetLastPostprocessRT(), fx, cmd);

    // Step 2: draw the fluid render target on top as a canvas-sized overlay.
    // gFluid.renderOutput was generated earlier by FluidSimulation3D::RenderVolume()
    // with the fluid_raymarching.hlsl compute shader. That pass ray-marches the
    // 3D simulation textures, shades special regions such as the obstacle and
    // smoke source spheres, and writes one 2D RGBA result per screen pixel.
    // The RGB written by that shader is accumulated in premultiplied form, so
    // here we composite the overlay with premultiplied alpha blending.
    //
    // This keeps the volumetric rendering path simple, but it also means the
    // fluid overlay is composited after the main 3D scene instead of being
    // depth-integrated into the scene rasterization path.
    if (gFluid.IsReady() && gFluid.renderOverlay && gFluid.renderOutput.IsValid())
    {
        wi::image::Params fluidFx;
        fluidFx.enableFullScreen();
        fluidFx.blendFlag = wi::enums::BLENDMODE_PREMULTIPLIED;
        wi::image::Draw(&gFluid.renderOutput, fluidFx, cmd);
    }

    // Step 3: draw the GUI last so it stays above both scene and fluid overlay.
    RenderPath2D::Compose(cmd);
}

