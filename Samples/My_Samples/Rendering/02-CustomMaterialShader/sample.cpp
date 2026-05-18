#include "sample.h"

#include <wiECS.h>
#include <wiScene_Components.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>

namespace
{
	wi::ecs::Entity customCube = wi::ecs::INVALID_ENTITY;
	wi::ecs::Entity customMaterial = wi::ecs::INVALID_ENTITY;
	wi::ecs::Entity referenceCube = wi::ecs::INVALID_ENTITY;

	wi::graphics::Shader customPixelShader;
	int customShaderID = -1;

	wi::gui::Label titleLabel;
	wi::gui::Label statusLabel;

	float elapsedTime = 0.0f;

	uint32_t FloatBits(float value)
	{
		uint32_t bits = 0;
		std::memcpy(&bits, &value, sizeof(bits));
		return bits;
	}

	bool IsRegisteredCustomShaderAlive()
	{
		const wi::vector<wi::renderer::CustomShader>& shaders = wi::renderer::GetCustomShaders();
		return customShaderID >= 0 &&
			customShaderID < static_cast<int>(shaders.size()) &&
			shaders[customShaderID].name == "SampleCustomMaterial";
	}

	bool RegisterCustomMaterialShader()
	{
		if (IsRegisteredCustomShaderAlive())
		{
			return true;
		}
		if (!wi::initializer::IsInitializeFinished() || wi::renderer::IsPipelineCreationActive() > 0)
		{
			return false;
		}

		using namespace wi::graphics;
		using namespace wi::renderer;
		using namespace wi::enums;

		const std::string originalBinPath = GetShaderPath();
		const std::string originalSrcPath = GetShaderSourcePath();
		const std::string shaderBinPath = wi::helper::GetCurrentPath() + "/shaders/";

		SetShaderPath(shaderBinPath);
#ifdef SHADER_SOURCE_DIR
		SetShaderSourcePath(SHADER_SOURCE_DIR);
#else
		SetShaderSourcePath(shaderBinPath);
#endif

		const bool loaded = LoadShader(ShaderStage::PS, customPixelShader, "custom_materialPS.cso");

		SetShaderPath(originalBinPath);
		SetShaderSourcePath(originalSrcPath);

		if (!loaded)
		{
			statusLabel.SetText("custom_materialPS.hlsl failed to compile");
			return false;
		}

		GraphicsDevice* device = wi::graphics::GetDevice();
		CustomShader customShader;
		customShader.name = "SampleCustomMaterial";
		customShader.filterMask = wi::enums::FILTER_OPAQUE;

		PipelineStateDesc desc = {};
		desc.vs = GetShader(VSTYPE_OBJECT_PREPASS);
		desc.ps = GetShader(PSTYPE_OBJECT_PREPASS);
		desc.rs = GetRasterizerState(RSTYPE_FRONT);
		desc.dss = GetDepthStencilState(DSSTYPE_DEFAULT);
		desc.pt = PrimitiveTopology::TRIANGLELIST;
		if (!device->CreatePipelineState(&desc, &customShader.pso[RENDERPASS_PREPASS]))
		{
			statusLabel.SetText("custom prepass PSO creation failed");
			return false;
		}

		desc = {};
		desc.vs = GetShader(VSTYPE_OBJECT_COMMON);
		desc.ps = &customPixelShader;
		desc.rs = GetRasterizerState(RSTYPE_FRONT);
		desc.bs = GetBlendState(BSTYPE_OPAQUE);
		desc.dss = GetDepthStencilState(DSSTYPE_DEPTHREADEQUAL);
		desc.pt = PrimitiveTopology::TRIANGLELIST;
		if (!device->CreatePipelineState(&desc, &customShader.pso[RENDERPASS_MAIN]))
		{
			statusLabel.SetText("custom main PSO creation failed");
			return false;
		}

		desc = {};
		desc.vs = GetShader(VSTYPE_SHADOW);
		desc.rs = GetRasterizerState(RSTYPE_SHADOW);
		desc.bs = GetBlendState(BSTYPE_OPAQUE);
		desc.dss = GetDepthStencilState(DSSTYPE_SHADOW);
		desc.pt = PrimitiveTopology::TRIANGLELIST;
		if (!device->CreatePipelineState(&desc, &customShader.pso[RENDERPASS_SHADOW]))
		{
			statusLabel.SetText("custom shadow PSO creation failed");
			return false;
		}

		customShaderID = RegisterCustomShader(customShader);

		if (wi::scene::MaterialComponent* material = wi::scene::GetScene().materials.GetComponent(customMaterial))
		{
			material->SetCustomShaderID(customShaderID);
			material->SetDirty();
		}

		statusLabel.SetText("custom HLSL material shader registered");
		return true;
	}

	void UpdateCustomMaterialParameters()
	{
		wi::scene::MaterialComponent* material = wi::scene::GetScene().materials.GetComponent(customMaterial);
		if (material == nullptr)
		{
			return;
		}

		const float stripeCount = 4.0f;
		const float warpAmount = 0.04f;
		const float scrollSpeed = 0.15f;

		material->userdata = uint4(
			FloatBits(stripeCount),
			FloatBits(warpAmount),
			FloatBits(scrollSpeed),
			0);
		material->SetDirty();
	}

	void InitScene()
	{
		using namespace wi::scene;

		Scene& scene = GetScene();

		WeatherComponent& weather = scene.weathers.Create(wi::ecs::CreateEntity());
		weather.ambient = XMFLOAT3(0.18f, 0.20f, 0.24f);
		weather.horizon = XMFLOAT3(0.12f, 0.16f, 0.20f);
		weather.zenith = XMFLOAT3(0.02f, 0.03f, 0.05f);

		const wi::ecs::Entity lightEntity = scene.Entity_CreateLight("KeyLight", XMFLOAT3(0, 4, -2));
		LightComponent& light = *scene.lights.GetComponent(lightEntity);
		light.SetType(LightComponent::DIRECTIONAL);
		light.intensity = 3.0f;
		light.SetCastShadow(true);
		TransformComponent* lightTransform = scene.transforms.GetComponent(lightEntity);
		lightTransform->RotateRollPitchYaw(XMFLOAT3(
			wi::math::DegreesToRadians(55.0f),
			wi::math::DegreesToRadians(-35.0f),
			0.0f));
		lightTransform->UpdateTransform();

		customCube = scene.Entity_CreateCube("custom_hlsl_material_cube");
		customMaterial = customCube;
		TransformComponent* customTransform = scene.transforms.GetComponent(customCube);
		customTransform->Translate(XMFLOAT3(-1.35f, 0.0f, 4.0f));
		customTransform->RotateRollPitchYaw(XMFLOAT3(0.0f, wi::math::DegreesToRadians(25.0f), 0.0f));
		customTransform->Scale(XMFLOAT3(0.9f, 0.9f, 0.9f));
		customTransform->UpdateTransform();

		MaterialComponent* material = scene.materials.GetComponent(customMaterial);
		material->shaderType = MaterialComponent::SHADERTYPE_UNLIT;
		material->SetBaseColor(XMFLOAT4(0.95f, 0.22f, 0.09f, 1.0f));
		material->SetRoughness(0.75f);
		UpdateCustomMaterialParameters();

		referenceCube = scene.Entity_CreateCube("regular_wicked_material_cube");
		TransformComponent* referenceTransform = scene.transforms.GetComponent(referenceCube);
		referenceTransform->Translate(XMFLOAT3(1.35f, 0.0f, 4.0f));
		referenceTransform->RotateRollPitchYaw(XMFLOAT3(0.0f, wi::math::DegreesToRadians(-25.0f), 0.0f));
		referenceTransform->Scale(XMFLOAT3(0.9f, 0.9f, 0.9f));
		referenceTransform->UpdateTransform();

		MaterialComponent* referenceMaterial = scene.materials.GetComponent(referenceCube);
		referenceMaterial->shaderType = MaterialComponent::SHADERTYPE_UNLIT;
		referenceMaterial->SetBaseColor(XMFLOAT4(0.18f, 0.34f, 0.92f, 1.0f));

		CameraComponent& camera = GetCamera();
		camera.CreatePerspective(1280.0f, 720.0f, 0.1f, 100.0f, XM_PIDIV4);

		TransformComponent cameraTransform;
		cameraTransform.ClearTransform();
		cameraTransform.Translate(XMFLOAT3(0.0f, 1.0f, -2.8f));
		cameraTransform.RotateRollPitchYaw(XMFLOAT3(wi::math::DegreesToRadians(8.0f), 0.0f, 0.0f));
		cameraTransform.UpdateTransform();
		camera.TransformCamera(cameraTransform);
	}
}

void SampleRenderPath::Load()
{
	titleLabel.Create("TitleLabel");
	titleLabel.SetText("Custom HLSL Material Shader");
	titleLabel.SetSize(XMFLOAT2(420, 24));
	titleLabel.font.params.h_align = wi::font::WIFALIGN_CENTER;
	GetGUI().AddWidget(&titleLabel);

	statusLabel.Create("StatusLabel");
	statusLabel.SetText("waiting for renderer shaders...");
	statusLabel.SetSize(XMFLOAT2(420, 42));
	statusLabel.font.params.h_align = wi::font::WIFALIGN_CENTER;
	GetGUI().AddWidget(&statusLabel);

	InitScene();
	RenderPath3D::Load();
}

void SampleRenderPath::ResizeLayout()
{
	RenderPath3D::ResizeLayout();

	const float screenW = GetLogicalWidth();
	titleLabel.SetPos(XMFLOAT2(screenW * 0.5f - titleLabel.scale.x * 0.5f, 18.0f));
	statusLabel.SetPos(XMFLOAT2(screenW * 0.5f - statusLabel.scale.x * 0.5f, 46.0f));
}

void SampleRenderPath::Update(float dt)
{
	elapsedTime += dt;

	RegisterCustomMaterialShader();
	UpdateCustomMaterialParameters();

	wi::scene::Scene& scene = wi::scene::GetScene();
	if (wi::scene::TransformComponent* transform = scene.transforms.GetComponent(customCube))
	{
		transform->RotateRollPitchYaw(XMFLOAT3(0.0f, dt * 0.55f, 0.0f));
	}
	if (wi::scene::TransformComponent* transform = scene.transforms.GetComponent(referenceCube))
	{
		transform->RotateRollPitchYaw(XMFLOAT3(0.0f, -dt * 0.35f, 0.0f));
	}

	RenderPath3D::Update(dt);
}
