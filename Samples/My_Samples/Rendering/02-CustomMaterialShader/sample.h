#pragma once

#include "stdafx.h"

class SampleRenderPath : public wi::RenderPath3D
{
public:
	void Load() override;
	void Update(float dt) override;
	void ResizeLayout() override;
};

class SampleApp : public wi::Application
{
public:
	SampleRenderPath renderPath;

	void Initialize() override
	{
		wi::Application::Initialize();

		infoDisplay.active = true;
		infoDisplay.watermark = true;
		infoDisplay.fpsinfo = true;
		infoDisplay.resolution = true;
		infoDisplay.device_name = true;

		renderPath.init(canvas);
		renderPath.Load();
		ActivatePath(&renderPath);
	}
};
