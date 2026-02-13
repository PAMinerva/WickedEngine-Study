#pragma once

#include <wiRenderPath3D.h>
#include <wiApplication.h>

class SampleRenderPath : public wi::RenderPath3D
{
public:
	void ResizeLayout() override;

	void Load() override;
	void Update(float dt) override;
	void FixedUpdate() override;
};

class SampleApp : public wi::Application
{
public:
    SampleRenderPath renderer;
	void Initialize() override;

    // const char *title {"01 - Hello Triangle"};
    // const char* GetTitle() {return title;};
};
