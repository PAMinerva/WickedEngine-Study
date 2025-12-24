#pragma once

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
};
