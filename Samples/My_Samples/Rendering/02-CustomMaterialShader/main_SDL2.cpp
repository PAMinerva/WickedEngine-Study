#include "sample.h"
#include "stdafx.h"
#include "sdl2.h"

SampleApp app;

int main(int argc, char* argv[])
{
	wi::arguments::Parse(argc, argv);

	sdl2::sdlsystem_ptr_t system = sdl2::make_sdlsystem(SDL_INIT_EVERYTHING | SDL_INIT_EVENTS);
	if (*system)
	{
		wilog_error("Error creating SDL2 system");
	}

	sdl2::window_ptr_t window = sdl2::make_window(
		"02 - Custom HLSL Material Shader",
		SDL_WINDOWPOS_CENTERED,
		SDL_WINDOWPOS_CENTERED,
		1280,
		800,
		SDL_WINDOW_SHOWN | SDL_WINDOW_VULKAN | SDL_WINDOW_RESIZABLE);
	if (!window)
	{
		wilog_error("Error creating window");
	}

	app.SetWindow(window.get());

	bool quit = false;
	while (!quit)
	{
		app.Run();

		SDL_Event event;
		while (SDL_PollEvent(&event))
		{
			switch (event.type)
			{
			case SDL_QUIT:
				quit = true;
				break;
			case SDL_WINDOWEVENT:
				switch (event.window.event)
				{
				case SDL_WINDOWEVENT_CLOSE:
					quit = true;
					break;
				case SDL_WINDOWEVENT_RESIZED:
					app.SetWindow(app.window);
					break;
				case SDL_WINDOWEVENT_FOCUS_LOST:
					app.is_window_active = false;
					break;
				case SDL_WINDOWEVENT_FOCUS_GAINED:
					app.is_window_active = true;
					if (wi::shadercompiler::GetRegisteredShaderCount() > 0 && !wi::renderer::IsPipelineCreationActive())
					{
						std::thread([] {
							if (wi::shadercompiler::CheckRegisteredShadersOutdated())
							{
								wi::eventhandler::Subscribe_Once(wi::eventhandler::EVENT_THREAD_SAFE_POINT, [](uint64_t) {
									wi::renderer::ReloadShaders();
								});
							}
						}).detach();
					}
					break;
				default:
					break;
				}
				break;
			default:
				break;
			}

			wi::input::sdlinput::ProcessEvent(event);
		}
	}

	SDL_Quit();
	return 0;
}
