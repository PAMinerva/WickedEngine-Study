#include "stdafx.h"
#include "sample.h"

#ifdef SDL2
#include <SDL2/SDL.h>

#define MAX_LOADSTRING 100

// Global Variables:
wi::Application application;					// Application class

int sdl_loop(wi::Application &application)
{
    SDL_Event event;

    bool quit = false;
    while (!quit)
    {
        SDL_PumpEvents();
        application.Run();

        while( SDL_PollEvent(&event)) 
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
                        application.SetWindow(application.window);
                        break;
                    default:
                        break;
                    }
                default:
                    break;
            }
        }
    }

    return 0;
}

int main(const int argc, const char* argv[])
{
    wi::Application application;
    #ifdef WickedEngine_SHADER_DIR
    wi::renderer::SetShaderSourcePath(WickedEngine_SHADER_DIR);
    #endif

    application.infoDisplay.active = true;
    application.infoDisplay.watermark = true;
    application.infoDisplay.resolution = true;
    application.infoDisplay.fpsinfo = true;
    application.infoDisplay.colorspace = true;
    application.infoDisplay.device_name = true;

    sdl2::sdlsystem_ptr_t system = sdl2::make_sdlsystem(SDL_INIT_EVERYTHING | SDL_INIT_EVENTS);
    sdl2::window_ptr_t window = sdl2::make_window(
            "Hello Triangle",
            SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
            1920, 1080,
            SDL_WINDOW_SHOWN | SDL_WINDOW_VULKAN | SDL_WINDOW_ALLOW_HIGHDPI);

    SDL_Event event;

    application.SetWindow(window.get());

    Sample sample;  // Sample class

    struct DPI_INFO
    {
        float ddpi = 0;
        float hdpi = 0;
        float vdpi = 0;
    } dpi_info;

    SDL_GetDisplayDPI(0, &dpi_info.ddpi, &dpi_info.hdpi, &dpi_info.vdpi);
    application.canvas.dpi = dpi_info.hdpi;

    application.Initialize(); // Perform graphics initialization for the application (create device, swapchain, etc.)
    application.ActivatePath(sample.GetRenderPath3D());

    int ret = sdl_loop(application);

    SDL_Quit();

    return ret;
}
#endif