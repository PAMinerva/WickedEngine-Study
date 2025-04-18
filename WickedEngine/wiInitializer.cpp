#include "wiInitializer.h"
#include "WickedEngine.h"

#include <string>
#include <thread>
#include <atomic>

namespace wi::initializer
{
	static std::atomic_bool initializationStarted{ false };
	static wi::jobsystem::context ctx;
	static wi::Timer timer;
	static std::atomic_bool systems[INITIALIZED_SYSTEM_COUNT]{};

	void InitializeComponentsImmediate()
	{
		if (IsInitializeFinished())
			return;
		if (!initializationStarted.load())
		{
			InitializeComponentsAsync();
		}
		WaitForInitializationsToFinish();
	}
	void InitializeComponentsAsync()
	{
		if (IsInitializeFinished())
			return;
		timer.record();

		initializationStarted.store(true);

#if defined(PLATFORM_WINDOWS_DESKTOP)
		static constexpr const char* platform_string = "Windows";
#elif defined(PLATFORM_LINUX)
		static constexpr const char* platform_string = "Linux";
#elif defined(PLATFORM_PS5)
		static constexpr const char* platform_string = "PS5";
#elif defined(PLATFORM_XBOX)
		static constexpr const char* platform_string = "Xbox";
#endif // PLATFORM

		wilog("\n[wi::initializer] Initializing Wicked Engine, please wait...\nVersion: %s\nPlatform: %s", wi::version::GetVersionString(), platform_string)
		
		size_t shaderdump_count = wi::renderer::GetShaderDumpCount();
		if (shaderdump_count > 0)
		{
			wi::backlog::post("\nEmbedded shaders found: " + std::to_string(shaderdump_count));
		}
		else
		{
			wi::backlog::post("\nNo embedded shaders found, shaders will be compiled at runtime if needed.\n\tShader source path: " + wi::renderer::GetShaderSourcePath() + "\n\tShader binary path: " + wi::renderer::GetShaderPath());
		}

		wi::backlog::post("");
		wi::jobsystem::Initialize();

		wi::backlog::post("");
		wi::jobsystem::Execute(ctx, [](wi::jobsystem::JobArgs args) { wi::font::Initialize(); systems[INITIALIZED_SYSTEM_FONT].store(true); });
		wi::jobsystem::Execute(ctx, [](wi::jobsystem::JobArgs args) { wi::image::Initialize(); systems[INITIALIZED_SYSTEM_IMAGE].store(true); });
		wi::jobsystem::Execute(ctx, [](wi::jobsystem::JobArgs args) { wi::input::Initialize(); systems[INITIALIZED_SYSTEM_INPUT].store(true); });
		wi::jobsystem::Execute(ctx, [](wi::jobsystem::JobArgs args) { wi::renderer::Initialize(); systems[INITIALIZED_SYSTEM_RENDERER].store(true); });
		wi::jobsystem::Execute(ctx, [](wi::jobsystem::JobArgs args) { wi::texturehelper::Initialize(); systems[INITIALIZED_SYSTEM_TEXTUREHELPER].store(true); });
		wi::jobsystem::Execute(ctx, [](wi::jobsystem::JobArgs args) { wi::HairParticleSystem::Initialize(); systems[INITIALIZED_SYSTEM_HAIRPARTICLESYSTEM].store(true); });
		wi::jobsystem::Execute(ctx, [](wi::jobsystem::JobArgs args) { wi::EmittedParticleSystem::Initialize(); systems[INITIALIZED_SYSTEM_EMITTEDPARTICLESYSTEM].store(true); });
		wi::jobsystem::Execute(ctx, [](wi::jobsystem::JobArgs args) { wi::Ocean::Initialize(); systems[INITIALIZED_SYSTEM_OCEAN].store(true); });
		wi::jobsystem::Execute(ctx, [](wi::jobsystem::JobArgs args) { wi::gpusortlib::Initialize(); systems[INITIALIZED_SYSTEM_GPUSORTLIB].store(true); });
		wi::jobsystem::Execute(ctx, [](wi::jobsystem::JobArgs args) { wi::GPUBVH::Initialize(); systems[INITIALIZED_SYSTEM_GPUBVH].store(true); });
		wi::jobsystem::Execute(ctx, [](wi::jobsystem::JobArgs args) { wi::physics::Initialize(); systems[INITIALIZED_SYSTEM_PHYSICS].store(true); });
		wi::jobsystem::Execute(ctx, [](wi::jobsystem::JobArgs args) { wi::TrailRenderer::Initialize(); systems[INITIALIZED_SYSTEM_TRAILRENDERER].store(true); });

		// Initialize these immediately:
		wi::lua::Initialize(); systems[INITIALIZED_SYSTEM_LUA].store(true);
		wi::audio::Initialize(); systems[INITIALIZED_SYSTEM_AUDIO].store(true);

		std::thread([] {
			wi::jobsystem::Wait(ctx);
			wilog("\n[wi::initializer] Wicked Engine Initialized (%d ms)", (int)std::round(timer.elapsed()));
		}).detach();

	}

	bool IsInitializeFinished(INITIALIZED_SYSTEM system)
	{
		if (system == INITIALIZED_SYSTEM_COUNT)
		{
			return initializationStarted.load() && !wi::jobsystem::IsBusy(ctx);
		}
		else
		{
			return systems[system].load();
		}
	}

	void WaitForInitializationsToFinish()
	{
		wi::jobsystem::Wait(ctx);
	}
}
