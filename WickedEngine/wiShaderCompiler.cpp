#include "wiShaderCompiler.h"
#include "wiBacklog.h"
#include "wiPlatform.h"
#include "wiHelper.h"
#include "wiArchive.h"
#include "wiUnorderedSet.h"

#include <mutex>

#ifdef PLATFORM_WINDOWS_DESKTOP
#define SHADERCOMPILER_ENABLED
#define SHADERCOMPILER_ENABLED_DXCOMPILER
#define SHADERCOMPILER_ENABLED_D3DCOMPILER
#include <wrl/client.h>
using namespace Microsoft::WRL;
#endif // _WIN32

#ifdef PLATFORM_LINUX
#define SHADERCOMPILER_ENABLED
#define SHADERCOMPILER_ENABLED_DXCOMPILER
#define __RPC_FAR
#define ComPtr CComPtr
#include "Utility/dxc/WinAdapter.h"
#endif // PLATFORM_LINUX

#ifdef SHADERCOMPILER_ENABLED_DXCOMPILER
#include "Utility/dxc/dxcapi.h"
#endif // SHADERCOMPILER_ENABLED_DXCOMPILER

#ifdef SHADERCOMPILER_ENABLED_D3DCOMPILER
#include <d3dcompiler.h>
#endif // SHADERCOMPILER_ENABLED_D3DCOMPILER

#if __has_include("wiShaderCompiler_XBOX.h")
#include "wiShaderCompiler_XBOX.h"
#define SHADERCOMPILER_XBOX_INCLUDED
#endif // __has_include("wiShaderCompiler_XBOX.h")

#if __has_include("wiShaderCompiler_PS5.h") && !defined(PLATFORM_PS5)
#include "wiShaderCompiler_PS5.h"
#define SHADERCOMPILER_PS5_INCLUDED
#endif // __has_include("wiShaderCompiler_PS5.h") && !PLATFORM_PS5

using namespace wi::graphics;

namespace wi::shadercompiler
{

#ifdef SHADERCOMPILER_ENABLED_DXCOMPILER
	struct InternalState_DXC
	{
		// To make it more convenient for you to use GetProcAddress to call DxcCreateInstance,
		// the DxcCreateInstanceProc typedef is provided by the dxcapi.h header file.
		// While the DxcCreateInstance function is similar to CoCreateInstance, there's no COM involvement.
		DxcCreateInstanceProc DxcCreateInstance = nullptr;

		InternalState_DXC(const std::string& modifier = "")
		{
#ifdef _WIN32
			// Load the dxcompiler.dll library
			const std::string library = "dxcompiler" + modifier + ".dll";
			HMODULE dxcompiler = wiLoadLibrary(library.c_str());
#elif defined(PLATFORM_LINUX)
			const std::string library = "./libdxcompiler" + modifier + ".so";
			HMODULE dxcompiler = wiLoadLibrary(library.c_str());
#endif
			if (dxcompiler != nullptr)
			{
				// Load the DxcCreateInstance function from the library
				DxcCreateInstance = (DxcCreateInstanceProc)wiGetProcAddress(dxcompiler, "DxcCreateInstance");
				if (DxcCreateInstance != nullptr)
				{
					// Retrieve the DirectX Shader Compiler Interface using the DxcCreateInstance function with the CLSID_DxcCompiler CLSID
					// as suggested in the documentation of the IDxcCompiler3 interface
					ComPtr<IDxcCompiler3> dxcCompiler;
					HRESULT hr = DxcCreateInstance(CLSID_DxcCompiler, IID_PPV_ARGS(&dxcCompiler));
					assert(SUCCEEDED(hr));
					// Get compiler version by querying the IDxcVersionInfo interface from the IDxcCompiler3 interface
					ComPtr<IDxcVersionInfo> info;
					hr = dxcCompiler->QueryInterface(IID_PPV_ARGS(&info));
					assert(SUCCEEDED(hr));
					uint32_t minor = 0;
					uint32_t major = 0;
					hr = info->GetVersion(&major, &minor);
					assert(SUCCEEDED(hr));
					wi::backlog::post("wi::shadercompiler: loaded " + library + " (version: " + std::to_string(major) + "." + std::to_string(minor) + ")");
				}
			}
			else
			{
				wi::backlog::post("wi::shadercompiler: could not load library " + library, wi::backlog::LogLevel::Error);
#ifdef PLATFORM_LINUX
				wi::backlog::post(dlerror(), wi::backlog::LogLevel::Error); // print dlopen() error detail: https://linux.die.net/man/3/dlerror
#endif // PLATFORM_LINUX
			}

		}
	};
	inline InternalState_DXC& dxc_compiler()
	{
		static InternalState_DXC internal_state;
		return internal_state;
	}
	inline InternalState_DXC& dxc_compiler_xs()
	{
		static InternalState_DXC internal_state("_xs");
		return internal_state;
	}

	void Compile_DXCompiler(const CompilerInput& input, CompilerOutput& output)
	{
		InternalState_DXC& compiler_internal = input.format == ShaderFormat::HLSL6_XS ? dxc_compiler_xs() : dxc_compiler();
		if (compiler_internal.DxcCreateInstance == nullptr)
		{
			return;
		}

		ComPtr<IDxcUtils> dxcUtils;
		ComPtr<IDxcCompiler3> dxcCompiler;

		HRESULT hr = compiler_internal.DxcCreateInstance(CLSID_DxcUtils, IID_PPV_ARGS(&dxcUtils));
		assert(SUCCEEDED(hr));
		hr = compiler_internal.DxcCreateInstance(CLSID_DxcCompiler, IID_PPV_ARGS(&dxcCompiler));
		assert(SUCCEEDED(hr));

		if (dxcCompiler == nullptr)
		{
			return;
		}

		wi::vector<uint8_t> shadersourcedata;
		if (!wi::helper::FileRead(input.shadersourcefilename, shadersourcedata))
		{
			return;
		}

		// https://github.com/microsoft/DirectXShaderCompiler/wiki/Using-dxc.exe-and-dxcompiler.dll#dxcompiler-dll-interface

		wi::vector<std::wstring> args = {
			//L"-res-may-alias", // Assume that unordered access views (UAVs) and shader resource views (SRVs) may alias
			//L"-flegacy-macro-expansion",
			//L"-no-legacy-cbuf-layout",
			//L"-pack-optimized", // this has problem with tessellation shaders: https://github.com/microsoft/DirectXShaderCompiler/issues/3362
			//L"-all-resources-bound",
			//L"-Gis", // Force IEEE strictness
			//L"-Gec", // Enable backward compatibility mode
			//L"-Ges", // Enable strict mode
			//L"-O0", // Optimization Level 0
			//L"-enable-16bit-types",
			L"-Wno-conversion", // Suppress warnings about implicit conversions that may change the value
		};

		if (has_flag(input.flags, Flags::DISABLE_OPTIMIZATION))
		{
			// Optimization Level 0: Disables optimizations
			// Additional options to include debug info in compiled shaders
			// Required to debug shaders using PIX or RenderDoc
			/*args.push_back(L"-Od");
			args.push_back(L"-Zi");
			args.push_back(L"-Qembed_debug");*/

			if (input.format == ShaderFormat::HLSL6 || input.format == ShaderFormat::HLSL6_XS || input.format == ShaderFormat::HLSL5)
			{
				// Optimization Level 0: Disables optimizations
				// Additional options to include debug info in compiled shaders
				// Required to debug shaders using PIX or RenderDoc
				args.push_back(L"-Od");
				args.push_back(L"-Zi");
				args.push_back(L"-Qembed_debug");
			}
			else if (input.format == ShaderFormat::SPIRV)
			{
				// generate NonSemantic.Shader.DebugInfo.100 extended instructions
				// to include debug info in compiled shaders and support source-level
				// shader debugging with tools such as RenderDoc
				args.push_back(L"-fspv-debug=vulkan-with-source");
			}
		}

		switch (input.format)
		{
		case ShaderFormat::HLSL6:
		case ShaderFormat::HLSL6_XS:
			// rootsig-define is used to define a default root signature for shaders that have one defined in the shader source; see globals.hlsli
			args.push_back(L"-rootsig-define"); args.push_back(L"WICKED_ENGINE_DEFAULT_ROOTSIGNATURE");
			if (has_flag(input.flags, Flags::STRIP_REFLECTION))
			{
				// By default, reflection data is stored in the shader. The -Qstrip_reflect flag can be used to remove it.
				args.push_back(L"-Qstrip_reflect"); // only valid in HLSL6 compiler
			}
			break;
		case ShaderFormat::SPIRV: // https://github.com/microsoft/DirectXShaderCompiler/blob/main/docs/SPIR-V.rst
			args.push_back(L"-spirv");
			args.push_back(L"-fspv-target-env=vulkan1.3");
			args.push_back(L"-fvk-use-dx-layout");
			args.push_back(L"-fvk-use-dx-position-w");
			//args.push_back(L"-fvk-b-shift"); args.push_back(L"0"); args.push_back(L"0"); // this is implicit
			args.push_back(L"-fvk-t-shift"); args.push_back(L"1000"); args.push_back(L"0");
			args.push_back(L"-fvk-u-shift"); args.push_back(L"2000"); args.push_back(L"0");
			args.push_back(L"-fvk-s-shift"); args.push_back(L"3000"); args.push_back(L"0");
			break;
		default:
			assert(0);
			return;
		}

		ShaderModel minshadermodel = input.minshadermodel;

		// global minshadermodel override:
		//minshadermodel = std::max(minshadermodel, ShaderModel::SM_6_2);

		if (input.format == ShaderFormat::HLSL6_XS)
		{
			minshadermodel = ShaderModel::SM_6_6;
		}

		args.push_back(L"-T"); // Set target profile
		switch (input.stage)
		{
		case ShaderStage::MS:
			switch (minshadermodel)
			{
			default:
				args.push_back(L"ms_6_5");
				break;
			case ShaderModel::SM_6_6:
				args.push_back(L"ms_6_6");
				break;
			case ShaderModel::SM_6_7:
				args.push_back(L"ms_6_7");
				break;
			}
			break;
		case ShaderStage::AS:
			switch (minshadermodel)
			{
			default:
				args.push_back(L"as_6_5");
				break;
			case ShaderModel::SM_6_6:
				args.push_back(L"as_6_6");
				break;
			case ShaderModel::SM_6_7:
				args.push_back(L"as_6_7");
				break;
			}
			break;
		case ShaderStage::VS:
			switch (minshadermodel)
			{
			default:
				args.push_back(L"vs_6_0");
				break;
			case ShaderModel::SM_6_1:
				args.push_back(L"vs_6_1");
				break;
			case ShaderModel::SM_6_2:
				args.push_back(L"vs_6_2");
				break;
			case ShaderModel::SM_6_3:
				args.push_back(L"vs_6_3");
				break;
			case ShaderModel::SM_6_4:
				args.push_back(L"vs_6_4");
				break;
			case ShaderModel::SM_6_5:
				args.push_back(L"vs_6_5");
				break;
			case ShaderModel::SM_6_6:
				args.push_back(L"vs_6_6");
				break;
			case ShaderModel::SM_6_7:
				args.push_back(L"vs_6_7");
				break;
			}
			break;
		case ShaderStage::HS:
			switch (minshadermodel)
			{
			default:
				args.push_back(L"hs_6_0");
				break;
			case ShaderModel::SM_6_1:
				args.push_back(L"hs_6_1");
				break;
			case ShaderModel::SM_6_2:
				args.push_back(L"hs_6_2");
				break;
			case ShaderModel::SM_6_3:
				args.push_back(L"hs_6_3");
				break;
			case ShaderModel::SM_6_4:
				args.push_back(L"hs_6_4");
				break;
			case ShaderModel::SM_6_5:
				args.push_back(L"hs_6_5");
				break;
			case ShaderModel::SM_6_6:
				args.push_back(L"hs_6_6");
				break;
			case ShaderModel::SM_6_7:
				args.push_back(L"hs_6_7");
				break;
			}
			break;
		case ShaderStage::DS:
			switch (minshadermodel)
			{
			default:
				args.push_back(L"ds_6_0");
				break;
			case ShaderModel::SM_6_1:
				args.push_back(L"ds_6_1");
				break;
			case ShaderModel::SM_6_2:
				args.push_back(L"ds_6_2");
				break;
			case ShaderModel::SM_6_3:
				args.push_back(L"ds_6_3");
				break;
			case ShaderModel::SM_6_4:
				args.push_back(L"ds_6_4");
				break;
			case ShaderModel::SM_6_5:
				args.push_back(L"ds_6_5");
				break;
			case ShaderModel::SM_6_6:
				args.push_back(L"ds_6_6");
				break;
			case ShaderModel::SM_6_7:
				args.push_back(L"ds_6_7");
				break;
			}
			break;
		case ShaderStage::GS:
			switch (minshadermodel)
			{
			default:
				args.push_back(L"gs_6_0");
				break;
			case ShaderModel::SM_6_1:
				args.push_back(L"gs_6_1");
				break;
			case ShaderModel::SM_6_2:
				args.push_back(L"gs_6_2");
				break;
			case ShaderModel::SM_6_3:
				args.push_back(L"gs_6_3");
				break;
			case ShaderModel::SM_6_4:
				args.push_back(L"gs_6_4");
				break;
			case ShaderModel::SM_6_5:
				args.push_back(L"gs_6_5");
				break;
			case ShaderModel::SM_6_6:
				args.push_back(L"gs_6_6");
				break;
			case ShaderModel::SM_6_7:
				args.push_back(L"gs_6_7");
				break;
			}
			break;
		case ShaderStage::PS:
			switch (minshadermodel)
			{
			default:
				args.push_back(L"ps_6_0");
				break;
			case ShaderModel::SM_6_1:
				args.push_back(L"ps_6_1");
				break;
			case ShaderModel::SM_6_2:
				args.push_back(L"ps_6_2");
				break;
			case ShaderModel::SM_6_3:
				args.push_back(L"ps_6_3");
				break;
			case ShaderModel::SM_6_4:
				args.push_back(L"ps_6_4");
				break;
			case ShaderModel::SM_6_5:
				args.push_back(L"ps_6_5");
				break;
			case ShaderModel::SM_6_6:
				args.push_back(L"ps_6_6");
				break;
			case ShaderModel::SM_6_7:
				args.push_back(L"ps_6_7");
				break;
			}
			break;
		case ShaderStage::CS:
			switch (minshadermodel)
			{
			default:
				args.push_back(L"cs_6_0");
				break;
			case ShaderModel::SM_6_1:
				args.push_back(L"cs_6_1");
				break;
			case ShaderModel::SM_6_2:
				args.push_back(L"cs_6_2");
				break;
			case ShaderModel::SM_6_3:
				args.push_back(L"cs_6_3");
				break;
			case ShaderModel::SM_6_4:
				args.push_back(L"cs_6_4");
				break;
			case ShaderModel::SM_6_5:
				args.push_back(L"cs_6_5");
				break;
			case ShaderModel::SM_6_6:
				args.push_back(L"cs_6_6");
				break;
			case ShaderModel::SM_6_7:
				args.push_back(L"cs_6_7");
				break;
			}
			break;
		case ShaderStage::LIB:
			switch (minshadermodel)
			{
			default:
				args.push_back(L"lib_6_5");
				break;
			case ShaderModel::SM_6_6:
				args.push_back(L"lib_6_6");
				break;
			case ShaderModel::SM_6_7:
				args.push_back(L"lib_6_7");
				break;
			}
			break;
		default:
			assert(0);
			return;
		}

		for (auto& x : input.defines) // Define macros
		{
			args.push_back(L"-D");
			wi::helper::StringConvert(x, args.emplace_back());
		}

		for (auto& x : input.include_directories) // Include directories
		{
			args.push_back(L"-I");
			wi::helper::StringConvert(x, args.emplace_back());
		}

#ifdef SHADERCOMPILER_XBOX_INCLUDED
		if (input.format == ShaderFormat::HLSL6_XS)
		{
			wi::shadercompiler::xbox::AddArguments(input, args);
		}
#endif // SHADERCOMPILER_XBOX_INCLUDED

		// Entry point parameter:
		std::wstring wentry;
		wi::helper::StringConvert(input.entrypoint, wentry);
		args.push_back(L"-E");
		args.push_back(wentry.c_str());

		// Add source file name as last parameter. This will be displayed in error messages
		std::wstring wsource;
		wi::helper::StringConvert(wi::helper::GetFileNameFromPath(input.shadersourcefilename), wsource);
		args.push_back(wsource.c_str());

		DxcBuffer Source;
		Source.Ptr = shadersourcedata.data();
		Source.Size = shadersourcedata.size();
		Source.Encoding = DXC_CP_ACP; // Assume BOM (Byte Order Mark) says UTF8 or UTF16 or this is ANSI text. See http://en.wikipedia.org/wiki/Byte_order_mark

		// IDxcIncludeHandler is an interface for handling include directives.
		// To customize the handling of include directives, you can provide an implementation of this interface.
		struct IncludeHandler : public IDxcIncludeHandler
		{
			const CompilerInput* input = nullptr;
			CompilerOutput* output = nullptr;
			ComPtr<IDxcIncludeHandler> dxcIncludeHandler;

			// The LoadSource method is overrided.
			// It is called by the compiler to load the source for an include file.
			HRESULT STDMETHODCALLTYPE LoadSource(
				_In_z_ LPCWSTR pFilename,                                 // Candidate filename.
				_COM_Outptr_result_maybenull_ IDxcBlob** ppIncludeSource  // Resultant source object for included file, nullptr if not found.
			) override
			{
				HRESULT hr = dxcIncludeHandler->LoadSource(pFilename, ppIncludeSource);
				if (SUCCEEDED(hr))
				{
					// Save the path to the included file (usually globals.hlsli) in the dependencies list of the output
					std::string& filename = output->dependencies.emplace_back();
					wi::helper::StringConvert(pFilename, filename);
				}
				return hr;
			}
			HRESULT STDMETHODCALLTYPE QueryInterface(
				/* [in] */ REFIID riid,
				/* [iid_is][out] */ _COM_Outptr_ void __RPC_FAR* __RPC_FAR* ppvObject) override
			{
				return dxcIncludeHandler->QueryInterface(riid, ppvObject);
			}

			ULONG STDMETHODCALLTYPE AddRef(void) override
			{
				return 0;
			}
			ULONG STDMETHODCALLTYPE Release(void) override
			{
				return 0;
			}
		} includehandler;
		includehandler.input = &input;
		includehandler.output = &output;

		hr = dxcUtils->CreateDefaultIncludeHandler(&includehandler.dxcIncludeHandler);
		assert(SUCCEEDED(hr));

		wi::vector<const wchar_t*> args_raw; // constant version of args; that is, args_raw is a vector of pointers to constant strings
		args_raw.reserve(args.size());
		for (auto& x : args)
		{
			args_raw.push_back(x.c_str());
		}

		ComPtr<IDxcResult> pResults;
		hr = dxcCompiler->Compile(
			&Source,						// Source buffer.
			args_raw.data(),			// Array of pointers to arguments.
			(uint32_t)args.size(),		// Number of arguments.
			&includehandler,		// User-provided interface to handle #include directives (optional).
			IID_PPV_ARGS(&pResults)	// Compiler output status, buffer, and errors.
		);
		assert(SUCCEEDED(hr));

		ComPtr<IDxcBlobUtf8> pErrors = nullptr;
		hr = pResults->GetOutput(DXC_OUT_ERRORS, IID_PPV_ARGS(&pErrors), nullptr);
		assert(SUCCEEDED(hr));
		if (pErrors != nullptr && pErrors->GetStringLength() != 0)
		{
			output.error_message = pErrors->GetStringPointer();
		}

		HRESULT hrStatus;
		hr = pResults->GetStatus(&hrStatus);
		assert(SUCCEEDED(hr));
		if (FAILED(hrStatus))
		{
			return;
		}

		ComPtr<IDxcBlob> pShader = nullptr;
		hr = pResults->GetOutput(DXC_OUT_OBJECT, IID_PPV_ARGS(&pShader), nullptr);
		assert(SUCCEEDED(hr));
		if (pShader != nullptr)
		{
			output.dependencies.push_back(input.shadersourcefilename);
			output.shaderdata = (const uint8_t*)pShader->GetBufferPointer();
			output.shadersize = pShader->GetBufferSize();

			// keep the blob alive == keep shader pointer valid!
			auto internal_state = std::make_shared<ComPtr<IDxcBlob>>();
			*internal_state = pShader;
			output.internal_state = internal_state;
		}

		if (input.format == ShaderFormat::HLSL6)
		{
			// When a shader is compiled with HLSL6, a hash of the shader is stored in the shader.
			// The blob content returned using IDxcResult::GetOutput(DXC_OUT_SHADER_HASH, ...) maps to the DxcShaderHash structure
			ComPtr<IDxcBlob> pHash = nullptr;
			hr = pResults->GetOutput(DXC_OUT_SHADER_HASH, IID_PPV_ARGS(&pHash), nullptr);
			assert(SUCCEEDED(hr));
			if (pHash != nullptr)
			{
				DxcShaderHash* pHashBuf = (DxcShaderHash*)pHash->GetBufferPointer();
				for (int i = 0; i < _countof(pHashBuf->HashDigest); i++)
				{
					output.shaderhash.push_back(pHashBuf->HashDigest[i]);
				}
			}
		}
	}
#endif // SHADERCOMPILER_ENABLED_DXCOMPILER

#ifdef SHADERCOMPILER_ENABLED_D3DCOMPILER
	struct InternalState_D3DCompiler
	{
		using PFN_D3DCOMPILE = decltype(&D3DCompile);
		PFN_D3DCOMPILE D3DCompile = nullptr;

		InternalState_D3DCompiler()
		{
			if (D3DCompile != nullptr)
			{
				return; // already initialized
			}

			HMODULE d3dcompiler = wiLoadLibrary("d3dcompiler_47.dll");
			if (d3dcompiler != nullptr)
			{
				D3DCompile = (PFN_D3DCOMPILE)wiGetProcAddress(d3dcompiler, "D3DCompile");
				if (D3DCompile != nullptr)
				{
					wi::backlog::post("wi::shadercompiler: loaded d3dcompiler_47.dll");
				}
			}
		}
	};
	inline InternalState_D3DCompiler& d3d_compiler()
	{
		static InternalState_D3DCompiler internal_state;
		return internal_state;
	}

	void Compile_D3DCompiler(const CompilerInput& input, CompilerOutput& output)
	{
		if (d3d_compiler().D3DCompile == nullptr)
		{
			return;
		}

		if (input.minshadermodel > ShaderModel::SM_5_0)
		{
			output.error_message = "SHADERFORMAT_HLSL5 cannot support specified minshadermodel!";
			return;
		}

		wi::vector<uint8_t> shadersourcedata;
		if (!wi::helper::FileRead(input.shadersourcefilename, shadersourcedata))
		{
			return;
		}

		D3D_SHADER_MACRO defines[] = {
			"HLSL5", "1",
			"DISABLE_WAVE_INTRINSICS", "1",
			NULL, NULL,
		};

		const char* target = nullptr;
		switch (input.stage)
		{
		default:
		case ShaderStage::MS:
		case ShaderStage::AS:
		case ShaderStage::LIB:
			// not applicable
			return;
		case ShaderStage::VS:
			target = "vs_5_0";
			break;
		case ShaderStage::HS:
			target = "hs_5_0";
			break;
		case ShaderStage::DS:
			target = "ds_5_0";
			break;
		case ShaderStage::GS:
			target = "gs_5_0";
			break;
		case ShaderStage::PS:
			target = "ps_5_0";
			break;
		case ShaderStage::CS:
			target = "cs_5_0";
			break;
		}

		struct IncludeHandler : public ID3DInclude
		{
			const CompilerInput* input = nullptr;
			CompilerOutput* output = nullptr;
			wi::vector<wi::vector<uint8_t>> filedatas;

			HRESULT Open(D3D_INCLUDE_TYPE IncludeType, LPCSTR pFileName, LPCVOID pParentData, LPCVOID* ppData, UINT* pBytes) override
			{
				for (auto& x : input->include_directories)
				{
					std::string filename = x + pFileName;
					if (!wi::helper::FileExists(filename))
						continue;
					wi::vector<uint8_t>& filedata = filedatas.emplace_back();
					if (wi::helper::FileRead(filename, filedata))
					{
						output->dependencies.push_back(filename);
						*ppData = filedata.data();
						*pBytes = (UINT)filedata.size();
						return S_OK;
					}
				}
				return E_FAIL;
			}

			HRESULT Close(LPCVOID pData) override
			{
				return S_OK;
			}
		} includehandler;
		includehandler.input = &input;
		includehandler.output = &output;

		// https://docs.microsoft.com/en-us/windows/win32/direct3dhlsl/d3dcompile-constants
		UINT Flags1 = 0;
		if (has_flag(input.flags, Flags::DISABLE_OPTIMIZATION))
		{
			Flags1 |= D3DCOMPILE_SKIP_OPTIMIZATION;
		}


		ComPtr<ID3DBlob> code;
		ComPtr<ID3DBlob> errors;
		HRESULT hr = d3d_compiler().D3DCompile(
			shadersourcedata.data(),
			shadersourcedata.size(),
			input.shadersourcefilename.c_str(),
			defines,
			&includehandler, //D3D_COMPILE_STANDARD_FILE_INCLUDE,
			input.entrypoint.c_str(),
			target,
			Flags1,
			0,
			&code,
			&errors
		);

		if (errors)
		{
			output.error_message = (const char*)errors->GetBufferPointer();
		}

		if (SUCCEEDED(hr))
		{
			output.dependencies.push_back(input.shadersourcefilename);
			output.shaderdata = (const uint8_t*)code->GetBufferPointer();
			output.shadersize = code->GetBufferSize();

			// keep the blob alive == keep shader pointer valid!
			auto internal_state = std::make_shared<ComPtr<ID3D10Blob>>();
			*internal_state = code;
			output.internal_state = internal_state;
		}
	}
#endif // SHADERCOMPILER_ENABLED_D3DCOMPILER

	void Compile(const CompilerInput& input, CompilerOutput& output)
	{
		output = CompilerOutput();

#ifdef SHADERCOMPILER_ENABLED
		switch (input.format)
		{
		default:
			break;

#ifdef SHADERCOMPILER_ENABLED_DXCOMPILER
		case ShaderFormat::HLSL6:
		case ShaderFormat::SPIRV:
		case ShaderFormat::HLSL6_XS:
			Compile_DXCompiler(input, output);
			break;
#endif // SHADERCOMPILER_ENABLED_DXCOMPILER

#ifdef SHADERCOMPILER_ENABLED_D3DCOMPILER
		case ShaderFormat::HLSL5:
			Compile_D3DCompiler(input, output);
			break;
#endif // SHADERCOMPILER_ENABLED_D3DCOMPILER

#ifdef SHADERCOMPILER_PS5_INCLUDED
		case ShaderFormat::PS5:
			wi::shadercompiler::ps5::Compile(input, output);
			break;
#endif // SHADERCOMPILER_PS5_INCLUDED

		}
#endif // SHADERCOMPILER_ENABLED
	}

	constexpr const char* shadermetaextension = "wishadermeta";
	bool SaveShaderAndMetadata(const std::string& shaderfilename, const CompilerOutput& output)
	{
#ifdef SHADERCOMPILER_ENABLED
		wi::helper::DirectoryCreate(wi::helper::GetDirectoryFromPath(shaderfilename));

		wi::Archive dependencyLibrary(wi::helper::ReplaceExtension(shaderfilename, shadermetaextension), false);
		if (dependencyLibrary.IsOpen())
		{
			std::string rootdir = dependencyLibrary.GetSourceDirectory();
			wi::vector<std::string> dependencies = output.dependencies;
			for (auto& x : dependencies)
			{
				wi::helper::MakePathRelative(rootdir, x);
			}
			dependencyLibrary << dependencies;
		}

		if (wi::helper::FileWrite(shaderfilename, output.shaderdata, output.shadersize))
		{
			return true;
		}
#endif // SHADERCOMPILER_ENABLED

		return false;
	}
	bool IsShaderOutdated(const std::string& shaderfilename)
	{
#ifdef SHADERCOMPILER_ENABLED
		std::string filepath = shaderfilename;
		wi::helper::MakePathAbsolute(filepath);
		if (!wi::helper::FileExists(filepath))
		{
			return true; // no shader file = outdated shader, apps can attempt to rebuild it
		}

		// wishadermeta file contains the list of dependencies of the shader:
		// usually the hlsl file itself containing the shader code plus its included files like globals.hlsli, ShaderInterop.h, etc.
		// no wishadermeta file = no dependencies = no rebuild needed
		std::string dependencylibrarypath = wi::helper::ReplaceExtension(shaderfilename, shadermetaextension);
		if (!wi::helper::FileExists(dependencylibrarypath))
		{
			return false; // no metadata file = no dependency, up to date (for example packaged builds)
		}

		const uint64_t tim = wi::helper::FileTimestamp(filepath);

		wi::Archive dependencyLibrary(dependencylibrarypath);
		if (dependencyLibrary.IsOpen())
		{
			std::string rootdir = dependencyLibrary.GetSourceDirectory();
			wi::vector<std::string> dependencies;
			dependencyLibrary >> dependencies;

			for (auto& x : dependencies)
			{
				std::string dependencypath = rootdir + x;
				wi::helper::MakePathAbsolute(dependencypath);
				if (wi::helper::FileExists(dependencypath))
				{
					const uint64_t dep_tim = wi::helper::FileTimestamp(dependencypath);

					// if any dependency is newer than the shader, the shader is outdated
					if (tim < dep_tim)
					{
						return true;
					}
				}
			}
		}
#endif // SHADERCOMPILER_ENABLED

		return false;
	}

	std::mutex locker;
	wi::unordered_set<std::string> registered_shaders;
	void RegisterShader(const std::string& shaderfilename)
	{
#ifdef SHADERCOMPILER_ENABLED
		std::scoped_lock lock(locker);
		registered_shaders.insert(shaderfilename); // insert the shader name to the unordered set only if it's not already there
#endif // SHADERCOMPILER_ENABLED
	}
	size_t GetRegisteredShaderCount()
	{
		std::scoped_lock lock(locker);
		return registered_shaders.size();
	}
	bool CheckRegisteredShadersOutdated()
	{
#ifdef SHADERCOMPILER_ENABLED
		std::scoped_lock lock(locker);
		for (auto& x : registered_shaders)
		{
			if (IsShaderOutdated(x))
			{
				return true;
			}
		}
#endif // SHADERCOMPILER_ENABLED
		return false;
	}
}
