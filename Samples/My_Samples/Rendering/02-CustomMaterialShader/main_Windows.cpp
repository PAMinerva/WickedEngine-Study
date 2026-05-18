#include "stdafx.h"
#include "sample.h"

#include <windows.h>

namespace
{
	SampleApp app;

	LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
	{
		switch (message)
		{
		case WM_SIZE:
		case WM_DPICHANGED:
			if (app.is_window_active)
			{
				app.SetWindow(hWnd);
			}
			break;
		case WM_INPUT:
			wi::input::rawinput::ParseMessage(reinterpret_cast<void*>(lParam));
			break;
		case WM_KILLFOCUS:
			app.is_window_active = false;
			break;
		case WM_SETFOCUS:
			app.is_window_active = true;
			break;
		case WM_DESTROY:
			PostQuitMessage(0);
			break;
		default:
			return DefWindowProc(hWnd, message, wParam, lParam);
		}
		return 0;
	}
}

int APIENTRY wWinMain(HINSTANCE hInstance, HINSTANCE, LPWSTR commandLine, int nCmdShow)
{
	wi::arguments::Parse(commandLine);

	const wchar_t* className = L"CustomMaterialShaderSample";
	WNDCLASSEXW wc = {};
	wc.cbSize = sizeof(WNDCLASSEXW);
	wc.style = CS_HREDRAW | CS_VREDRAW;
	wc.lpfnWndProc = WndProc;
	wc.hInstance = hInstance;
	wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
	wc.lpszClassName = className;
	RegisterClassExW(&wc);

	HWND hWnd = CreateWindowW(
		className,
		L"02 - Custom HLSL Material Shader",
		WS_OVERLAPPEDWINDOW,
		CW_USEDEFAULT,
		0,
		1280,
		800,
		nullptr,
		nullptr,
		hInstance,
		nullptr);

	if (hWnd == nullptr)
	{
		return 1;
	}

	ShowWindow(hWnd, nCmdShow);
	UpdateWindow(hWnd);

	app.allow_hdr = false;
	app.SetWindow(hWnd);

	MSG msg = {};
	while (msg.message != WM_QUIT)
	{
		if (PeekMessage(&msg, nullptr, 0, 0, PM_REMOVE))
		{
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
		else
		{
			app.Run();
		}
	}

	return static_cast<int>(msg.wParam);
}
