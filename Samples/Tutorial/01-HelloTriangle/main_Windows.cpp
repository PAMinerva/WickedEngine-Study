#include "stdafx.h"
#include "sample.h"

#define MAX_LOADSTRING 100

// Global Variables:
HINSTANCE hInst;                                // current instance
WCHAR szTitle[MAX_LOADSTRING];                  // The title bar text
WCHAR szWindowClass[MAX_LOADSTRING];            // the main window class name
wi::Application application;					// Application class

// Forward declarations of functions included in this code module:
BOOL                CreateWnd(WCHAR szTitle[], WCHAR szWindowClass[], HINSTANCE hInstance, int nCmdShow);
LRESULT CALLBACK    WndProc(HWND, UINT, WPARAM, LPARAM);

_Use_decl_annotations_
int APIENTRY wWinMain(_In_ HINSTANCE hInstance,
                     _In_opt_ HINSTANCE hPrevInstance,
                     _In_ LPWSTR    lpCmdLine,
                     _In_ int       nCmdShow)
{
    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);

    BOOL dpi_success = SetProcessDpiAwarenessContext(DPI_AWARENESS_CONTEXT_PER_MONITOR_AWARE_V2);
    assert(dpi_success);

    wi::arguments::Parse(lpCmdLine); // if you wish to use command line arguments, here is a good place to parse them...

    // Initialize global strings
    //wi::helper::StringConvert(sample.GetTitle(), szTitle, MAX_LOADSTRING);
    wi::helper::StringConvert("Hello Triangle", szTitle, MAX_LOADSTRING);
    wi::helper::StringConvert("LVWndClass", szWindowClass, MAX_LOADSTRING);

    // Create a window
    if (!CreateWnd(szTitle, szWindowClass, hInstance, nCmdShow))
    {
        return FALSE;
    }

    hInst = hInstance; // Store instance handle in our global variable

    // just show some basic info:
	application.infoDisplay.active = true;
	application.infoDisplay.watermark = true;
	application.infoDisplay.resolution = true;
	application.infoDisplay.fpsinfo = true;

    Sample sample;  // Sample class

    application.Initialize(); // Perform graphics initialization for the application (create device, swapchain, etc.)
    application.ActivatePath(sample.GetRenderPath3D());


 	MSG msg = { 0 };
	while (msg.message != WM_QUIT)
	{
		if (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE)) {
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
		else {

			application.Run(); // run the update - render loop (mandatory)

		}
	}

    return (int) msg.wParam;
}

//
//  FUNCTION: WndProc(HWND, UINT, WPARAM, LPARAM)
//
//  PURPOSE:  Processes messages for the main window.
//
//  WM_COMMAND  - process the application menu
//  WM_PAINT    - Paint the main window
//  WM_DESTROY  - post a quit message and return
//
//
LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    // Retrieve application from LPCREATESTRUCT
    //wi::Application* application = reinterpret_cast<wi::Application*>(GetWindowLongPtr(hWnd, GWLP_USERDATA));

    switch (message)
    {
    case WM_CREATE:
        {
            // Store application passed in to CreateWindow.
            LPCREATESTRUCT pCreateStruct = reinterpret_cast<LPCREATESTRUCT>(lParam);
            SetWindowLongPtr(hWnd, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(pCreateStruct->lpCreateParams));
        }
        return 0;
    case WM_SIZE:
    case WM_DPICHANGED:
        if (application.is_window_active)
            application.SetWindow(hWnd);
        break;
    case WM_KILLFOCUS:
        application.is_window_active = false;
        break;
    case WM_SETFOCUS:
        application.is_window_active = true;
        break;
    case WM_PAINT:
        {
            PAINTSTRUCT ps;
            HDC hdc = BeginPaint(hWnd, &ps);
            // TODO: Add any drawing code that uses hdc here...
            EndPaint(hWnd, &ps);
        }
        break;
    case WM_DESTROY:
        PostQuitMessage(0);
        break;
    default:
        return DefWindowProcW(hWnd, message, wParam, lParam);
    }
    return 0;
}

BOOL CreateWnd(WCHAR szTitle[], WCHAR szWindowClass[], HINSTANCE hInstance, int nCmdShow)
{
    // Register a window class
    WNDCLASSEXW wcex {0};

    wcex.cbSize = sizeof(WNDCLASSEXW);

    wcex.style          = CS_HREDRAW | CS_VREDRAW;
    wcex.lpfnWndProc    = WndProc;
    wcex.cbClsExtra     = 0;
    wcex.cbWndExtra     = 0;
    wcex.hInstance      = hInstance;
    //wcex.hIcon          = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_TEMPLATEWINDOWS));
    wcex.hCursor        = LoadCursor(nullptr, IDC_ARROW);
    wcex.hbrBackground  = (HBRUSH)(COLOR_WINDOW+1);
    //wcex.lpszMenuName   = MAKEINTRESOURCEW(IDC_TEMPLATEWINDOWS);
    wcex.lpszClassName  = szWindowClass;
    //wcex.hIconSm        = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

    if (!RegisterClassExW(&wcex))
    {
        DWORD error = GetLastError();
        return false;
    }

    // Create a window
    HWND hWnd = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW,
    CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, nullptr, nullptr, hInstance, nullptr);

    if (!hWnd)
    {
        return false;
    }

    ShowWindow(hWnd, nCmdShow);
    UpdateWindow(hWnd);


    application.SetWindow(hWnd); // assign window handle (mandatory)

    return TRUE;
}