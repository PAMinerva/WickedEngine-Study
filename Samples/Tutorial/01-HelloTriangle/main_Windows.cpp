#include "stdafx.h"
#include "sample.h"

#include "resource.h"

#define MAX_LOADSTRING 100

// Global Variables:
HINSTANCE hInst;                                // current instance
WCHAR szTitle[MAX_LOADSTRING];                  // The title bar text
WCHAR szWindowClass[MAX_LOADSTRING];            // the main window class name
SampleApp application;							// Application class

// Forward declarations of functions included in this code module:
BOOL                CreateWnd(WCHAR szTitle[], WCHAR szWindowClass[], HINSTANCE hInstance, int nCmdShow);
LRESULT CALLBACK    WndProc(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK    About(HWND, UINT, WPARAM, LPARAM);

_Use_decl_annotations_
int APIENTRY wWinMain(_In_ HINSTANCE hInstance,
                     _In_opt_ HINSTANCE hPrevInstance,
                     _In_ LPWSTR    lpCmdLine,
                     _In_ int       nCmdShow)
{
    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);

	/* Sets the current process to a specified dots per inch(dpi) awareness context to update desktop applications to handle
	display scale factor (dots per inch, or DPI) changes dynamically, allowing their applications to be crisp on any display they're rendered on.
	As display technology has progressed, display panel manufacturers have packed an increasing number of pixels into each unit of physical space on their panels.
	This has resulted in the dots per inch(DPI) of modern display panels being much higher than they have historically been.
	In the past, most displays had 96 pixels per linear inch of physical space(96 DPI); in 2017, displays with nearly 300 DPI or higher are readily available.
	Most legacy desktop UI frameworks have built - in assumptions that the display DPI will not change during the lifetime of the process.
	This assumption no longer holds true, with display DPIs commonly changing several times throughout an application process's lifetime.
	Some common scenarios where the display scale factor/DPI changes are:

	- Multiple - monitor setups where each display has a different scale factor and the application is moved from one display to another(such as a 4K and a 1080p display)
	- Docking and undocking a high DPI laptop with a low - DPI external display(or vice versa)
	- Connecting via Remote Desktop from a high DPI laptop/tablet to a low - DPI device (or vice versa)
	- Making a display - scale - factor settings change while applications are running

	In these scenarios, UWP applications redraw themselves for the new DPI automatically.
	By default, and without additional developer work, desktop applications do not.
	Desktop applications that don't do this extra work to respond to DPI changes may appear blurry or incorrectly-sized to the user.

	Desktop applications must tell Windows if they support DPI scaling. By default, the system considers desktop applications DPI unaware and bitmap-stretches their windows.
	By setting a DPI awareness mode, applications can explicitly tell Windows how they wish to handle DPI scaling.

	It is recommended that desktop applications be updated to use per-monitor DPI awareness mode, allowing them to immediately render correctly whenever the DPI changes.
	When an application reports to Windows that it wants to run in this mode, Windows will not bitmap stretch the application when the DPI changes, instead sending WM_DPICHANGED
	to the application window. It is then the complete responsibility of the application to handle resizing itself for the new DPI.
	Most UI frameworks used by desktop applications (Windows common controls (comctl32), Windows Forms, Windows Presentation Framework, etc.) do not support automatic DPI scaling,
	requiring developers to resize and reposition the contents of their windows themselves.
	There are two versions of Per-Monitor awareness that an application can register itself as: version 1 and version 2 (PMv2). Registering a process as running in PMv2 awareness mode results in:

	- The application being notified when the DPI changes (both the top-level and child HWNDs)
	- The application seeing the raw pixels of each display
	- The application never being bitmap scaled by Windows
	- Automatic non-client area (window caption, scroll bars, etc.) DPI scaling by Windows
	- Win32 dialogs (from CreateDialog) automatically DPI scaled by Windows
	- Theme-drawn bitmap assets in common controls (checkboxes, button backgrounds, etc.) being automatically rendered at the appropriate DPI scale factor

	When running in Per-Monitor v2 Awareness mode, applications are notified when their DPI has changed.
	If an application does not resize itself for the new DPI, the application UI will appear too small or too large (depending on the difference in the previous and new DPI values). */
    BOOL dpi_success = SetProcessDpiAwarenessContext(DPI_AWARENESS_CONTEXT_PER_MONITOR_AWARE_V2);
    assert(dpi_success);

	// Saves the lpCmdLine string into wi::arguments::params, which is an unordered_set<string>
    wi::arguments::Parse(lpCmdLine); // if you wish to use command line arguments, here is a good place to parse them...

    // Initialize global strings
    wi::helper::StringConvert("Hello Triangle", szTitle, MAX_LOADSTRING);
    wi::helper::StringConvert("LVWndClass", szWindowClass, MAX_LOADSTRING);

	// disable HDR rendering (even if supported by the display).
	// Otherwise, the HDR color space could be selected even if HDR is disabled in Windows display settings, which is not desirable.
	application.allow_hdr = false;

	// Create a window and call Application::SetWindow, which creates a graphics device and swapchain for the window
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
    application.infoDisplay.colorspace = true;
    application.infoDisplay.device_name = true;

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

BOOL CreateWnd(WCHAR szTitle[], WCHAR szWindowClass[], HINSTANCE hInstance, int nCmdShow)
{
	// Register a window class
	WNDCLASSEXW wcex{ 0 };

	wcex.cbSize = sizeof(WNDCLASSEXW);

	wcex.style = CS_HREDRAW | CS_VREDRAW;
	wcex.lpfnWndProc = WndProc;
	wcex.cbClsExtra = 0;
	wcex.cbWndExtra = 0;
	wcex.hInstance = hInstance;
	wcex.hIcon = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_TEMPLATEWINDOWS));
	wcex.hCursor = LoadCursor(nullptr, IDC_ARROW);
	wcex.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
	wcex.lpszMenuName = MAKEINTRESOURCEW(IDC_TEMPLATEWINDOWS);
	wcex.lpszClassName = szWindowClass;
	wcex.hIconSm = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

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

	// assign window handle (mandatory)
	// create a graphics device and swapchain for the window
	// set shader path to the shader folder containing the compiled shader binaries
	application.SetWindow(hWnd);

	return TRUE;
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
	case WM_COMMAND:
		{
			int wmId = LOWORD(wParam);
			// Parse the menu selections:
			switch (wmId)
			{
			case IDM_ABOUT:
				DialogBox(hInst, MAKEINTRESOURCE(IDD_ABOUTBOX), hWnd, About);
				break;
			case IDM_EXIT:
				DestroyWindow(hWnd);
				break;
			default:
				return DefWindowProc(hWnd, message, wParam, lParam);
			}
		}
		break;
    /*case WM_CREATE:
        {
            // Store application passed in to CreateWindow.
            LPCREATESTRUCT pCreateStruct = reinterpret_cast<LPCREATESTRUCT>(lParam);
            SetWindowLongPtr(hWnd, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(pCreateStruct->lpCreateParams));
        }*/
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

// Message handler for about box.
INT_PTR CALLBACK About(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
	UNREFERENCED_PARAMETER(lParam);
	switch (message)
	{
	case WM_INITDIALOG:
		return (INT_PTR)TRUE;

	case WM_COMMAND:
		if (LOWORD(wParam) == IDOK || LOWORD(wParam) == IDCANCEL)
		{
			EndDialog(hDlg, LOWORD(wParam));
			return (INT_PTR)TRUE;
		}
		break;
	}
	return (INT_PTR)FALSE;
}
