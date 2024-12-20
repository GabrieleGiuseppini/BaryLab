/***************************************************************************************
 * Original Author:     Gabriele Giuseppini
 * Created:             2020-05-15
 * Copyright:           Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "MainFrame.h"

#include "StandardSystemPaths.h"

#include <GameOpenGL/GameOpenGL.h>

#include <GameCore/GameException.h>
#include <GameCore/Log.h>
#include <GameCore/Utils.h>
#include <GameCore/Version.h>

#include <wx/intl.h>
#include <wx/msgdlg.h>
#include <wx/panel.h>
#include <wx/settings.h>
#include <wx/sizer.h>
#include <wx/string.h>
#include <wx/tooltip.h>

#include <cassert>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <map>
#include <sstream>
#include <thread>

#ifdef _MSC_VER
 // Nothing to do here - we use RC files
#else
#include "Resources/SLabBBB.xpm"
#endif

long const ID_MAIN_CANVAS = wxNewId();

long const ID_LOAD_MESH_MENUITEM = wxNewId();
long const ID_RESET_MENUITEM = wxNewId();
long const ID_QUIT_MENUITEM = wxNewId();

long const ID_ZOOM_IN_MENUITEM = wxNewId();
long const ID_ZOOM_OUT_MENUITEM = wxNewId();
long const ID_RESET_VIEW_MENUITEM = wxNewId();

long const ID_OPEN_SETTINGS_WINDOW_MENUITEM = wxNewId();
long const ID_RELOAD_LAST_MODIFIED_SETTINGS_MENUITEM = wxNewId();
long const ID_OPEN_LOG_WINDOW_MENUITEM = wxNewId();
long const ID_FULL_SCREEN_MENUITEM = wxNewId();
long const ID_NORMAL_SCREEN_MENUITEM = wxNewId();

long const ID_ABOUT_MENUITEM = wxNewId();

long const ID_SIMULATION_TIMER = wxNewId();

MainFrame::MainFrame(wxApp * mainApp)
    : mIsMouseCapturedByGLCanvas(false)
    , mMainApp(mainApp)
    , mLabController()
    , mSettingsManager()
    , mToolController()
{
    Create(
        nullptr,
        wxID_ANY,
        std::string(APPLICATION_NAME_WITH_SHORT_VERSION),
        wxDefaultPosition,
        wxDefaultSize,
        wxDEFAULT_FRAME_STYLE | wxMAXIMIZE,
        _T("Main Frame"));

    SetIcon(wxICON(BBB_SLAB_ICON));
    SetBackgroundColour(wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE));
    Centre();

    Bind(wxEVT_CLOSE_WINDOW, &MainFrame::OnMainFrameClose, this);

    mMainPanel = new wxPanel(this, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxWANTS_CHARS);
    mMainPanel->Bind(wxEVT_CHAR_HOOK, &MainFrame::OnKeyDown, this);

    //
    // Build main GL canvas and activate GL context
    //

    // Note: Using the wxWidgets 3.1 style does not work on OpenGL 4 drivers; it forces a 1.1.0 context

    int mainGLCanvasAttributes[] =
    {
        WX_GL_RGBA,
        WX_GL_DOUBLEBUFFER,
        WX_GL_DEPTH_SIZE,      16,
        WX_GL_STENCIL_SIZE,    1,

        // Cannot specify CORE_PROFILE or else wx tries OpenGL 3.0 and fails if it's not supported
        //WX_GL_CORE_PROFILE,

        // Useless to specify version as Glad will always take the max
        //WX_GL_MAJOR_VERSION,    GameOpenGL::MinOpenGLVersionMaj,
        //WX_GL_MINOR_VERSION,    GameOpenGL::MinOpenGLVersionMin,

        0, 0
    };

    mMainGLCanvas = std::make_unique<wxGLCanvas>(
        mMainPanel,
        ID_MAIN_CANVAS,
        mainGLCanvasAttributes,
        wxDefaultPosition,
        wxDefaultSize,
        0L,
        _T("Main GL Canvas"));

    mMainGLCanvas->Connect(wxEVT_PAINT, (wxObjectEventFunction)&MainFrame::OnMainGLCanvasPaint, 0, this);
    mMainGLCanvas->Connect(wxEVT_SIZE, (wxObjectEventFunction)&MainFrame::OnMainGLCanvasResize, 0, this);
    mMainGLCanvas->Connect(wxEVT_LEFT_DOWN, (wxObjectEventFunction)&MainFrame::OnMainGLCanvasLeftDown, 0, this);
    mMainGLCanvas->Connect(wxEVT_LEFT_UP, (wxObjectEventFunction)&MainFrame::OnMainGLCanvasLeftUp, 0, this);
    mMainGLCanvas->Connect(wxEVT_RIGHT_DOWN, (wxObjectEventFunction)&MainFrame::OnMainGLCanvasRightDown, 0, this);
    mMainGLCanvas->Connect(wxEVT_RIGHT_UP, (wxObjectEventFunction)&MainFrame::OnMainGLCanvasRightUp, 0, this);
    mMainGLCanvas->Connect(wxEVT_MOTION, (wxObjectEventFunction)&MainFrame::OnMainGLCanvasMouseMove, 0, this);
    mMainGLCanvas->Connect(wxEVT_MOUSEWHEEL, (wxObjectEventFunction)&MainFrame::OnMainGLCanvasMouseWheel, 0, this);
    mMainGLCanvas->Connect(wxEVT_MOUSE_CAPTURE_LOST, (wxObjectEventFunction)&MainFrame::OnMainGLCanvasCaptureMouseLost, 0, this);

    // Take context for this canvas
    mMainGLCanvasContext = std::make_unique<wxGLContext>(mMainGLCanvas.get());

    // Activate context
    mMainGLCanvasContext->SetCurrent(*mMainGLCanvas);


    //
    // Layout panel
    //

    mMainPanelVSizer = new wxBoxSizer(wxVERTICAL);

    // Top
    {
        mMainPanelTopHSizer = new wxBoxSizer(wxHORIZONTAL);

        // Control toolbar
        mControlToolbar = new ControlToolbar(mMainPanel);
        mControlToolbar->Connect(ControlToolbar::ID_MOVE_PARTICLE, ControlToolbar::wxEVT_TOOLBAR_ACTION, (wxObjectEventFunction)&MainFrame::OnMoveParticle, 0, this);
        mControlToolbar->Connect(ControlToolbar::ID_SET_PARTICLE_TRAJECTORY, ControlToolbar::wxEVT_TOOLBAR_ACTION, (wxObjectEventFunction)&MainFrame::OnSetParticleTrajectory, 0, this);
        mControlToolbar->Connect(ControlToolbar::ID_SET_ORIGIN_TRIANGLE, ControlToolbar::wxEVT_TOOLBAR_ACTION, (wxObjectEventFunction)&MainFrame::OnSetOriginTriangle, 0, this);
        mControlToolbar->Connect(ControlToolbar::ID_SELECT_PARTICLE, ControlToolbar::wxEVT_TOOLBAR_ACTION, (wxObjectEventFunction)&MainFrame::OnSelectParticle, 0, this);
        mControlToolbar->Connect(ControlToolbar::ID_MOVE_VERTEX, ControlToolbar::wxEVT_TOOLBAR_ACTION, (wxObjectEventFunction)&MainFrame::OnMoveVertex, 0, this);
        mControlToolbar->Connect(ControlToolbar::ID_ROTATE_MESH_BY_POSITION, ControlToolbar::wxEVT_TOOLBAR_ACTION, (wxObjectEventFunction)&MainFrame::OnRotateMeshByPosition, 0, this);
        mControlToolbar->Connect(ControlToolbar::ID_ROTATE_MESH_BY_PARTICLE, ControlToolbar::wxEVT_TOOLBAR_ACTION, (wxObjectEventFunction)&MainFrame::OnRotateMeshByParticle, 0, this);

        mControlToolbar->Connect(ControlToolbar::ID_ADD_FURNITURE_PARTICLE_NPC, ControlToolbar::wxEVT_TOOLBAR_ACTION, (wxObjectEventFunction)&MainFrame::OnAddFurnitureParticleNpc, 0, this);
        mControlToolbar->Connect(ControlToolbar::ID_ADD_FURNITURE_QUAD_NPC, ControlToolbar::wxEVT_TOOLBAR_ACTION, (wxObjectEventFunction)&MainFrame::OnAddFurnitureQuadNpc, 0, this);
        mControlToolbar->Connect(ControlToolbar::ID_ADD_HUMAN_NPC, ControlToolbar::wxEVT_TOOLBAR_ACTION, (wxObjectEventFunction)&MainFrame::OnAddHumanNpc, 0, this);
        mControlToolbar->Connect(ControlToolbar::ID_MOVE_NPC, ControlToolbar::wxEVT_TOOLBAR_ACTION, (wxObjectEventFunction)&MainFrame::OnMoveNpc, 0, this);
        mControlToolbar->Connect(ControlToolbar::ID_REMOVE_NPC, ControlToolbar::wxEVT_TOOLBAR_ACTION, (wxObjectEventFunction)&MainFrame::OnRemoveNpc, 0, this);

        mControlToolbar->Connect(ControlToolbar::ID_SIMULATION_CONTROL_PLAY, ControlToolbar::wxEVT_TOOLBAR_ACTION, (wxObjectEventFunction)&MainFrame::OnSimulationControlPlay, 0, this);
        mControlToolbar->Connect(ControlToolbar::ID_SIMULATION_CONTROL_PAUSE, ControlToolbar::wxEVT_TOOLBAR_ACTION, (wxObjectEventFunction)&MainFrame::OnSimulationControlPause, 0, this);
        mControlToolbar->Connect(ControlToolbar::ID_SIMULATION_CONTROL_STEP, ControlToolbar::wxEVT_TOOLBAR_ACTION, (wxObjectEventFunction)&MainFrame::OnSimulationControlStep, 0, this);
        mControlToolbar->Connect(ControlToolbar::ID_ACTION_RESET, ControlToolbar::wxEVT_TOOLBAR_ACTION, (wxObjectEventFunction)&MainFrame::OnResetMenuItemSelected, 0, this);
        mControlToolbar->Connect(ControlToolbar::ID_ACTION_LOAD_MESH, ControlToolbar::wxEVT_TOOLBAR_ACTION, (wxObjectEventFunction)&MainFrame::OnLoadMeshMenuItemSelected, 0, this);
        mControlToolbar->Connect(ControlToolbar::ID_ACTION_SETTINGS, ControlToolbar::wxEVT_TOOLBAR_ACTION, (wxObjectEventFunction)&MainFrame::OnOpenSettingsWindowMenuItemSelected, 0, this);
        mControlToolbar->Connect(ControlToolbar::ID_VIEW_CONTROL_GRID, ControlToolbar::wxEVT_TOOLBAR_ACTION, (wxObjectEventFunction)&MainFrame::OnViewControlGridToggled, 0, this);
        mControlToolbar->Bind(EVT_MESH_TRANSFORMATION_CHANGED, &MainFrame::OnMeshTransformationChanged, this);
        mControlToolbar->Bind(EVT_HUMAN_NPC_PANIC_LEVEL_CHANGED, &MainFrame::OnHumanNpcPanicLevelChanged, this);

        mMainPanelTopHSizer->Add(
            mControlToolbar,
            0,                  // Use own horizontal size
            wxEXPAND,           // Expand vertically
            0);                 // Border

        // Canvas
        mMainPanelTopHSizer->Add(
            mMainGLCanvas.get(),
            1,                  // Occupy all available horizontal space
            wxEXPAND,           // Expand also vertically
            0);                 // Border

        mMainPanelVSizer->Add(
            mMainPanelTopHSizer,
            1,                  // Occupy all available vertical space
            wxEXPAND,           // Expand also horizontally
            0);                 // Border
    }

    // Bottom
    {
        // Probe toolbar
        mProbeToolbar = new ProbeToolbar(mMainPanel);

        mMainPanelVSizer->Add(
            mProbeToolbar,
            0,                  // Own height
            wxEXPAND);          // Expand horizontally
    }

    //
    // Build menu
    //

    {
        wxMenuBar * mainMenuBar = new wxMenuBar();


        // File

        wxMenu * fileMenu = new wxMenu();

        wxMenuItem * loadMeshMenuItem = new wxMenuItem(fileMenu, ID_LOAD_MESH_MENUITEM, _("Load Mesh\tCtrl+O"), wxEmptyString, wxITEM_NORMAL);
        fileMenu->Append(loadMeshMenuItem);
        Connect(ID_LOAD_MESH_MENUITEM, wxEVT_COMMAND_MENU_SELECTED, (wxObjectEventFunction)&MainFrame::OnLoadMeshMenuItemSelected);

        wxMenuItem * resetMenuItem = new wxMenuItem(fileMenu, ID_RESET_MENUITEM, _("Reset\tCtrl+R"), wxEmptyString, wxITEM_NORMAL);
        fileMenu->Append(resetMenuItem);
        Connect(ID_RESET_MENUITEM, wxEVT_COMMAND_MENU_SELECTED, (wxObjectEventFunction)&MainFrame::OnResetMenuItemSelected);

        fileMenu->Append(new wxMenuItem(fileMenu, wxID_SEPARATOR));

        wxMenuItem * quitMenuItem = new wxMenuItem(fileMenu, ID_QUIT_MENUITEM, _("Quit\tAlt-F4"), _("Quit the application"), wxITEM_NORMAL);
        fileMenu->Append(quitMenuItem);
        Connect(ID_QUIT_MENUITEM, wxEVT_COMMAND_MENU_SELECTED, (wxObjectEventFunction)&MainFrame::OnQuit);

        mainMenuBar->Append(fileMenu, _("&File"));


        // Controls

        wxMenu * controlsMenu = new wxMenu();

        wxMenuItem * zoomInMenuItem = new wxMenuItem(controlsMenu, ID_ZOOM_IN_MENUITEM, _("Zoom In\t+"), wxEmptyString, wxITEM_NORMAL);
        controlsMenu->Append(zoomInMenuItem);
        Connect(ID_ZOOM_IN_MENUITEM, wxEVT_COMMAND_MENU_SELECTED, (wxObjectEventFunction)&MainFrame::OnZoomInMenuItemSelected);

        wxMenuItem * zoomOutMenuItem = new wxMenuItem(controlsMenu, ID_ZOOM_OUT_MENUITEM, _("Zoom Out\t-"), wxEmptyString, wxITEM_NORMAL);
        controlsMenu->Append(zoomOutMenuItem);
        Connect(ID_ZOOM_OUT_MENUITEM, wxEVT_COMMAND_MENU_SELECTED, (wxObjectEventFunction)&MainFrame::OnZoomOutMenuItemSelected);

        controlsMenu->Append(new wxMenuItem(controlsMenu, wxID_SEPARATOR));

        wxMenuItem * resetViewMenuItem = new wxMenuItem(controlsMenu, ID_RESET_VIEW_MENUITEM, _("Reset View\tHOME"), wxEmptyString, wxITEM_NORMAL);
        controlsMenu->Append(resetViewMenuItem);
        Connect(ID_RESET_VIEW_MENUITEM, wxEVT_COMMAND_MENU_SELECTED, (wxObjectEventFunction)&MainFrame::OnResetViewMenuItemSelected);

        mainMenuBar->Append(controlsMenu, _("Controls"));


        // Options

        wxMenu * optionsMenu = new wxMenu();

        wxMenuItem * openSettingsWindowMenuItem = new wxMenuItem(optionsMenu, ID_OPEN_SETTINGS_WINDOW_MENUITEM, _("Simulation Settings...\tCtrl+S"), wxEmptyString, wxITEM_NORMAL);
        optionsMenu->Append(openSettingsWindowMenuItem);
        Connect(ID_OPEN_SETTINGS_WINDOW_MENUITEM, wxEVT_COMMAND_MENU_SELECTED, (wxObjectEventFunction)&MainFrame::OnOpenSettingsWindowMenuItemSelected);

        mReloadLastModifiedSettingsMenuItem = new wxMenuItem(optionsMenu, ID_RELOAD_LAST_MODIFIED_SETTINGS_MENUITEM, _("Reload Last-Modified Simulation Settings\tCtrl+D"), wxEmptyString, wxITEM_NORMAL);
        optionsMenu->Append(mReloadLastModifiedSettingsMenuItem);
        Connect(ID_RELOAD_LAST_MODIFIED_SETTINGS_MENUITEM, wxEVT_COMMAND_MENU_SELECTED, (wxObjectEventFunction)&MainFrame::OnReloadLastModifiedSettingsMenuItem);

        optionsMenu->Append(new wxMenuItem(optionsMenu, wxID_SEPARATOR));

        wxMenuItem * openLogWindowMenuItem = new wxMenuItem(optionsMenu, ID_OPEN_LOG_WINDOW_MENUITEM, _("Open Log Window\tCtrl+L"), wxEmptyString, wxITEM_NORMAL);
        optionsMenu->Append(openLogWindowMenuItem);
        Connect(ID_OPEN_LOG_WINDOW_MENUITEM, wxEVT_COMMAND_MENU_SELECTED, (wxObjectEventFunction)&MainFrame::OnOpenLogWindowMenuItemSelected);

        optionsMenu->Append(new wxMenuItem(optionsMenu, wxID_SEPARATOR));

        mFullScreenMenuItem = new wxMenuItem(optionsMenu, ID_FULL_SCREEN_MENUITEM, _("Full Screen\tF11"), wxEmptyString, wxITEM_NORMAL);
        optionsMenu->Append(mFullScreenMenuItem);
        Connect(ID_FULL_SCREEN_MENUITEM, wxEVT_COMMAND_MENU_SELECTED, (wxObjectEventFunction)&MainFrame::OnFullScreenMenuItemSelected);
        mFullScreenMenuItem->Enable(!StartInFullScreenMode);

        mNormalScreenMenuItem = new wxMenuItem(optionsMenu, ID_NORMAL_SCREEN_MENUITEM, _("Normal Screen\tESC"), wxEmptyString, wxITEM_NORMAL);
        optionsMenu->Append(mNormalScreenMenuItem);
        Connect(ID_NORMAL_SCREEN_MENUITEM, wxEVT_COMMAND_MENU_SELECTED, (wxObjectEventFunction)&MainFrame::OnNormalScreenMenuItemSelected);
        mNormalScreenMenuItem->Enable(StartInFullScreenMode);

        mainMenuBar->Append(optionsMenu, _("Options"));


        // Help

        wxMenu * helpMenu = new wxMenu();

        wxMenuItem * aboutMenuItem = new wxMenuItem(helpMenu, ID_ABOUT_MENUITEM, _("About\tF2"), _("Show credits and other I'vedunnit stuff"), wxITEM_NORMAL);
        helpMenu->Append(aboutMenuItem);
        Connect(ID_ABOUT_MENUITEM, wxEVT_COMMAND_MENU_SELECTED, (wxObjectEventFunction)&MainFrame::OnAboutMenuItemSelected);

        mainMenuBar->Append(helpMenu, _("Help"));

        SetMenuBar(mainMenuBar);
    }

    //
    // Finalize frame
    //

    {
        mMainPanel->SetSizer(mMainPanelVSizer);
        mMainPanel->Layout();

        auto * wholeSizer = new wxBoxSizer(wxVERTICAL);
        wholeSizer->Add(mMainPanel, 1, wxEXPAND, 0);
        this->SetSizer(wholeSizer);
    }


    //
    // Initialize tooltips
    //

    wxToolTip::Enable(true);
    wxToolTip::SetDelay(200);


    //
    // PostInitialize
    //

    Show(true);

    if (StartInFullScreenMode)
        ShowFullScreen(true, wxFULLSCREEN_NOBORDER);


    //
    // Initialize timers
    //

    mSimulationTimer = std::make_unique<wxTimer>(this, ID_SIMULATION_TIMER);
    Connect(ID_SIMULATION_TIMER, wxEVT_TIMER, (wxObjectEventFunction)&MainFrame::OnSimulationTimer);
    mSimulationTimer->Start(0);
}

MainFrame::~MainFrame()
{
}

//
// App event handlers
//

void MainFrame::OnMainFrameClose(wxCloseEvent & /*event*/)
{
    if (!!mSimulationTimer)
        mSimulationTimer->Stop();

    if (!!mSettingsManager)
        mSettingsManager->SaveLastModifiedSettings();

    Destroy();
}

void MainFrame::OnQuit(wxCommandEvent & /*event*/)
{
    Close();
}

void MainFrame::OnKeyDown(wxKeyEvent & event)
{
    assert(!!mLabController);

    if (event.GetKeyCode() == WXK_LEFT)
    {
        // Left
        mLabController->Pan(vec2f(-20.0, 0.0f));
    }
    else if (event.GetKeyCode() == WXK_UP)
    {
        // Up
        mLabController->Pan(vec2f(00.0f, -20.0f));
    }
    else if (event.GetKeyCode() == WXK_RIGHT)
    {
        // Right
        mLabController->Pan(vec2f(20.0f, 0.0f));
    }
    else if (event.GetKeyCode() == WXK_DOWN)
    {
        // Down
        mLabController->Pan(vec2f(0.0f, 20.0f));
    }
    else if (event.GetKeyCode() == WXK_DELETE)
    {
        auto dummy = wxCommandEvent(wxEVT_NULL);
        OnRemoveNpc(dummy);
        mControlToolbar->ReconciliateUIWithTool(mToolController->GetTool());
    }
    else if (event.GetKeyCode() == '3')
    {
        auto dummy = wxCommandEvent(wxEVT_NULL);
        OnMoveNpc(dummy);
        mControlToolbar->ReconciliateUIWithTool(mToolController->GetTool());
    }
    else if (event.GetKeyCode() == '/')
    {
        // Query

        assert(!!mToolController);

        vec2f screenCoords = mToolController->GetMouseScreenCoordinates();
        vec2f worldCoords = mLabController->ScreenToWorld(screenCoords);

        LogMessage(worldCoords.toString(), ":");

        //mLabController->QueryNearestNpcParticleAt(screenCoords);
        mLabController->QueryPointAt(screenCoords);
    }
    else if (event.GetKeyCode() == 'F')
    {
        // Flip

        mLabController->FlipCurrentlySelectedHuman();
    }
    else if (event.GetKeyCode() == 'V')
    {
        // Video

        mLabController->DoStepForVideo();
    }
    else if (mControlToolbar->ProcessKeyDown(event.GetKeyCode(), event.GetModifiers()))
    {
        // Processed
    }
    else
    {
        // Keep processing it
        event.Skip();
    }
}

//
// Main canvas event handlers
//

void MainFrame::OnMainGLCanvasPaint(wxPaintEvent & event)
{
    if (mLabController)
    {
        mLabController->Render();

        assert(mMainGLCanvas);
        mMainGLCanvas->SwapBuffers();
    }

    event.Skip();
}

void MainFrame::OnMainGLCanvasResize(wxSizeEvent & event)
{
    LogMessage("OnMainGLCanvasResize: ", event.GetSize().GetX(), "x", event.GetSize().GetY());

    if (mLabController
        && event.GetSize().GetX() > 0
        && event.GetSize().GetY() > 0)
    {
        mLabController->SetCanvasSize(
            event.GetSize().GetX(),
            event.GetSize().GetY());
    }

    if (mProbeToolbar)
    {
        mProbeToolbar->Refresh();
    }

    event.Skip();
}

void MainFrame::OnMainGLCanvasLeftDown(wxMouseEvent & /*event*/)
{
    if (mToolController)
    {
        mToolController->OnLeftMouseDown();
    }

    // Hang on to the mouse for as long as the button is pressed
    if (!mIsMouseCapturedByGLCanvas)
    {
        mMainGLCanvas->CaptureMouse();
        mIsMouseCapturedByGLCanvas = true;
    }
}

void MainFrame::OnMainGLCanvasLeftUp(wxMouseEvent & /*event*/)
{
    // We can now release the mouse
    if (mIsMouseCapturedByGLCanvas)
    {
        mMainGLCanvas->ReleaseMouse();
        mIsMouseCapturedByGLCanvas = false;
    }

    if (mToolController)
    {
        mToolController->OnLeftMouseUp();
    }
}

void MainFrame::OnMainGLCanvasRightDown(wxMouseEvent & /*event*/)
{
    if (mToolController)
    {
        mToolController->OnRightMouseDown();
    }

    // Hang on to the mouse for as long as the button is pressed
    if (!mIsMouseCapturedByGLCanvas)
    {
        mMainGLCanvas->CaptureMouse();
        mIsMouseCapturedByGLCanvas = true;
    }
}

void MainFrame::OnMainGLCanvasRightUp(wxMouseEvent & /*event*/)
{
    // We can now release the mouse
    if (mIsMouseCapturedByGLCanvas)
    {
        mMainGLCanvas->ReleaseMouse();
        mIsMouseCapturedByGLCanvas = false;
    }

    if (mToolController)
    {
        mToolController->OnRightMouseUp();
    }
}

void MainFrame::OnMainGLCanvasMouseMove(wxMouseEvent & event)
{
    if (mToolController)
    {
        mToolController->OnMouseMove(event.GetX(), event.GetY());
    }
}

void MainFrame::OnMainGLCanvasMouseWheel(wxMouseEvent & event)
{
    if (mLabController)
    {
        mLabController->AdjustZoom(powf(1.002f, event.GetWheelRotation()));
    }
}

void MainFrame::OnMainGLCanvasCaptureMouseLost(wxMouseCaptureLostEvent & /*event*/)
{
    if (mToolController)
    {
        mToolController->UnsetTool();
    }
}

//
// Menu event handlers
//

void MainFrame::OnLoadMeshMenuItemSelected(wxCommandEvent & /*event*/)
{
    if (!mFileOpenDialog)
    {
        mFileOpenDialog = std::make_unique<wxFileDialog>(
            this,
            L"Select Object",
            ResourceLocator::GetInstalledMeshesFolderPath().string(),
            wxEmptyString,
            L"Object files (*.png)|*.png",
            wxFD_OPEN | wxFD_FILE_MUST_EXIST,
            wxDefaultPosition,
            wxDefaultSize,
            _T("File Open Dialog"));
    }

    assert(!!mFileOpenDialog);

    if (mFileOpenDialog->ShowModal() == wxID_OK)
    {
        std::string const filepath = mFileOpenDialog->GetPath().ToStdString();

        assert(!!mLabController);
        try
        {
            mLabController->LoadShip(filepath);
        }
        catch (std::exception const & ex)
        {
            OnError(ex.what(), false);
        }
    }
}

void MainFrame::OnResetMenuItemSelected(wxCommandEvent & /*event*/)
{
    assert(!!mLabController);
    mLabController->Reset();

    mControlToolbar->ReconciliateUI(
        mLabController->GetSimulationControlState(),
        mLabController->IsViewGridEnabled());
}

void MainFrame::OnResetViewMenuItemSelected(wxCommandEvent & /*event*/)
{
    assert(!!mLabController);

    mLabController->ResetPan();
    mLabController->ResetZoom();
}

void MainFrame::OnZoomInMenuItemSelected(wxCommandEvent & /*event*/)
{
    assert(!!mLabController);
    mLabController->AdjustZoom(1.05f);
}

void MainFrame::OnZoomOutMenuItemSelected(wxCommandEvent & /*event*/)
{
    assert(!!mLabController);
    mLabController->AdjustZoom(1.0f / 1.05f);
}

////////////////////////////////////////////////////////////////////////////

void MainFrame::OnOpenSettingsWindowMenuItemSelected(wxCommandEvent & /*event*/)
{
    if (!mSettingsDialog)
    {
        mSettingsDialog = std::make_unique<SettingsDialog>(
            this,
            mSettingsManager,
            mLabController);
    }

    mSettingsDialog->Open();
}

void MainFrame::OnReloadLastModifiedSettingsMenuItem(wxCommandEvent & /*event*/)
{
    // Load last-modified settings
    try
    {
        assert(!!mSettingsManager);
        mSettingsManager->EnforceDefaultsAndLastModifiedSettings();
    }
    catch (std::exception const & exc)
    {
        OnError("Could not load last-modified settings: " + std::string(exc.what()), false);

        // Disable menu item
        mReloadLastModifiedSettingsMenuItem->Enable(false);
    }
}

void MainFrame::OnOpenLogWindowMenuItemSelected(wxCommandEvent & /*event*/)
{
    if (!mLoggingDialog)
    {
        mLoggingDialog = std::make_unique<LoggingDialog>(this);
    }

    mLoggingDialog->Open();
}

void MainFrame::OnFullScreenMenuItemSelected(wxCommandEvent & /*event*/)
{
    mFullScreenMenuItem->Enable(false);
    mNormalScreenMenuItem->Enable(true);

    this->ShowFullScreen(true, wxFULLSCREEN_NOBORDER);
}

void MainFrame::OnNormalScreenMenuItemSelected(wxCommandEvent & /*event*/)
{
    mFullScreenMenuItem->Enable(true);
    mNormalScreenMenuItem->Enable(false);

    this->ShowFullScreen(false);
}

void MainFrame::OnAboutMenuItemSelected(wxCommandEvent & /*event*/)
{
    if (!mAboutDialog)
    {
        mAboutDialog = std::make_unique<AboutDialog>(this);
    }

    mAboutDialog->Open();
}

void MainFrame::OnMoveParticle(wxCommandEvent & /*event*/)
{
    assert(!!mToolController);
    mToolController->SetTool(ToolType::MoveParticle);
}

void MainFrame::OnSetParticleTrajectory(wxCommandEvent & /*event*/)
{
    assert(!!mToolController);
    mToolController->SetTool(ToolType::SetParticleTrajectory);
}

void MainFrame::OnSetOriginTriangle(wxCommandEvent & /*event*/)
{
    assert(!!mToolController);
    mToolController->SetTool(ToolType::SetOriginTriangle);
}

void MainFrame::OnSelectParticle(wxCommandEvent & /*event*/)
{
    assert(!!mToolController);
    mToolController->SetTool(ToolType::SelectParticle);
}

void MainFrame::OnMoveVertex(wxCommandEvent & /*event*/)
{
    assert(!!mToolController);
    mToolController->SetTool(ToolType::MoveVertex);
}

void MainFrame::OnRotateMeshByPosition(wxCommandEvent & /*event*/)
{
    assert(!!mToolController);
    mToolController->SetTool(ToolType::RotateMeshByPosition);
}

void MainFrame::OnRotateMeshByParticle(wxCommandEvent & /*event*/)
{
    assert(!!mToolController);
    mToolController->SetTool(ToolType::RotateMeshByParticle);
}

void MainFrame::OnAddFurnitureParticleNpc(wxCommandEvent & /*event*/)
{
    assert(!!mToolController);
    mToolController->SetTool(ToolType::AddFurnitureParticleNpc);
}

void MainFrame::OnAddFurnitureQuadNpc(wxCommandEvent & /*event*/)
{
    assert(!!mToolController);
    mToolController->SetTool(ToolType::AddFurnitureQuadNpc);
}

void MainFrame::OnAddHumanNpc(wxCommandEvent & /*event*/)
{
    assert(!!mToolController);
    mToolController->SetHumanNpcPlaceTool(0); // Whatever is this NPC
}

void MainFrame::OnMoveNpc(wxCommandEvent & /*event*/)
{
    assert(!!mToolController);
    mToolController->SetTool(ToolType::MoveNpc);
}

void MainFrame::OnRemoveNpc(wxCommandEvent & /*event*/)
{
    assert(!!mToolController);
    mToolController->SetTool(ToolType::RemoveNpc);
}

void MainFrame::OnSimulationControlPlay(wxCommandEvent & /*event*/)
{
    assert(!!mLabController);
    mLabController->SetSimulationControlState(SimulationControlStateType::Play);
}

void MainFrame::OnSimulationControlPause(wxCommandEvent & /*event*/)
{
    assert(!!mLabController);
    mLabController->SetSimulationControlState(SimulationControlStateType::Paused);
}

void MainFrame::OnSimulationControlStep(wxCommandEvent & /*event*/)
{
    assert(!!mLabController);
    mLabController->SetSimulationControlPulse();
}

void MainFrame::OnViewControlGridToggled(wxCommandEvent & event)
{
    assert(!!mLabController);
    mLabController->SetViewGridEnabled(event.GetInt() != 0);
}

void MainFrame::OnMeshTransformationChanged(ControlToolbar::meshTransformationChangedEvent & event)
{
    assert(!!mLabController);
    mLabController->SetShipVelocity(event.GetMeshVelocity());
    mLabController->SetWavesAmplitude(event.GetWavesAmplitude());
    mLabController->SetWavesSpeed(event.GetWavesSpeed());
}

void MainFrame::OnHumanNpcPanicLevelChanged(ControlToolbar::humanNpcPanicLevelChangedEvent & event)
{
    assert(!!mLabController);
    mLabController->SetNpcPanicLevelForAllHumans(event.GetPanicLevel());
}

void MainFrame::OnSimulationTimer(wxTimerEvent & /*event*/)
{
    //
    // Complete initialization, if still not done
    //
    // We do it here to make sure we get the final canvas size
    //

    if (!mLabController)
    {
        try
        {
            FinishInitialization();
        }
        catch (GameException const & e)
        {
            mSimulationTimer->Stop(); // Stop looping and allow Die() to finish

            OnError(std::string(e.what()), true);
            return;
        }
    }

    //
    // Update tools
    //

    assert(!!mToolController);

    if (wxGetKeyState(WXK_SHIFT))
    {
        if (!mToolController->IsShiftKeyDown())
            mToolController->OnShiftKeyDown();
    }
    else
    {
        if (mToolController->IsShiftKeyDown())
            mToolController->OnShiftKeyUp();
    }

    mToolController->Update();


    //
    // Update
    //

    assert(!!mLabController);

    mLabController->Update();


    //
    // Render
    //

    assert(!!mMainGLCanvas);

    mMainGLCanvas->Refresh();


    //
    // Publish perf
    //

    mLabController->PublishPerf();


    //
    // Update probe toolbar
    //

    assert(!!mProbeToolbar);

    mProbeToolbar->Update();

    // TODOTEST
////#ifndef _DEBUG
////    static int TODO = 0;
////    ++TODO;
////    if (TODO > 10 * 64)
////    {
////        Close();
////    }
////#endif
}

/////////////////////////////////////////////////////////////////////////////////////////////////

void MainFrame::FinishInitialization()
{
    //
    // Create Simulation Controller
    //

    assert(!mLabController);

    try
    {
        mLabController = LabController::Create(
            mMainGLCanvas->GetSize().x,
            mMainGLCanvas->GetSize().y);
    }
    catch (std::exception const & e)
    {
        throw GameException("Error during initialization of simulation controller: " + std::string(e.what()));
    }

    //
    // Create Settings Manager
    //

    mSettingsManager = std::make_shared<SettingsManager>(
        mLabController,
        StandardSystemPaths::GetInstance().GetUserSettingsRootFolderPath());

    // Enable "Reload Last Modified Settings" menu if we have last-modified settings
    mReloadLastModifiedSettingsMenuItem->Enable(mSettingsManager->HasLastModifiedSettingsPersisted());

    //
    // Create Tool Controller
    //

    ToolType constexpr InitialToolType = ToolType::AddHumanNpc;

    try
    {
        mToolController = std::make_unique<ToolController>(
            InitialToolType,
            mMainGLCanvas.get(),
            mLabController);
    }
    catch (std::exception const & e)
    {
        throw GameException("Error during initialization of tool controller: " + std::string(e.what()));
    }

    mControlToolbar->ReconciliateUIWithTool(InitialToolType);

    //
    // Register event handlers
    //

    mLabController->RegisterBLabEventHandler(mProbeToolbar);
    mLabController->RegisterNpcGameEventHandler(mProbeToolbar);

    //
    // Load initial ship
    //

    mLabController->LoadShip(ResourceLocator::GetDefaultMeshDefinitionFilePath());

    //
    // Reconciliate UI
    //

    mControlToolbar->ReconciliateUI(
        mLabController->GetSimulationControlState(),
        mLabController->IsViewGridEnabled());
}

void MainFrame::OnError(
    std::string const & message,
    bool die)
{
    //
    // Show message
    //

    wxMessageBox(message, wxT("Simulation Disaster"), wxICON_ERROR);

    if (die)
    {
        //
        // Exit
        //

        this->Destroy();
    }
}