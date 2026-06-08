//==============================================================================
// main.cc — pBN fully hand-made multi-channel simulation
//==============================================================================
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "PhysicsListWithTALYS.hh"

#include <cstdlib>
#include <string>

int main(int argc, char** argv)
{
    G4UIExecutive* ui = nullptr;
    if (argc == 1) ui = new G4UIExecutive(argc, argv);

    auto* runManager = G4RunManagerFactory::CreateRunManager(
        G4RunManagerType::Default);

    int nThreads = 4;
    for (int i = 1; i < argc-1; ++i)
        if (std::string(argv[i]) == "-t")
            nThreads = std::atoi(argv[i+1]);
    runManager->SetNumberOfThreads(nThreads);

    auto* det = new DetectorConstruction();
    runManager->SetUserInitialization(det);

    // Single data directory — all channel XS + DDX subdirs live here
    runManager->SetUserInitialization(new PhysicsListWithTALYS("data"));
    runManager->SetUserInitialization(new ActionInitialization(det));

    auto* visManager = new G4VisExecutive("quiet");
    visManager->Initialize();

    auto* uiMgr = G4UImanager::GetUIpointer();
    if (ui) {
        uiMgr->ApplyCommand("/control/execute vis.mac");
        ui->SessionStart();
        delete ui;
    } else {
        uiMgr->ApplyCommand("/control/execute " + G4String(argv[argc-1]));
    }

    delete visManager;
    delete runManager;
    return 0;
}
