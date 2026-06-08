//==============================================================================
// ActionInitialization.cc
//==============================================================================
#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

ActionInitialization::ActionInitialization(const DetectorConstruction* dc)
: G4VUserActionInitialization(), fDetector(dc)
{}

void ActionInitialization::BuildForMaster() const
{
    SetUserAction(new RunAction());
}

void ActionInitialization::Build() const
{
    SetUserAction(new PrimaryGeneratorAction());

    RunAction*      ra = new RunAction();
    EventAction*    ea = new EventAction(ra);
    SteppingAction* sa = new SteppingAction(ea, fDetector);

    ea->SetSteppingAction(sa);

    SetUserAction(ra);
    SetUserAction(ea);
    SetUserAction(sa);
}
