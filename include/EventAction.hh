#pragma once
// EventAction.hh
#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "NeutronChannel.hh"
#include <array>

class RunAction;
class SteppingAction;
class G4Event;

class EventAction : public G4UserEventAction {
public:
    explicit EventAction(RunAction* ra);
    ~EventAction() override;

    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction  (const G4Event*) override;

    void RecordNeutron(G4double kineticE, G4double theta, G4double phi,
                       const G4ThreeVector& origin,
                       const G4String& process, NeutronChannel channel);
    void RecordAlpha  (G4double kineticE, NeutronChannel channel);

    void SetSteppingAction(SteppingAction* sa) { fSteppingAction = sa; }

private:
    RunAction*     fRunAction     = nullptr;
    SteppingAction* fSteppingAction = nullptr;

    G4int  fNeutronCount = 0;
    G4int  fAlphaCount   = 0;
    G4int  fDiagCount    = 0;
    std::array<G4int, 8> fChannelNeutrons{};
    std::array<G4int, 8> fChannelAlphas{};
};
