#pragma once
// SteppingAction.hh
#include "G4UserSteppingAction.hh"
#include "NeutronChannel.hh"
#include <array>
#include <unordered_set>

class EventAction;
class DetectorConstruction;
class G4Step;
class G4Track;

// The process name registered in PhysicsListWithTALYS for 11B(p,n)
static const G4String kB11TALYSProcessName = "B11_pn_TALYS";

class SteppingAction : public G4UserSteppingAction {
public:
    SteppingAction(EventAction* ea, const DetectorConstruction* dc);
    ~SteppingAction() override;

    void UserSteppingAction(const G4Step* step) override;
    void NewEvent();
    void PrintSummary() const;

private:
    EventAction*                 fEventAction;
    const DetectorConstruction*  fDetector;

    std::unordered_set<G4int>    fRecordedTracks;
    G4long  fTotalSteps   = 0;
    G4long  fNeutronsSeen = 0;
    G4long  fAlphasSeen   = 0;
    G4int   fDiagPrinted  = 0;

    std::array<G4long, 8> fChannelNeutrons{};
    std::array<G4long, 8> fChannelAlphas{};

    NeutronChannel ClassifySecondary(const G4Track*, const G4Step*) const;
};
