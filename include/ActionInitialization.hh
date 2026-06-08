#pragma once
// ActionInitialization.hh
#include "G4VUserActionInitialization.hh"
class DetectorConstruction;

class ActionInitialization : public G4VUserActionInitialization {
public:
    explicit ActionInitialization(const DetectorConstruction* dc);
    void BuildForMaster() const override;
    void Build()          const override;
private:
    const DetectorConstruction* fDetector;
};
