#pragma once
// DetectorConstruction.hh
#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    DetectorConstruction();
    ~DetectorConstruction() override;
    G4VPhysicalVolume* Construct()          override;
    void               ConstructSDandField() override;

    const G4LogicalVolume* GetCatcherLV()      const { return fCatcherLV; }
    const G4LogicalVolume* GetScoringSphLV()   const { return fScoringSphLV; }
    const G4Material*      GetCatcherMaterial() const { return fCatcherMat; }

private:
    void DefineMaterials();
    G4LogicalVolume* fCatcherLV    = nullptr;
    G4LogicalVolume* fScoringSphLV = nullptr;
    G4Material*      fCatcherMat   = nullptr;
    G4Material*      fWorldMat     = nullptr;
};
