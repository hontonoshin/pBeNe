#pragma once
// PrimaryGeneratorAction.hh
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UImessenger.hh"
#include "G4SystemOfUnits.hh"

class G4ParticleGun;
class G4Event;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIdirectory;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;
    void GeneratePrimaries(G4Event* event) override;

    // Setters called by messenger
    void SetMono(G4double e) { fEMono = e; }
    void SetTp  (G4double t) { fTp   = t; }
    void SetEmax(G4double e) { fEmax = e; }
    void SetMode(const G4String& m) { fMode = m; }

    // Getters for RunAction filename
    G4String GetMode()  const { return fMode; }
    G4double GetEMono() const { return fEMono; }
    G4double GetTp()    const { return fTp; }

private:
    class Messenger;
    G4ParticleGun* fParticleGun = nullptr;
    Messenger*     fMessenger   = nullptr;

    G4String fMode  = "mono";
    G4double fEMono = 9.0 * MeV;
    G4double fTp    = 1.5 * MeV;
    G4double fEmax  = 8.0 * MeV;

    G4double SampleTNSA() const;
};
