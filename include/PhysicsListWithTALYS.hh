//==============================================================================
// PhysicsListWithTALYS.hh  — added fOwner flag for MT-safe deletion
//==============================================================================
#pragma once
#include "G4VModularPhysicsList.hh"
#include "G4String.hh"

class MultiChannelSampler;

class PhysicsListWithTALYS : public G4VModularPhysicsList
{
public:
    explicit PhysicsListWithTALYS(const G4String& dataDir = "data");
    ~PhysicsListWithTALYS() override;

    void ConstructParticle() override;
    void ConstructProcess()  override;
    void SetCuts()           override;

private:
    void AddHandMadeProcesses();
    void AddNeutronProcesses();

    G4String             fDataDir;
    MultiChannelSampler* fSampler = nullptr;
    bool                 fOwner   = false;  // true only for master instance
};
