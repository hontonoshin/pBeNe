//==============================================================================
// PrimaryGeneratorAction.cc — fixed: fills h1_proton_energy histogram
//==============================================================================
#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4AnalysisManager.hh"
#include "Randomize.hh"
#include <cmath>

class PrimaryGeneratorAction::Messenger : public G4UImessenger
{
public:
    explicit Messenger(PrimaryGeneratorAction* g) : fGen(g)
    {
        fDir = new G4UIdirectory("/beam/");
        fDir->SetGuidance("Proton beam settings");
        auto cmd = [&](const char* name, const char* hint) {
            auto* c = new G4UIcmdWithADoubleAndUnit(name, this);
            c->SetGuidance(hint);
            c->SetParameterName("val", false);
            c->SetRange("val>0");
            c->SetDefaultUnit("MeV");
            c->SetUnitCandidates("keV MeV");
            c->AvailableForStates(G4State_PreInit, G4State_Idle);
            return c;
        };
        fEnergyCmd = cmd("/beam/energy", "Mono proton energy");
        fTpCmd     = cmd("/beam/Tp",     "TNSA temperature T_p");
        fEmaxCmd   = cmd("/beam/Emax",   "TNSA high-energy cutoff");
        fModeCmd = new G4UIcmdWithAString("/beam/mode", this);
        fModeCmd->SetGuidance("Beam mode: mono | spectrum");
        fModeCmd->SetParameterName("mode", false);
        fModeCmd->SetCandidates("mono spectrum");
        fModeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    }
    ~Messenger() {
        delete fEnergyCmd; delete fTpCmd;
        delete fEmaxCmd;   delete fModeCmd; delete fDir;
    }
    void SetNewValue(G4UIcommand* cmd, G4String val) override {
        if      (cmd == fEnergyCmd) { fGen->SetMono(fEnergyCmd->GetNewDoubleValue(val)); fGen->SetMode("mono"); }
        else if (cmd == fTpCmd)       fGen->SetTp  (fTpCmd->GetNewDoubleValue(val));
        else if (cmd == fEmaxCmd)     fGen->SetEmax(fEmaxCmd->GetNewDoubleValue(val));
        else if (cmd == fModeCmd)     fGen->SetMode(val);
    }
private:
    PrimaryGeneratorAction*    fGen;
    G4UIdirectory*             fDir;
    G4UIcmdWithADoubleAndUnit* fEnergyCmd;
    G4UIcmdWithADoubleAndUnit* fTpCmd;
    G4UIcmdWithADoubleAndUnit* fEmaxCmd;
    G4UIcmdWithAString*        fModeCmd;
};

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction()
{
    fParticleGun = new G4ParticleGun(1);
    auto* proton = G4ParticleTable::GetParticleTable()->FindParticle("proton");
    fParticleGun->SetParticleDefinition(proton);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
    fMessenger = new Messenger(this);
    G4cout << "[PrimaryGeneratorAction] default: mono " << fEMono/MeV << " MeV\n";
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{ delete fMessenger; delete fParticleGun; }

G4double PrimaryGeneratorAction::SampleTNSA() const
{
    const G4double norm = 1.0 - std::exp(-fEmax / fTp);
    const G4double u    = G4UniformRand();
    const G4double E    = -fTp * std::log(1.0 - u * norm);
    return std::min(E, fEmax);
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
    const G4double energy = (fMode == "spectrum") ? SampleTNSA() : fEMono;

    fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., -0.76*mm));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
    fParticleGun->SetParticleEnergy(energy);
    fParticleGun->GeneratePrimaryVertex(event);

    // FIX: fill proton energy histogram (id=4, booked in RunAction)
    G4AnalysisManager::Instance()->FillH1(4, energy / MeV);
}
