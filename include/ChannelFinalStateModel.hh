#pragma once
//==============================================================================
// ChannelFinalStateModel.hh
// G4HadronicInteraction that samples (En, μ, φ) from TALYS DDX tables
// and constructs the two-body final state with momentum conservation.
//==============================================================================
#include "G4HadronicInteraction.hh"
#include "G4HadFinalState.hh"
#include "G4ThreeVector.hh"
#include "NeutronChannel.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleDefinition.hh"

class MultiChannelSampler;

class ChannelFinalStateModel : public G4HadronicInteraction {
public:
    ChannelFinalStateModel(MultiChannelSampler* sampler,
                           NeutronChannel       channel,
                           G4int projZ,  G4int projA,
                           G4int targZ,  G4int targA,
                           G4int ejectZ, G4int ejectA,
                           G4int recoilZ,G4int recoilA,
                           G4double Qvalue_MeV);

    G4bool IsApplicable(const G4HadProjectile& proj,
                        G4Nucleus& target) override;

    G4HadFinalState* ApplyYourself(const G4HadProjectile& projectile,
                                   G4Nucleus& target) override;

private:
    MultiChannelSampler* fSampler;
    NeutronChannel       fChannel;

    G4int    fProjZ,   fProjA;
    G4int    fTargZ,   fTargA;
    G4int    fEjectZ,  fEjectA;
    G4int    fRecoilZ, fRecoilA;
    G4double fQvalue;

    G4HadFinalState fLocalState;

    G4ThreeVector DirectionFromMuPhi(const G4ThreeVector& axis,
                                     G4double mu, G4double phi) const;
                                     


private:
    G4ParticleDefinition* GetParticleDef(
        G4int Z,
        G4int A,
        G4double excitation = 0.0
    ) const;

    G4bool AddTwoBodyFinalStateCM(
        const G4LorentzVector& p4tot,
        const G4ThreeVector& projDir,
        G4ParticleDefinition* ejectDef,
        G4ParticleDefinition* recoilDef
    );

};
