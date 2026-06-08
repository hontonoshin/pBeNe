#pragma once
// TALYSFinalStateModel.hh — legacy single-channel (B11 p,n) model kept for back-compat
#include "G4HadronicInteraction.hh"
#include "G4HadFinalState.hh"
#include "G4ThreeVector.hh"

class TALYSNeutronSampler;

class TALYSFinalStateModel : public G4HadronicInteraction {
public:
    explicit TALYSFinalStateModel(TALYSNeutronSampler* sampler);

    G4bool IsApplicable(const G4HadProjectile& projectile,
                        G4Nucleus& targetNucleus) override;

    G4HadFinalState* ApplyYourself(const G4HadProjectile& projectile,
                                   G4Nucleus& target) override;

private:
    TALYSNeutronSampler* fSampler;
    G4HadFinalState      fLocalState;

    G4ThreeVector DirectionFromMuPhi(const G4ThreeVector& axis,
                                     G4double mu, G4double phi) const;
};
