#pragma once
//==============================================================================
// ChannelCrossSectionData.hh
// G4VCrossSectionDataSet wrapper — returns the TENDL/JENDL XS for one channel.
//==============================================================================
#include "G4VCrossSectionDataSet.hh"
#include "NeutronChannel.hh"

class MultiChannelSampler;

class ChannelCrossSectionData : public G4VCrossSectionDataSet {
public:
    ChannelCrossSectionData(MultiChannelSampler* sampler,
                            NeutronChannel       channel,
                            G4int projZ,
                            G4int targZ, G4int targA);

    G4bool IsIsoApplicable(const G4DynamicParticle*, G4int Z, G4int A,
                           const G4Element*, const G4Material*) override;

    G4bool IsElementApplicable(const G4DynamicParticle*, G4int Z,
                               const G4Material*) override;

    G4double GetIsoCrossSection(const G4DynamicParticle*, G4int Z, G4int A,
                                const G4Isotope*, const G4Element*,
                                const G4Material*) override;

    G4double GetElementCrossSection(const G4DynamicParticle*, G4int Z,
                                    const G4Material*) override;

private:
    MultiChannelSampler* fSampler;
    NeutronChannel       fChannel;
    G4int fProjZ, fTargZ, fTargA;
    
    
    
private:
    G4bool IsCorrectProjectile(const G4DynamicParticle* dp) const;
};


