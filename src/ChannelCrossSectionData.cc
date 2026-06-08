#include "ChannelCrossSectionData.hh"
#include "MultiChannelSampler.hh"

#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4Material.hh"

#include <algorithm>

ChannelCrossSectionData::ChannelCrossSectionData(
    MultiChannelSampler* sampler,
    NeutronChannel       channel,
    G4int projZ,
    G4int targZ,
    G4int targA)
: G4VCrossSectionDataSet("XS_" + ChannelTag(channel)),
  fSampler(sampler),
  fChannel(channel),
  fProjZ(projZ),
  fTargZ(targZ),
  fTargA(targA)
{}

G4bool ChannelCrossSectionData::IsCorrectProjectile(
    const G4DynamicParticle* dp) const
{
    if (!dp) return false;

    auto* def = dp->GetDefinition();
    if (!def) return false;

    if (fProjZ == 1 && def == G4Proton::ProtonDefinition()) return true;
    if (fProjZ == 2 && def == G4Alpha::AlphaDefinition())  return true;

    return false;
}

G4bool ChannelCrossSectionData::IsIsoApplicable(
    const G4DynamicParticle* dp,
    G4int,
    G4int,
    const G4Element*,
    const G4Material*)
{
    return IsCorrectProjectile(dp);
}

G4bool ChannelCrossSectionData::IsElementApplicable(
    const G4DynamicParticle* dp,
    G4int,
    const G4Material*)
{
    return IsCorrectProjectile(dp);
}

G4double ChannelCrossSectionData::GetIsoCrossSection(
    const G4DynamicParticle* dp,
    G4int Z,
    G4int A,
    const G4Isotope*,
    const G4Element*,
    const G4Material*)
{
    if (!IsCorrectProjectile(dp)) return 0.0;

    if (Z != fTargZ || A != fTargA) {
        return 0.0;
    }

    const G4double Einc_MeV = dp->GetKineticEnergy() / MeV;
    const G4double xs_mb = fSampler->GetXS_mb(fChannel, Einc_MeV);

    return xs_mb * millibarn;
}

G4double ChannelCrossSectionData::GetElementCrossSection(
    const G4DynamicParticle* dp,
    G4int Z,
    const G4Material* mat)
{
    if (!IsCorrectProjectile(dp)) return 0.0;

    if (Z != fTargZ) {
        return 0.0;
    }

    const G4double Einc_MeV = dp->GetKineticEnergy() / MeV;
    const G4double xs_mb = fSampler->GetXS_mb(fChannel, Einc_MeV);

    if (!mat) {
        return xs_mb * millibarn;
    }

    G4double isotopeFraction = 0.0;

    for (G4int i = 0; i < (G4int)mat->GetNumberOfElements(); ++i) {
        const G4Element* el = mat->GetElement(i);
        if (!el) continue;
        if ((G4int)el->GetZ() != fTargZ) continue;

        const G4double elementMassFraction = mat->GetFractionVector()[i];

        const G4int nIso = el->GetNumberOfIsotopes();

        if (nIso <= 0) {
            isotopeFraction += elementMassFraction;
            continue;
        }

        for (G4int j = 0; j < nIso; ++j) {
            const G4Isotope* iso = el->GetIsotope(j);
            if (!iso) continue;

            if (iso->GetZ() == fTargZ && iso->GetN() == fTargA) {
                isotopeFraction +=
                    elementMassFraction *
                    el->GetRelativeAbundanceVector()[j];
            }
        }
    }

    isotopeFraction = std::max(0.0, std::min(1.0, isotopeFraction));

    return isotopeFraction * xs_mb * millibarn;
}
