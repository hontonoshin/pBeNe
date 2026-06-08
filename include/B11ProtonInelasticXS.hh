#pragma once
//==============================================================================
// B11ProtonInelasticXS.hh
// External TALYS/TENDL 11B(p,n)11C cross section data set.
// Only xsMode="pn" is supported in this strict-external-data build.
//==============================================================================
#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"
#include "G4Proton.hh"
#include "G4SystemOfUnits.hh"
#include <vector>

class B11ProtonInelasticXS : public G4VCrossSectionDataSet {
public:
    explicit B11ProtonInelasticXS(const G4String& dataFile = "",
                                  const G4String& xsMode  = "pn");

    // G4VCrossSectionDataSet interface
    void BuildPhysicsTable(const G4ParticleDefinition& part) override;
    void DumpPhysicsTable (const G4ParticleDefinition& part) override;

    G4bool IsElementApplicable(const G4DynamicParticle*, G4int,
                               const G4Material*) override;
    G4bool IsIsoApplicable    (const G4DynamicParticle*, G4int, G4int,
                               const G4Element*, const G4Material*) override;

    G4double GetIsoCrossSection(const G4DynamicParticle*, G4int Z, G4int A,
                                const G4Isotope*, const G4Element*,
                                const G4Material*) override;
    G4double GetElementCrossSection(const G4DynamicParticle*, G4int Z,
                                    const G4Material*) override;

    // Public helper — returns cross section in mb for a given lab KE
    G4double GetCrossSection_mb(G4double Ep_MeV) const;

private:
    G4String fDataFile;
    G4String fXSMode;
    bool     fDataLoaded = false;

    static constexpr G4double fEthreshold = 3.02 * MeV;
    static constexpr G4double fEmax       = 200.0 * MeV;

    std::vector<G4double> fEp_pn, fXS_pn;
    std::vector<G4double> fEp_np, fXS_np;
    std::vector<G4double> fEp_na, fXS_na;
    std::vector<G4double> fEp_2n, fXS_2n;

    void LoadData(const G4String& fileBase);
    void LoadChannel(const G4String& fileName,
                     std::vector<G4double>& E,
                     std::vector<G4double>& xs);
    G4double Interpolate(const std::vector<G4double>& E,
                         const std::vector<G4double>& xs,
                         G4double Ep_MeV) const;
};
