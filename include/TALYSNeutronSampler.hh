#pragma once
//==============================================================================
// TALYSNeutronSampler.hh
// Loads TALYS/TENDL double-differential (DDX) neutron tables and
// importance-samples (En, μ, φ) triples for one channel.
//==============================================================================
#include "G4Types.hh"
#include "G4String.hh"
#include <map>
#include <vector>

// ── Return type ───────────────────────────────────────────────────────────────
struct TALYSNeutronSample {
    G4double En_MeV = 0.;
    G4double mu     = 0.;   // cos(θ) in beam frame
    G4double phi    = 0.;   // azimuthal [0, 2π]
    bool     valid  = false;
};

// ── Sampler ───────────────────────────────────────────────────────────────────
class TALYSNeutronSampler {
public:
    TALYSNeutronSampler();

    // ── Loading ──────────────────────────────────────────────────────────────
    /// Load a single DDX file: columns  En_MeV  mu  d²σ/dEn/dΩ [mb/MeV/sr]
    void LoadDDXFile(const G4String& fileName, G4double Ep_MeV);

    /// Scan a directory for DDX files matching "{filePrefix}*MeV.dat".
    /// If filePrefix is empty, falls back to legacy "B11_pn_ddx_".
    void LoadDDXDirectory(const G4String& dirPath,
                          const G4String& filePrefix = "");

    /// Load 1-D energy spectrum (angle-integrated): En_MeV  dσ/dEn [mb/MeV]
    void LoadSpectrumFile(const G4String& fileName, G4double Ep_MeV);
    void LoadSpectraDirectory(const G4String& dirPath);

    // ── Sampling ─────────────────────────────────────────────────────────────
    TALYSNeutronSample Sample(G4double Ep_MeV) const;
    G4double SampleEnergy(G4double Ep_MeV) const;
    G4double MeanEnergy  (G4double Ep_MeV) const;

    // ── Accessors ─────────────────────────────────────────────────────────────
    int  NDDX()     const { return (int)fDDXTables.size(); }
    int  NSpectra() const { return (int)fSpectra.size(); }
    bool HasData()  const { return NDDX() > 0 || NSpectra() > 0; }

    void PrintLibrary() const;

private:
    // ── DDX table ─────────────────────────────────────────────────────────────
    struct DDXPoint {
        G4double En, mu, w;
        G4double eLo, eHi, muLo, muHi;
    };
    struct DDXTable {
        G4double Ep_MeV      = 0.;
        G4double totalWeight = 0.;
        G4double meanEn      = 0.;
        std::vector<DDXPoint>  points;
        std::vector<G4double>  cdf;
    };

    // ── 1-D spectrum CDF ──────────────────────────────────────────────────────
    struct SpectrumCDF {
        G4double Ep_MeV  = 0.;
        G4double totalXS = 0.;
        G4double meanEn  = 0.;
        std::vector<G4double> En;
        std::vector<G4double> cdf;
    };

    std::map<G4double, DDXTable>    fDDXTables;
    std::map<G4double, SpectrumCDF> fSpectra;

    DDXTable    BuildDDXTable(G4double Ep_MeV,
                              const std::vector<G4double>& en,
                              const std::vector<G4double>& mu,
                              const std::vector<G4double>& ddx) const;

    SpectrumCDF BuildCDF(G4double Ep_MeV,
                         const std::vector<G4double>& en,
                         const std::vector<G4double>& xs) const;

    TALYSNeutronSample SampleFromDDX(const DDXTable& table) const;
    TALYSNeutronSample SampleDDXInterpolated(const DDXTable& lo,
                                             const DDXTable& hi,
                                             G4double alpha) const;
    G4double SampleFromCDF(const SpectrumCDF& spec) const;
    G4double SampleInterpolated(const SpectrumCDF& lo,
                                const SpectrumCDF& hi,
                                G4double alpha) const;
    TALYSNeutronSample FallbackTwoBodyApprox(G4double Ep_MeV) const;
};
