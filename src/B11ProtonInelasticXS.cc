//==============================================================================
// B11ProtonInelasticXS.cc — external TALYS/TENDL 11B(p,n)11C XS only
//==============================================================================
#include "B11ProtonInelasticXS.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Exception.hh"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

//──────────────────────────────────────────────────────────────────────────────
B11ProtonInelasticXS::B11ProtonInelasticXS(const G4String& dataFile,
                                           const G4String& xsMode)
: G4VCrossSectionDataSet("B11_ProtonInelasticXS_TENDL2025"),
  fDataFile(dataFile),
  fXSMode(xsMode)
{
    SetMinKinEnergy(fEthreshold);
    SetMaxKinEnergy(fEmax);
}

//──────────────────────────────────────────────────────────────────────────────
void B11ProtonInelasticXS::LoadChannel(const G4String& fileName,
                                       std::vector<G4double>& E,
                                       std::vector<G4double>& xs)
{
    E.clear(); xs.clear();
    std::ifstream in(fileName.c_str());
    if (!in.good()) return;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        double ep, x;
        if (iss >> ep >> x) { E.push_back(ep); xs.push_back(std::max(0.0,x)); }
    }
    if (!E.empty())
        G4cout << "[B11ProtonInelasticXS] Loaded " << E.size()
               << " pts from " << fileName << G4endl;
}

//──────────────────────────────────────────────────────────────────────────────
void B11ProtonInelasticXS::LoadData(const G4String& fileBase)
{
    // Strict TALYS/TENDL-file mode.
    // No embedded fallback tables are allowed, because the calculation must be
    // reproducible from the external TALYS/TENDL data files only.
    fEp_pn.clear();
    fXS_pn.clear();
    fEp_np.clear();
    fXS_np.clear();
    fEp_na.clear();
    fXS_na.clear();
    fEp_2n.clear();
    fXS_2n.clear();

    if (fXSMode != "pn") {
        G4Exception("B11ProtonInelasticXS::LoadData",
                    "B11_XS_MODE_001",
                    FatalException,
                    ("This TALYS-only configuration supports only xsMode='pn'. Requested mode: "
                     + fXSMode).c_str());
    }

    G4String xsFile = fileBase;
    if (xsFile.empty()) {
        xsFile = "data/B11_pn_xs.dat";
    }

    LoadChannel(xsFile, fEp_pn, fXS_pn);

    if (fEp_pn.empty() || fXS_pn.empty()) {
        G4Exception("B11ProtonInelasticXS::LoadData",
                    "B11_XS_FILE_001",
                    FatalException,
                    ("Could not load external TALYS/TENDL 11B(p,n)11C cross-section file: "
                     + xsFile).c_str());
    }

    fDataLoaded = true;

    G4cout << "[B11ProtonInelasticXS] Loaded external TALYS/TENDL file only.\n"
           << "  mode    : " << fXSMode << "\n"
           << "  file    : " << xsFile << "\n"
           << "  points  : " << fEp_pn.size() << "\n"
           << "  peak mb : "
           << *std::max_element(fXS_pn.begin(), fXS_pn.end()) << "\n"
           << G4endl;
}

//──────────────────────────────────────────────────────────────────────────────
G4double B11ProtonInelasticXS::Interpolate(const std::vector<G4double>& E,
                                            const std::vector<G4double>& xs,
                                            G4double Ep_MeV) const
{
    if (E.empty()) return 0.0;
    if (Ep_MeV <= E.front()) return 0.0;
    if (Ep_MeV >= E.back())  return xs.back();

    auto it  = std::upper_bound(E.begin(), E.end(), Ep_MeV);
    int  idx = std::max(1, std::min((int)(it-E.begin()), (int)E.size()-1));

    const double E0=E[idx-1], E1=E[idx], X0=xs[idx-1], X1=xs[idx];
    const double t = (Ep_MeV-E0)/(E1-E0);
    if (X0>0 && X1>0 && E0>3.5) {
        double lt = std::log(Ep_MeV/E0)/std::log(E1/E0);
        return std::exp(std::log(X0)+lt*(std::log(X1)-std::log(X0)));
    }
    return std::max(0.0, X0+t*(X1-X0));
}

//──────────────────────────────────────────────────────────────────────────────
G4double B11ProtonInelasticXS::GetCrossSection_mb(G4double Ep_MeV) const
{
    if (!fDataLoaded) {
        const_cast<B11ProtonInelasticXS*>(this)->LoadData(
            fDataFile.empty() ? "data/B11_pn_xs.dat" : fDataFile);
    }
    return Interpolate(fEp_pn, fXS_pn, Ep_MeV);
}

//──────────────────────────────────────────────────────────────────────────────
void B11ProtonInelasticXS::BuildPhysicsTable(const G4ParticleDefinition& part)
{
    if (part.GetParticleName() != "proton") return;
    if (!fDataLoaded) LoadData(fDataFile.empty() ? "data/" : fDataFile);
    G4cout << "\n[B11ProtonInelasticXS] External TALYS/TENDL 11B(p,n)11C XS ready.\n";
}

//──────────────────────────────────────────────────────────────────────────────
// KEY FIX: check Z==5 and A==11, not just proton definition
//──────────────────────────────────────────────────────────────────────────────
G4bool B11ProtonInelasticXS::IsElementApplicable(
    const G4DynamicParticle* dp, G4int, const G4Material*)
{
    // The TALYS process is attached globally to protons.
    // Therefore Geant4 may ask this dataset for materials outside the 11B target
    // such as world/vacuum materials.  Return true for protons so the data store
    // always has an applicable dataset; the actual cross section is returned as
    // zero for every element except 11B in GetElementCrossSection/GetIsoCrossSection.
    return dp && dp->GetDefinition() == G4Proton::ProtonDefinition();
}

G4bool B11ProtonInelasticXS::IsIsoApplicable(
    const G4DynamicParticle* dp, G4int, G4int,
    const G4Element*, const G4Material*)
{
    // Same logic as IsElementApplicable: applicable to proton queries globally,
    // but numerically non-zero only for Z=5, A=11.
    return dp && dp->GetDefinition() == G4Proton::ProtonDefinition();
}

//──────────────────────────────────────────────────────────────────────────────
G4double B11ProtonInelasticXS::GetIsoCrossSection(
    const G4DynamicParticle* dp, G4int Z, G4int A,
    const G4Isotope*, const G4Element*, const G4Material*)
{
    if (!dp || dp->GetDefinition() != G4Proton::ProtonDefinition()) return 0.0;
    if (Z != 5 || A != 11) return 0.0;
    if (!fDataLoaded) LoadData(fDataFile.empty() ? "data/" : fDataFile);
    return GetCrossSection_mb(dp->GetKineticEnergy()/MeV) * millibarn;
}

G4double B11ProtonInelasticXS::GetElementCrossSection(
    const G4DynamicParticle* dp, G4int Z, const G4Material* mat)
{
    if (!dp || dp->GetDefinition() != G4Proton::ProtonDefinition()) return 0.0;
    if (Z != 5) return 0.0;
    return GetIsoCrossSection(dp, 5, 11, nullptr, nullptr, mat);
}

//──────────────────────────────────────────────────────────────────────────────
void B11ProtonInelasticXS::DumpPhysicsTable(const G4ParticleDefinition&)
{
    G4cout << "\n=== B11ProtonInelasticXS (" << fXSMode << ") ===\n"
           << std::setw(10) << "Ep [MeV]"
           << std::setw(12) << "s(p,n)"
           << std::setw(12) << "s(p,np)"
           << std::setw(12) << "s(p,na)"
           << std::setw(12) << "s(p,2n)"
           << std::setw(12) << "TOTAL [mb]\n"
           << std::string(70,'-') << "\n";
    for (double Ep : {4.,5.,6.,7.,8.,9.,10.,12.,15.,20.,25.,30.,40.,50.,70.}) {
        double pn = Interpolate(fEp_pn,fXS_pn,Ep);
        double np = fXSMode=="all_n" ? Interpolate(fEp_np,fXS_np,Ep) : 0.;
        double na = fXSMode=="all_n" ? Interpolate(fEp_na,fXS_na,Ep) : 0.;
        double tn = fXSMode=="all_n" ? Interpolate(fEp_2n,fXS_2n,Ep) : 0.;
        G4cout << std::setw(10) << Ep
               << std::setw(12) << pn
               << std::setw(12) << np
               << std::setw(12) << na
               << std::setw(12) << tn
               << std::setw(12) << pn+np+na+tn << "\n";
    }
    G4cout << G4endl;
}
