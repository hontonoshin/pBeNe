#include "TALYSNeutronSampler.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <dirent.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {
std::vector<double> UniqueSorted(std::vector<double> v)
{
    std::sort(v.begin(), v.end());
    std::vector<double> out;
    for (double x : v) {
        if (out.empty() || std::abs(x - out.back()) > 1e-10) out.push_back(x);
    }
    return out;
}

std::map<double, std::pair<double,double>> MakeBinEdges(const std::vector<double>& grid,
                                                        double lowClamp,
                                                        double highClamp)
{
    std::map<double, std::pair<double,double>> edges;
    if (grid.empty()) return edges;
    if (grid.size() == 1) {
        double half = 0.5;
        edges[grid[0]] = {std::max(lowClamp, grid[0]-half), std::min(highClamp, grid[0]+half)};
        return edges;
    }
    for (std::size_t i=0; i<grid.size(); ++i) {
        double lo = (i==0) ? grid[i] - 0.5*(grid[i+1]-grid[i]) : 0.5*(grid[i-1]+grid[i]);
        double hi = (i+1==grid.size()) ? grid[i] + 0.5*(grid[i]-grid[i-1]) : 0.5*(grid[i]+grid[i+1]);
        edges[grid[i]] = {std::max(lowClamp, lo), std::min(highClamp, hi)};
    }
    return edges;
}

bool HasPrefixSuffix(const std::string& name, const std::string& prefix, const std::string& suffix)
{
    return name.size() >= prefix.size()+suffix.size()
        && name.substr(0, prefix.size()) == prefix
        && name.substr(name.size()-suffix.size()) == suffix;
}

double ExtractEnergyFromName(const std::string& name,
                             const std::string& prefix,
                             const std::string& suffix)
{
    std::string e = name.substr(prefix.size(),
                                name.size() - prefix.size() - suffix.size());

    std::replace(e.begin(), e.end(), 'p', '.');

    return std::stod(e);
}
}

TALYSNeutronSampler::TALYSNeutronSampler() {}

void TALYSNeutronSampler::LoadDDXFile(const G4String& fileName, G4double Ep_MeV)
{
    std::ifstream in(fileName.c_str());
    if (!in.good()) {
        G4cerr << "[TALYSNeutronSampler] WARNING: cannot open DDX file '" << fileName << "'\n";
        return;
    }

    std::vector<G4double> en, mu, ddx;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        std::istringstream iss(line);
        double e, m, v;
        if (iss >> e >> m >> v) {
            if (e >= 0. && m >= -1.000001 && m <= 1.000001 && v > 0.) {
                en.push_back(e);
                mu.push_back(std::max(-1.0, std::min(1.0, m)));
                ddx.push_back(v);
            }
        }
    }

    if (en.size() < 2) {
        G4cerr << "[TALYSNeutronSampler] Too few DDX points in " << fileName << " — skipping.\n";
        return;
    }

    DDXTable table = BuildDDXTable(Ep_MeV, en, mu, ddx);
    if (table.cdf.empty() || table.totalWeight <= 0.) {
        G4cerr << "[TALYSNeutronSampler] Empty/zero DDX table in " << fileName << " — skipping.\n";
        return;
    }
    fDDXTables[Ep_MeV] = std::move(table);

    G4cout << "[TALYSNeutronSampler] Loaded DDX Ep=" << Ep_MeV
           << " MeV: " << fDDXTables[Ep_MeV].points.size() << " cells, "
           << "<En>=" << fDDXTables[Ep_MeV].meanEn << " MeV, "
           << "norm=" << fDDXTables[Ep_MeV].totalWeight << "\n";
}

//──────────────────────────────────────────────────────────────────────────────
// FIX: accept an explicit prefix so every channel scans its own subdir
// correctly. The original code hardcoded "B11_pn_ddx_" which broke all
// channels except B11(p,n).
//──────────────────────────────────────────────────────────────────────────────
void TALYSNeutronSampler::LoadDDXDirectory(const G4String& dirPath,
                                           const G4String& filePrefix)
{
    DIR* dir = opendir(dirPath.c_str());
    if (!dir) {
        G4cerr << "[TALYSNeutronSampler] Cannot open DDX directory: " << dirPath << G4endl;
        return;
    }

    // Use the caller-supplied prefix; fall back to legacy "B11_pn_ddx_"
    const std::string prefix = filePrefix.empty()
                               ? std::string("B11_pn_ddx_")
                               : std::string(filePrefix);
    const std::string suffix = "MeV.dat";
    int loaded = 0;
    struct dirent* ent;
    while ((ent = readdir(dir)) != nullptr) {
        std::string name(ent->d_name);
        if (!HasPrefixSuffix(name, prefix, suffix)) continue;
        double Ep = ExtractEnergyFromName(name, prefix, suffix);
        LoadDDXFile(std::string(dirPath) + "/" + name, Ep);
        ++loaded;
    }
    closedir(dir);
    G4cout << "[TALYSNeutronSampler] DDX scan (" << prefix << "*" << suffix
           << "): " << loaded << " file(s) in " << dirPath << G4endl;
}

TALYSNeutronSampler::DDXTable
TALYSNeutronSampler::BuildDDXTable(G4double Ep_MeV,
                                   const std::vector<G4double>& en,
                                   const std::vector<G4double>& mu,
                                   const std::vector<G4double>& ddx) const
{
    DDXTable table;
    table.Ep_MeV = Ep_MeV;

    auto eGrid  = UniqueSorted(std::vector<double>(en.begin(), en.end()));
    auto muGrid = UniqueSorted(std::vector<double>(mu.begin(), mu.end()));
    auto eEdges = MakeBinEdges(eGrid, 0.0, 1.0e9);
    auto mEdges = MakeBinEdges(muGrid, -1.0, 1.0);

    table.points.reserve(en.size());
    table.cdf.reserve(en.size());

    G4double total = 0.;
    G4double sumE  = 0.;

    for (std::size_t i=0; i<en.size(); ++i) {
        auto ee = eEdges[en[i]];
        auto me = mEdges[mu[i]];
        const double dE  = std::max(0.0, ee.second - ee.first);
        const double dmu = std::max(0.0, me.second - me.first);
        const double weight = ddx[i] * dE * (twopi * dmu);
        if (weight <= 0.) continue;

        DDXPoint p;
        p.En = en[i]; p.mu = mu[i]; p.w = weight;
        p.eLo = ee.first; p.eHi = ee.second;
        p.muLo = me.first; p.muHi = me.second;
        table.points.push_back(p);
        total += weight;
        sumE  += weight * en[i];
    }

    table.totalWeight = total;
    table.meanEn = (total > 0.) ? sumE/total : 0.;

    G4double cumul = 0.;
    for (const auto& p : table.points) {
        cumul += p.w / total;
        table.cdf.push_back(cumul);
    }
    if (!table.cdf.empty()) table.cdf.back() = 1.0;
    return table;
}

TALYSNeutronSample TALYSNeutronSampler::SampleFromDDX(const DDXTable& table) const
{
    TALYSNeutronSample s;
    if (table.points.empty() || table.cdf.empty()) return s;

    double r = G4UniformRand();
    auto it = std::lower_bound(table.cdf.begin(), table.cdf.end(), r);
    int idx = (int)(it - table.cdf.begin());
    idx = std::max(0, std::min(idx, (int)table.points.size()-1));
    const auto& p = table.points[idx];

    s.En_MeV = p.eLo + G4UniformRand() * std::max(0.0, p.eHi - p.eLo);
    s.mu     = p.muLo + G4UniformRand() * std::max(0.0, p.muHi - p.muLo);
    s.mu     = std::max(-1.0, std::min(1.0, s.mu));
    s.phi    = twopi * G4UniformRand();
    s.valid  = true;
    return s;
}

TALYSNeutronSample TALYSNeutronSampler::SampleDDXInterpolated(const DDXTable& lo,
                                                              const DDXTable& hi,
                                                              G4double alpha) const
{
    return (G4UniformRand() < alpha) ? SampleFromDDX(hi) : SampleFromDDX(lo);
}

void TALYSNeutronSampler::LoadSpectrumFile(const G4String& fileName, G4double Ep_MeV)
{
    std::ifstream in(fileName.c_str());
    if (!in.good()) return;
    std::vector<G4double> enVec, xsVec;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0]=='#') continue;
        std::istringstream iss(line);
        double en, xs;
        if (iss >> en >> xs && xs >= 0.) {
            enVec.push_back(en);
            xsVec.push_back(xs);
        }
    }
    if (enVec.size() < 2) return;
    SpectrumCDF cdf = BuildCDF(Ep_MeV, enVec, xsVec);
    fSpectra[Ep_MeV] = std::move(cdf);
    G4cout << "[TALYSNeutronSampler] Loaded 1D spectrum Ep=" << Ep_MeV
           << " MeV: " << enVec.size() << " points, <En>="
           << fSpectra[Ep_MeV].meanEn << " MeV\n";
}

void TALYSNeutronSampler::LoadSpectraDirectory(const G4String& dirPath)
{
    DIR* dir = opendir(dirPath.c_str());
    if (!dir) return;
    const std::string prefix = "B11_pn_nspec_";
    const std::string suffix = "MeV.dat";
    int loaded = 0;
    struct dirent* ent;
    while ((ent = readdir(dir)) != nullptr) {
        std::string name(ent->d_name);
        if (!HasPrefixSuffix(name, prefix, suffix)) continue;
        double Ep = ExtractEnergyFromName(name, prefix, suffix);
        LoadSpectrumFile(std::string(dirPath) + "/" + name, Ep);
        ++loaded;
    }
    closedir(dir);
    if (loaded > 0) G4cout << "[TALYSNeutronSampler] 1D spectrum scan: " << loaded << " file(s) matched.\n";
}

TALYSNeutronSampler::SpectrumCDF
TALYSNeutronSampler::BuildCDF(G4double Ep_MeV,
                              const std::vector<G4double>& en,
                              const std::vector<G4double>& xs) const
{
    SpectrumCDF cdf;
    cdf.Ep_MeV = Ep_MeV;
    const int N = (int)en.size();
    if (N < 2) return cdf;

    std::vector<G4double> dI(N-1, 0.);
    G4double total = 0.;
    G4double sumE = 0.;
    for (int i=0; i<N-1; ++i) {
        double dE = en[i+1] - en[i];
        if (dE <= 0.) continue;
        double area = 0.5*(xs[i] + xs[i+1]) * dE;
        dI[i] = std::max(0.0, area);
        total += dI[i];
        sumE += dI[i] * 0.5*(en[i] + en[i+1]);
    }
    cdf.totalXS = total;
    cdf.meanEn = (total > 0.) ? sumE/total : 0.;
    cdf.En.push_back(en[0]);
    cdf.cdf.push_back(0.);
    G4double cumul = 0.;
    for (int i=0; i<N-1; ++i) {
        cumul += (total > 0.) ? dI[i]/total : 0.;
        cdf.En.push_back(en[i+1]);
        cdf.cdf.push_back(cumul);
    }
    if (!cdf.cdf.empty()) cdf.cdf.back() = 1.;
    return cdf;
}

G4double TALYSNeutronSampler::SampleFromCDF(const SpectrumCDF& spec) const
{
    if (spec.cdf.size() < 2) return std::max(1.e-3, spec.meanEn);
    G4double r = G4UniformRand();
    auto it = std::lower_bound(spec.cdf.begin(), spec.cdf.end(), r);
    int idx = (int)(it - spec.cdf.begin());
    idx = std::max(1, std::min(idx, (int)spec.En.size()-1));
    double c0 = spec.cdf[idx-1], c1 = spec.cdf[idx];
    double e0 = spec.En[idx-1],  e1 = spec.En[idx];
    if (c1 <= c0) return e0;
    double t = (r-c0)/(c1-c0);
    return e0 + t*(e1-e0);
}

G4double TALYSNeutronSampler::SampleInterpolated(const SpectrumCDF& lo,
                                                 const SpectrumCDF& hi,
                                                 G4double alpha) const
{
    return (G4UniformRand() < alpha) ? SampleFromCDF(hi) : SampleFromCDF(lo);
}

TALYSNeutronSample TALYSNeutronSampler::FallbackTwoBodyApprox(G4double Ep_MeV) const
{
    TALYSNeutronSample s;
    const double Q = -2.764036;
    double available = std::max(1.e-3, Ep_MeV + Q);
    s.En_MeV = std::max(1.e-3, 0.80 * available);
    s.mu = 2.0*G4UniformRand() - 1.0;
    s.phi = twopi * G4UniformRand();
    s.valid = true;
    static bool warned = false;
    if (!warned) {
        G4cerr << "[TALYSNeutronSampler] WARNING: no DDX/spectrum tables loaded. "
               << "Using rough two-body placeholder. Do not use for final physics.\n";
        warned = true;
    }
    return s;
}

TALYSNeutronSample TALYSNeutronSampler::Sample(G4double Ep_MeV) const
{
    if (!fDDXTables.empty()) {
        auto hi = fDDXTables.lower_bound(Ep_MeV);
        if (hi == fDDXTables.begin()) return SampleFromDDX(hi->second);
        if (hi == fDDXTables.end()) return SampleFromDDX(std::prev(hi)->second);
        auto lo = std::prev(hi);
        double alpha = (Ep_MeV - lo->first) / (hi->first - lo->first);
        alpha = std::max(0.0, std::min(1.0, alpha));
        return SampleDDXInterpolated(lo->second, hi->second, alpha);
    }

    if (!fSpectra.empty()) {
        auto hi = fSpectra.lower_bound(Ep_MeV);
        G4double En = 0.;
        if (hi == fSpectra.begin()) En = SampleFromCDF(hi->second);
        else if (hi == fSpectra.end()) En = SampleFromCDF(std::prev(hi)->second);
        else {
            auto lo = std::prev(hi);
            double alpha = (Ep_MeV - lo->first) / (hi->first - lo->first);
            En = SampleInterpolated(lo->second, hi->second, alpha);
        }
        TALYSNeutronSample s;
        s.En_MeV = En;
        s.mu = 2.0*G4UniformRand() - 1.0;
        s.phi = twopi * G4UniformRand();
        s.valid = true;
        return s;
    }

    return FallbackTwoBodyApprox(Ep_MeV);
}

G4double TALYSNeutronSampler::SampleEnergy(G4double Ep_MeV) const
{
    return Sample(Ep_MeV).En_MeV;
}

G4double TALYSNeutronSampler::MeanEnergy(G4double Ep_MeV) const
{
    if (!fDDXTables.empty()) {
        auto hi = fDDXTables.lower_bound(Ep_MeV);
        if (hi == fDDXTables.begin()) return hi->second.meanEn;
        if (hi == fDDXTables.end()) return std::prev(hi)->second.meanEn;
        auto lo = std::prev(hi);
        double alpha = (Ep_MeV - lo->first) / (hi->first - lo->first);
        return (1.0-alpha)*lo->second.meanEn + alpha*hi->second.meanEn;
    }
    if (!fSpectra.empty()) {
        auto hi = fSpectra.lower_bound(Ep_MeV);
        if (hi == fSpectra.begin()) return hi->second.meanEn;
        if (hi == fSpectra.end()) return std::prev(hi)->second.meanEn;
        auto lo = std::prev(hi);
        double alpha = (Ep_MeV - lo->first) / (hi->first - lo->first);
        return (1.0-alpha)*lo->second.meanEn + alpha*hi->second.meanEn;
    }
    return 0.;
}

void TALYSNeutronSampler::PrintLibrary() const
{
    G4cout << "\n=== TALYSNeutronSampler DDX library ===\n";
    G4cout << std::setw(12) << "Ep [MeV]" << std::setw(14) << "<En> [MeV]"
           << std::setw(16) << "Norm [mb]" << std::setw(12) << "Cells" << "\n";
    G4cout << std::string(54, '-') << "\n";
    for (const auto& kv : fDDXTables) {
        G4cout << std::setw(12) << kv.first
               << std::setw(14) << kv.second.meanEn
               << std::setw(16) << kv.second.totalWeight
               << std::setw(12) << kv.second.points.size() << "\n";
    }
    if (fDDXTables.empty()) G4cout << "  No DDX tables loaded.\n";

    if (!fSpectra.empty()) {
        G4cout << "\n=== 1D fallback spectra ===\n";
        for (const auto& kv : fSpectra) {
            G4cout << "  Ep=" << kv.first << " MeV  <En>=" << kv.second.meanEn
                   << " MeV  bins=" << kv.second.En.size() << "\n";
        }
    }
    G4cout << G4endl;
}
