//==============================================================================
// MultiChannelSampler.cc
//
// FIX (2026-05-17): Missing XS data files for required channels now produce a
// G4Exception(FatalException) instead of silently registering an empty table
// and skipping the channel.  Previously a typo in a filename or a missing
// data/ subdirectory would cause B10(p,n) and all (α,n) channels to show zero
// yield across the entire scan with no error message.
//
// Channels that are marked "required" (see kRequiredChannels below) abort the
// run at LoadAll() time if their XS file cannot be read.  Optional channels
// (currently none, but easily extensible) still print a warning and continue.
//==============================================================================
#include "MultiChannelSampler.hh"

#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4Exception.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <set>

MultiChannelSampler::MultiChannelSampler(const G4String& dataDir)
: fDataDir(dataDir)
{}

//──────────────────────────────────────────────────────────────────────────────
void MultiChannelSampler::LoadAll()
{
    // Channel config: { enum, xs_filename, ddx_subdir, ddx_file_prefix, required }
    struct Cfg {
        NeutronChannel ch;
        std::string    xsFile;
        std::string    ddxDir;
        std::string    ddxPrefix;
        bool           required;   // if true, abort run when XS file is missing
    };

    const std::vector<Cfg> cfg = {
        { NeutronChannel::B11_pn,      "B11_pn_xs.dat",  "ddx",       "B11_pn_ddx_", true  },
        { NeutronChannel::B10_pn,      "B10_pn_xs.dat",  "ddx_B10pn", "B10_pn_ddx_", true  },
        { NeutronChannel::N14_pn,      "N14_pn_xs.dat",  "ddx_N14pn", "N14_pn_ddx_", true  },
        { NeutronChannel::B11_alpha_n, "B11_an_xs.dat",  "ddx_B11an", "B11_an_ddx_", true  },
        { NeutronChannel::B10_alpha_n, "B10_an_xs.dat",  "ddx_B10an", "B10_an_ddx_", true  },
        { NeutronChannel::N14_alpha_n, "N14_an_xs.dat",  "ddx_N14an", "N14_an_ddx_", true  },
        { NeutronChannel::B11_p2alpha, "B11_pa_xs.dat",  "ddx_B11pa", "B11_pa_ddx_", true  },
    };

    G4cout << "\n[MultiChannelSampler] Loading all channels from: "
           << fDataDir << "\n";

    for (const auto& c : cfg) {
        const int idx = static_cast<int>(c.ch);

        // XS file
        std::string xsPath = std::string(fDataDir) + "/" + c.xsFile;
        LoadXSFile(c.ch, xsPath, c.required);

        // DDX directory
        std::string ddxPath = std::string(fDataDir) + "/" + c.ddxDir;
        fSamplers[idx].LoadDDXDirectory(ddxPath, c.ddxPrefix);

        G4cout << "  " << std::left << std::setw(20) << ChannelName(c.ch)
               << "  XS pts=" << std::setw(4) << fXSTables[idx].size()
               << "  DDX tables=" << fSamplers[idx].NDDX()
               << "\n";
    }
    G4cout << "\n";
}

//──────────────────────────────────────────────────────────────────────────────
// FIX: added 'required' parameter.  When true and the file is absent/empty,
// throw a FatalException so the run aborts with a clear diagnostic instead of
// silently producing zero yield for that channel.
//──────────────────────────────────────────────────────────────────────────────
void MultiChannelSampler::LoadXSFile(NeutronChannel ch,
                                     const std::string& path,
                                     bool required)
{
    const int idx = static_cast<int>(ch);
    fXSTables[idx].clear();

    std::ifstream in(path);
    if (!in.good()) {
        if (required) {
            G4Exception("MultiChannelSampler::LoadXSFile",
                        "MULTICHAN_XS_MISSING",
                        FatalException,
                        ("Required XS file not found: " + path +
                         "\nCheck that the data/ directory is complete and "
                         "filenames match the channel configuration.").c_str());
        } else {
            G4cout << "  [MultiChannelSampler] Optional XS file not found (skipping): "
                   << path << "\n";
        }
        return;
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        double E, xs;
        if (iss >> E >> xs && E > 0 && xs >= 0)
            fXSTables[idx].push_back({E, xs});
    }

    if (fXSTables[idx].empty()) {
        if (required) {
            G4Exception("MultiChannelSampler::LoadXSFile",
                        "MULTICHAN_XS_EMPTY",
                        FatalException,
                        ("Required XS file is present but contains no valid data: "
                         + path).c_str());
        } else {
            G4cout << "  [MultiChannelSampler] Optional XS file empty (skipping): "
                   << path << "\n";
        }
    }
}

//──────────────────────────────────────────────────────────────────────────────
G4double MultiChannelSampler::InterpolateXS(NeutronChannel ch,
                                             G4double Ep_MeV) const
{
    const auto& table = fXSTables[static_cast<int>(ch)];
    if (table.empty()) return 0.;

    if (Ep_MeV <= table.front().E) return 0.;
    if (Ep_MeV >= table.back().E)  return table.back().xs;

    auto it = std::lower_bound(table.begin(), table.end(), Ep_MeV,
        [](const XSPoint& p, double e){ return p.E < e; });
    if (it == table.begin()) return 0.;
    const auto& hi = *it;
    const auto& lo = *std::prev(it);

    const double t = (Ep_MeV - lo.E) / (hi.E - lo.E);
    if (lo.xs > 0 && hi.xs > 0 && lo.E > 0) {
        double lt = std::log(Ep_MeV / lo.E) / std::log(hi.E / lo.E);
        return std::exp(std::log(lo.xs) + lt * (std::log(hi.xs) - std::log(lo.xs)));
    }
    return std::max(0., lo.xs + t * (hi.xs - lo.xs));
}

//──────────────────────────────────────────────────────────────────────────────
G4double MultiChannelSampler::GetXS_mb(NeutronChannel ch,
                                        G4double Ep_MeV) const
{
    return InterpolateXS(ch, Ep_MeV);
}

//──────────────────────────────────────────────────────────────────────────────
TALYSNeutronSample MultiChannelSampler::Sample(NeutronChannel ch,
                                                G4double Ep_MeV) const
{
    return fSamplers[static_cast<int>(ch)].Sample(Ep_MeV);
}

//──────────────────────────────────────────────────────────────────────────────
void MultiChannelSampler::PrintSummary() const
{
    G4cout << "\n=== MultiChannelSampler library ===\n"
           << std::left  << std::setw(22) << "Channel"
           << std::right << std::setw(10) << "XS pts"
           << std::setw(12) << "DDX tables"
           << std::setw(14) << "Peak XS [mb]"
           << "\n" << std::string(58,'-') << "\n";

    for (int i = 0; i < kNChannels; ++i) {
        auto ch = static_cast<NeutronChannel>(i);
        const auto& table = fXSTables[i];
        double peakXS = 0.;
        for (const auto& p : table) peakXS = std::max(peakXS, p.xs);

        G4cout << ChannelColour(ch)
               << std::left  << std::setw(22) << ChannelName(ch)
               << std::right
               << std::setw(10) << table.size()
               << std::setw(12) << fSamplers[i].NDDX()
               << std::setw(14) << std::fixed << std::setprecision(1) << peakXS
               << kReset << "\n";
    }
    G4cout << "\n";
}
