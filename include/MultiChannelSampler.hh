//==============================================================================
// MultiChannelSampler.hh  — thread-safe singleton, loaded once on master
//==============================================================================
#pragma once

#include "NeutronChannel.hh"
#include "TALYSNeutronSampler.hh"
#include "G4String.hh"
#include "G4ios.hh"
#include <array>
#include <vector>
#include <memory>

class MultiChannelSampler
{
public:
    explicit MultiChannelSampler(const G4String& dataDir = "data");

    // Called once on master thread
    void LoadAll();
    void PrintSummary() const;

    // Thread-safe read-only accessors (called from worker threads)
    TALYSNeutronSample Sample(NeutronChannel ch, G4double Ep_MeV) const;
    G4double GetXS_mb(NeutronChannel ch, G4double Ep_MeV) const;

    bool HasXS (NeutronChannel ch) const {
        return !fXSTables[static_cast<int>(ch)].empty();
    }
    bool HasDDX(NeutronChannel ch) const {
        return fSamplers[static_cast<int>(ch)].NDDX() > 0;
    }
    bool IsLoaded(NeutronChannel ch) const {
        return HasXS(ch) || HasDDX(ch);
    }
    int NDDX(NeutronChannel ch) const {
        return fSamplers[static_cast<int>(ch)].NDDX();
    }

private:
    G4double InterpolateXS(NeutronChannel ch, G4double Ep_MeV) const;
    void LoadXSFile(NeutronChannel ch, const std::string& path,
                bool required = false);

    G4String fDataDir;

    // Use unique_ptr array to avoid copy issues with MT
    // TALYSNeutronSampler is not copyable safely across threads
    std::array<TALYSNeutronSampler, kNChannels> fSamplers;

    struct XSPoint { G4double E, xs; };
    std::array<std::vector<XSPoint>, kNChannels> fXSTables;
};
