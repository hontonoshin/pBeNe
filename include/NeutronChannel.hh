#pragma once
//==============================================================================
// NeutronChannel.hh — enum + helpers for all reaction channels
//==============================================================================
#include "G4String.hh"
#include <array>
#include <string>

// ── Channel enumeration ───────────────────────────────────────────────────────
// Keep this in sync with MultiChannelSampler::LoadAll() cfg vector order.
enum class NeutronChannel : int {
    B11_pn      = 0,   // 11B(p,n)11C       TENDL-2025
    B10_pn      = 1,   // 10B(p,n)10C       TENDL-2025
    N14_pn      = 2,   // 14N(p,n)14O       TENDL-2025
    B11_alpha_n = 3,   // 11B(α,n)14N       JENDL-5
    B10_alpha_n = 4,   // 10B(α,n)13N       JENDL-5
    N14_alpha_n = 5,   // 14N(α,n)17F       JENDL-5
    B11_p2alpha = 6,   // 11B(p,2α)         TENDL-2025  [no neutron]
    OTHER       = 7,   // unidentified
};

static constexpr int kNChannels = 7;   // excludes OTHER

// ── ANSI colour codes for terminal output ─────────────────────────────────────
inline const char* kReset = "\033[0m";
inline G4String ChannelColour(NeutronChannel ch)
{
    switch (ch) {
        case NeutronChannel::B11_pn:      return "\033[1;32m"; // bright green
        case NeutronChannel::B10_pn:      return "\033[1;36m"; // bright cyan
        case NeutronChannel::N14_pn:      return "\033[1;34m"; // bright blue
        case NeutronChannel::B11_alpha_n: return "\033[1;33m"; // bright yellow
        case NeutronChannel::B10_alpha_n: return "\033[1;35m"; // bright magenta
        case NeutronChannel::N14_alpha_n: return "\033[0;33m"; // yellow
        case NeutronChannel::B11_p2alpha: return "\033[1;31m"; // bright red
        default:                          return "\033[0;37m"; // grey
    }
}

// ── Human-readable name ───────────────────────────────────────────────────────
inline G4String ChannelName(NeutronChannel ch)
{
    switch (ch) {
        case NeutronChannel::B11_pn:      return "11B(p,n)11C";
        case NeutronChannel::B10_pn:      return "10B(p,n)10C";
        case NeutronChannel::N14_pn:      return "14N(p,n)14O";
        case NeutronChannel::B11_alpha_n: return "11B(a,n)14N";
        case NeutronChannel::B10_alpha_n: return "10B(a,n)13N";
        case NeutronChannel::N14_alpha_n: return "14N(a,n)17F";
        case NeutronChannel::B11_p2alpha: return "11B(p,2a)";
        default:                          return "OTHER";
    }
}

// ── Short tag used in histogram / file names ──────────────────────────────────
inline G4String ChannelTag(NeutronChannel ch)
{
    switch (ch) {
        case NeutronChannel::B11_pn:      return "B11pn";
        case NeutronChannel::B10_pn:      return "B10pn";
        case NeutronChannel::N14_pn:      return "N14pn";
        case NeutronChannel::B11_alpha_n: return "B11an";
        case NeutronChannel::B10_alpha_n: return "B10an";
        case NeutronChannel::N14_alpha_n: return "N14an";
        case NeutronChannel::B11_p2alpha: return "B11pa";
        default:                          return "other";
    }
}

// ── Does this channel produce a neutron? ──────────────────────────────────────
inline bool ChannelProducesNeutron(NeutronChannel ch)
{
    return ch != NeutronChannel::B11_p2alpha && ch != NeutronChannel::OTHER;
}
