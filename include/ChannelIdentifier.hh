#pragma once
//==============================================================================
// ChannelIdentifier.hh
// Maps (processName, projectilePDG, targetZ, targetA) → NeutronChannel.
// Used in SteppingAction to classify G4-native secondaries.
//==============================================================================
#include "NeutronChannel.hh"
#include "G4String.hh"
#include <string>

inline NeutronChannel IdentifyChannel(const G4String& procName,
                                      G4int projPDG,
                                      G4int tZ, G4int tA)
{
    const bool isProton = (projPDG == 2212);
    const bool isAlpha  = (projPDG == 1000020040);

    // Hand-made TALYS process names (registered in PhysicsListWithTALYS)
    if (procName == "B11_pn_TALYS") return NeutronChannel::B11_pn;
    if (procName == "B10_pn_TALYS") return NeutronChannel::B10_pn;
    if (procName == "N14_pn_TALYS") return NeutronChannel::N14_pn;
    if (procName == "B11_an_TALYS") return NeutronChannel::B11_alpha_n;
    if (procName == "B10_an_TALYS") return NeutronChannel::B10_alpha_n;
    if (procName == "N14_an_TALYS") return NeutronChannel::N14_alpha_n;
    if (procName == "B11_pa_TALYS") return NeutronChannel::B11_p2alpha;

    // Fallback: classify by target nucleus
    if (isProton) {
        if (tZ == 5 && tA == 11) return NeutronChannel::B11_pn;
        if (tZ == 5 && tA == 10) return NeutronChannel::B10_pn;
        if (tZ == 7 && tA == 14) return NeutronChannel::N14_pn;
    }
    if (isAlpha) {
        if (tZ == 5 && tA == 11) return NeutronChannel::B11_alpha_n;
        if (tZ == 5 && tA == 10) return NeutronChannel::B10_alpha_n;
        if (tZ == 7 && tA == 14) return NeutronChannel::N14_alpha_n;
    }
    return NeutronChannel::OTHER;
}
