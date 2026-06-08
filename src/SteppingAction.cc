//==============================================================================
// SteppingAction.cc — multi-channel neutron and alpha scoring
//
// FIXES (2026-05-17):
//   1. Alpha scoring now uses EXACT process name "B11_pa_TALYS" only.
//      The previous broad proc.find("Inelastic") check matched every single
//      hand-made G4HadronInelasticProcess, causing alpha recoils from B11(p,n),
//      B10(p,n), N14(p,n), and all (α,n) channels to be mis-attributed to
//      B11_p2alpha — inflating that counter and corrupting the summary table.
//      Neutron recording is completely unchanged.
//   2. IsHandMadeProtonProcess() helper removed — it was only used for the
//      broken broad-match alpha scoring and is no longer needed.
//==============================================================================
#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "ChannelIdentifier.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4VProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4DynamicParticle.hh"
#include "G4Alpha.hh"

SteppingAction::SteppingAction(EventAction* ea, const DetectorConstruction* dc)
: G4UserSteppingAction(),
  fEventAction(ea),
  fDetector(dc)
{
    fChannelNeutrons.fill(0);
    fChannelAlphas.fill(0);
}

SteppingAction::~SteppingAction() {}

//──────────────────────────────────────────────────────────────────────────────
static std::pair<G4int,G4int> InferTargetZA(const G4Material* mat,
                                             G4int /*projectilePDG*/)
{
    if (!mat) return {0,0};
    G4int bestZ = 0, bestA = 0;
    G4double bestFrac = 0.;
    for (G4int ei = 0; ei < (G4int)mat->GetNumberOfElements(); ++ei) {
        const G4Element* el = mat->GetElement(ei);
        const G4double frac = mat->GetFractionVector()[ei];
        for (G4int ii = 0; ii < (G4int)el->GetNumberOfIsotopes(); ++ii) {
            const G4Isotope* iso = el->GetIsotope(ii);
            if (iso->GetZ() <= 1) continue;
            const G4double isoFrac = frac * el->GetRelativeAbundanceVector()[ii];
            if (isoFrac > bestFrac) {
                bestFrac = isoFrac;
                bestZ    = iso->GetZ();
                bestA    = iso->GetN();
            }
        }
    }
    return {bestZ, bestA};
}

//──────────────────────────────────────────────────────────────────────────────
NeutronChannel SteppingAction::ClassifySecondary(const G4Track* track,
                                                  const G4Step*  step) const
{
    const G4VProcess* creator = track->GetCreatorProcess();
    if (!creator) return NeutronChannel::OTHER;
    const G4String& procName = creator->GetProcessName();

    // All hand-made TALYS processes → classify by exact name
    if (procName == "B11_pn_TALYS") return NeutronChannel::B11_pn;
    if (procName == "B10_pn_TALYS") return NeutronChannel::B10_pn;
    if (procName == "N14_pn_TALYS") return NeutronChannel::N14_pn;
    if (procName == "B11_an_TALYS") return NeutronChannel::B11_alpha_n;
    if (procName == "B10_an_TALYS") return NeutronChannel::B10_alpha_n;
    if (procName == "N14_an_TALYS") return NeutronChannel::N14_alpha_n;
    if (procName == "B11_pa_TALYS") return NeutronChannel::B11_p2alpha;

    // Fallback: infer from material
    const G4Material* mat = step->GetPreStepPoint()->GetMaterial();
    G4int projPDG = 2212;
    if (procName.find("alpha") != std::string::npos ||
        procName.find("Alpha") != std::string::npos)
        projPDG = 1000020040;
    auto [tZ, tA] = InferTargetZA(mat, projPDG);
    return IdentifyChannel(procName, projPDG, tZ, tA);
}

//──────────────────────────────────────────────────────────────────────────────
void SteppingAction::UserSteppingAction(const G4Step* step)
{
    G4Track* track = step->GetTrack();
    ++fTotalSteps;

    if (track->GetParentID() == 0) return;        // skip primary
    if (track->GetCurrentStepNumber() != 1) return;

    const G4String& pname = track->GetDefinition()->GetParticleName();
    const G4int evtID = G4RunManager::GetRunManager()
                            ->GetCurrentEvent()->GetEventID();

    // Diagnostic: all secondaries, first 5 events
    if (evtID < 5) {
        const G4VProcess* proc = track->GetCreatorProcess();
        G4cout << "[DIAG evt=" << evtID << "] secondary: " << pname
               << "  E=" << track->GetKineticEnergy()/MeV << " MeV"
               << "  proc=" << (proc ? proc->GetProcessName() : "?")
               << "  parentID=" << track->GetParentID()
               << "  vol=" << step->GetPreStepPoint()->GetPhysicalVolume()->GetName()
               << G4endl;
    }

    // ── Score neutrons ────────────────────────────────────────────────────
    // UNCHANGED — neutron scoring is not affected by the alpha fix.
    if (pname == "neutron") {
        const G4int trackID = track->GetTrackID();
        if (fRecordedTracks.count(trackID)) return;
        fRecordedTracks.insert(trackID);

        NeutronChannel channel = ClassifySecondary(track, step);

        const G4double      kineticE = step->GetPreStepPoint()->GetKineticEnergy();
        const G4ThreeVector dir      = step->GetPreStepPoint()->GetMomentumDirection();
        const G4ThreeVector origin   = track->GetVertexPosition();
        const G4VProcess*   creator  = track->GetCreatorProcess();
        const G4String      procName = creator ? creator->GetProcessName() : "?";

        ++fNeutronsSeen;
        ++fChannelNeutrons[static_cast<int>(channel)];

        if (fNeutronsSeen <= 20) {
            G4cout << ChannelColour(channel)
                   << "[NEUTRON #" << fNeutronsSeen
                   << "] ch=" << ChannelName(channel)
                   << "  E=" << kineticE/MeV << " MeV"
                   << "  θ=" << dir.theta()/CLHEP::deg << "°"
                   << "  proc=" << procName
                   << kReset << G4endl;
        }

        fEventAction->RecordNeutron(kineticE, dir.theta(), dir.phi(),
                                    origin, procName, channel);
        return;
    }

    // ── Score alpha particles ─────────────────────────────────────────────
    // FIX: only score alphas that come EXCLUSIVELY from B11_pa_TALYS.
    //
    // Previously the check used proc.find("Inelastic") which matched every
    // G4HadronInelasticProcess in this physics list (B11_pn_TALYS, B10_pn_TALYS,
    // N14_pn_TALYS, B11_an_TALYS, B10_an_TALYS, N14_an_TALYS, B11_pa_TALYS) —
    // all of them are inelastic processes and all can produce alpha recoils.
    // This inflated fAlphasSeen and the B11_p2alpha channel counter.
    //
    // Now we require the exact process name. Alphas that are recoil fragments
    // from p,n or α,n reactions are intentionally NOT counted here; they are
    // secondary projectiles for the (α,n) channels and are handled separately
    // by Geant4 transporting them as new primaries.
    if (pname == "alpha") {
        const G4VProcess* creator = track->GetCreatorProcess();
        if (!creator) return;
        const G4String& proc = creator->GetProcessName();

        // EXACT match only — do not use substring search.
        if (proc == "B11_pa_TALYS") {
            ++fAlphasSeen;
            ++fChannelAlphas[static_cast<int>(NeutronChannel::B11_p2alpha)];
            fEventAction->RecordAlpha(track->GetKineticEnergy(),
                                      NeutronChannel::B11_p2alpha);

            if (fAlphasSeen <= 10) {
                G4cout << ChannelColour(NeutronChannel::B11_p2alpha)
                       << "[ALPHA #" << fAlphasSeen
                       << "] E=" << track->GetKineticEnergy()/MeV << " MeV"
                       << "  proc=" << proc
                       << kReset << G4endl;
            }
        }
        // All other alpha producers (recoils from p,n, α,n, etc.) are silently
        // ignored here. They are transported by Geant4 as tracks and will
        // naturally trigger the B11_an_TALYS / B10_an_TALYS / N14_an_TALYS
        // processes when they reach the catcher material.
    }
}

void SteppingAction::NewEvent()
{ fRecordedTracks.clear(); fDiagPrinted = 0; }

void SteppingAction::PrintSummary() const
{
    G4cout << "\n[SteppingAction] Total steps : " << fTotalSteps
           << "\n                 Neutrons    : " << fNeutronsSeen
           << "\n                 Alphas (B11_pa only): " << fAlphasSeen << "\n";
    for (int i = 0; i < kNChannels; ++i) {
        if (fChannelNeutrons[i] > 0) {
            auto ch = static_cast<NeutronChannel>(i);
            G4cout << ChannelColour(ch) << "  " << ChannelName(ch)
                   << ": " << fChannelNeutrons[i] << kReset << "\n";
        }
    }
    G4cout << G4endl;
}
