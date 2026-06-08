//==============================================================================
// EventAction.cc — per-channel histogram filling
//
// H1 layout (booked in RunAction):
//   id=0  h1_energy_all      all neutrons, energy
//   id=1  h1_theta_all       all neutrons, polar angle
//   id=2  h1_phi_all         all neutrons, azimuthal
//   id=3  h1_yield           1-bin yield summary (filled in EndOfRunAction)
//   id=4  h1_proton_energy   primary proton energy spectrum
//   -- per channel (ch=0..kNChannels-1): base = 5 + ch*3 --
//   id=base+0  h1_energy_<tag>
//   id=base+1  h1_theta_<tag>
//   id=base+2  h1_phi_<tag>
//
// H2 layout:
//   id=0  h2_E_theta_all   all neutrons
//   id=1+ch  h2_E_theta_<tag>  per channel
//
// Ntuple columns:
//   0: E_MeV        double
//   1: theta_deg    double
//   2: phi_deg      double
//   3: originX_mm   double
//   4: originY_mm   double
//   5: originZ_mm   double
//   6: process      string
//   7: channel_id   int
//   8: channel_name string
//   9: is_neutron   int    (1=neutron, 0=alpha for parasitic tracking)
//==============================================================================
#include "EventAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"

#include "G4Event.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

EventAction::EventAction(RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction)
{
    fChannelNeutrons.fill(0);
    fChannelAlphas.fill(0);
}

EventAction::~EventAction() {}

//──────────────────────────────────────────────────────────────────────────────
void EventAction::BeginOfEventAction(const G4Event*)
{
    fNeutronCount = 0;
    fAlphaCount   = 0;
    fDiagCount    = 0;
    fChannelNeutrons.fill(0);
    fChannelAlphas.fill(0);
    if (fSteppingAction) fSteppingAction->NewEvent();
}

//──────────────────────────────────────────────────────────────────────────────
void EventAction::RecordNeutron(G4double        kineticE,
                                 G4double        theta,
                                 G4double        phi,
                                 const G4ThreeVector& origin,
                                 const G4String& process,
                                 NeutronChannel  channel)
{
    const G4double E_MeV    = kineticE / MeV;
    G4double theta_deg = theta / deg;
    G4double phi_deg   = phi   / deg;
    if (phi_deg < 0.) phi_deg += 360.;

    const G4int chIdx = static_cast<int>(channel);
    auto* am = G4AnalysisManager::Instance();

    // All-channel combined H1
    am->FillH1(0, E_MeV);
    am->FillH1(1, theta_deg);
    am->FillH1(2, phi_deg);

    // Per-channel H1  (base = 5 + chIdx*3)
    const G4int base = 5 + chIdx * 3;
    am->FillH1(base + 0, E_MeV);
    am->FillH1(base + 1, theta_deg);
    am->FillH1(base + 2, phi_deg);

    // H2
    am->FillH2(0, E_MeV, theta_deg);
    am->FillH2(1 + chIdx, E_MeV, theta_deg);

    // Ntuple
    am->FillNtupleDColumn(0, E_MeV);
    am->FillNtupleDColumn(1, theta_deg);
    am->FillNtupleDColumn(2, phi_deg);
    am->FillNtupleDColumn(3, origin.x()/mm);
    am->FillNtupleDColumn(4, origin.y()/mm);
    am->FillNtupleDColumn(5, origin.z()/mm);
    am->FillNtupleSColumn(6, process);
    am->FillNtupleIColumn(7, chIdx);
    am->FillNtupleSColumn(8, ChannelName(channel));
    am->FillNtupleIColumn(9, 1);   // is_neutron=1
    am->AddNtupleRow();

    ++fNeutronCount;
    ++fChannelNeutrons[chIdx];

    if (fDiagCount < 10) {
        G4cout << ChannelColour(channel)
               << "[Neutron] " << ChannelName(channel)
               << "  E=" << E_MeV << " MeV"
               << "  θ=" << theta_deg << "°"
               << "  proc=" << process
               << kReset << G4endl;
        ++fDiagCount;
    }
}

//──────────────────────────────────────────────────────────────────────────────
void EventAction::RecordAlpha(G4double kineticE, NeutronChannel channel)
{
    // Record alpha from ¹¹B(p,2α) in the ntuple for parasitic-neutron analysis
    const G4double E_MeV = kineticE / MeV;
    const G4int chIdx    = static_cast<int>(channel);

    auto* am = G4AnalysisManager::Instance();
    am->FillNtupleDColumn(0, E_MeV);
    am->FillNtupleDColumn(1, -1.);   // theta undefined for alphas in this context
    am->FillNtupleDColumn(2, -1.);
    am->FillNtupleDColumn(3, 0.); am->FillNtupleDColumn(4, 0.); am->FillNtupleDColumn(5, 0.);
    am->FillNtupleSColumn(6, "protonInelastic");
    am->FillNtupleIColumn(7, chIdx);
    am->FillNtupleSColumn(8, ChannelName(channel));
    am->FillNtupleIColumn(9, 0);   // is_neutron=0 (alpha)
    am->AddNtupleRow();

    ++fAlphaCount;
    ++fChannelAlphas[chIdx];
}

//──────────────────────────────────────────────────────────────────────────────
void EventAction::EndOfEventAction(const G4Event*)
{
    fRunAction->AddNeutrons(fNeutronCount);
    fRunAction->AddProtons(1);
    fRunAction->AddAlphas(fAlphaCount);

    for (int i = 0; i < kNChannels; ++i) {
        fRunAction->AddChannelNeutrons(static_cast<NeutronChannel>(i),
                                       fChannelNeutrons[i]);
        fRunAction->AddChannelAlphas (static_cast<NeutronChannel>(i),
                                       fChannelAlphas[i]);
    }
}
