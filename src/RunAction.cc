#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4AccumulableManager.hh"
#include "G4ios.hh"

#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <string>

namespace {
std::string SanitizeForFileName(std::string s)
{
    for (char& c : s) {
        if (c == '.') c = 'p';
        else if (c == '+' || c == ' ' || c == '/' || c == '\\' || c == ':' || c == ';') c = '_';
    }
    return s;
}
}

//------------------------------------------------------------------------------
RunAction::RunAction()
: G4UserRunAction()
{
    auto* accMgr = G4AccumulableManager::Instance();

    accMgr->Register(fTotalNeutrons);
    accMgr->Register(fTotalProtons);
    accMgr->Register(fTotalAlphas);

    for (auto& acc : fChannelNeutrons) {
        acc = G4Accumulable<G4long>(0);
        accMgr->Register(acc);
    }

    for (auto& acc : fChannelAlphas) {
        acc = G4Accumulable<G4long>(0);
        accMgr->Register(acc);
    }

    auto* am = G4AnalysisManager::Instance();
    am->SetDefaultFileType("root");
    am->SetVerboseLevel(1);
    am->SetNtupleMerging(true);
}

RunAction::~RunAction() = default;

//------------------------------------------------------------------------------
G4String RunAction::BuildOutputFileName(const G4Run* run) const
{
    // Priority 1: explicit C++ setter.
    if (!fOutputFileName.empty()) {
        return fOutputFileName;
    }

    // Priority 2: shell/environment override. This is the safest method for
    // one-energy-per-process scans:
    //   PBN_OUTPUT_FILE=../results_scan/B11_multichannel_9p0MeV.root ./pBN ../run.mac
    if (const char* forced = std::getenv("PBN_OUTPUT_FILE")) {
        if (forced[0] != '\0') {
            return G4String(forced);
        }
    }

    // Priority 3: energy passed by shell. This does not depend on the state of
    // PrimaryGeneratorAction on the master thread.
    if (const char* eEnv = std::getenv("PBN_E_MEV")) {
        if (eEnv[0] != '\0') {
            std::ostringstream oss;
            oss << "B11_multichannel_"
                << SanitizeForFileName(eEnv)
                << "MeV_run" << run->GetRunID() << ".root";
            return G4String(oss.str());
        }
    }

    // Priority 4: try to get energy from PrimaryGeneratorAction. This may fail
    // in some MT/master setups, so it is intentionally only a fallback.
    G4double eMeV = -1.0;
    G4String mode = "unknown";

    const auto* gen = dynamic_cast<const PrimaryGeneratorAction*>(
        G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

    if (gen) {
        mode = gen->GetMode();
        if (mode == "mono") {
            eMeV = gen->GetEMono() / MeV;
        } else {
            eMeV = gen->GetTp() / MeV;
        }
    }

    std::ostringstream oss;
    if (eMeV > 0.0) {
        oss << "B11_multichannel_"
            << std::fixed << std::setprecision(2)
            << eMeV << "MeV_run" << run->GetRunID() << ".root";
    } else {
        oss << "B11_multichannel_unknownE_run" << run->GetRunID() << ".root";
    }

    return G4String(oss.str());
}

//------------------------------------------------------------------------------
void RunAction::BookHistogramsAndNtuple()
{
    auto* am = G4AnalysisManager::Instance();

    // H1: all-channel combined, ids 0..4.
    am->CreateH1("h1_energy_all",
                 "Neutron energy (all ch);E_{n} [MeV];Counts/event",
                 500, 0., 50.);                              // id=0

    am->CreateH1("h1_theta_all",
                 "Neutron polar angle (all ch);#theta [deg];Counts/event",
                 180, 0., 180.);                             // id=1

    am->CreateH1("h1_phi_all",
                 "Neutron azimuthal angle (all ch);#phi [deg];Counts/event",
                 360, 0., 360.);                             // id=2

    am->CreateH1("h1_yield",
                 "Neutron yield summary;bin;Y_{n/p}",
                 1, 0., 1.);                                 // id=3

    // Wider range for your 3–20 MeV scan.
    am->CreateH1("h1_proton_energy",
                 "Primary proton energy;E_{p} [MeV];Counts",
                 250, 0., 25.);                              // id=4

    // Per-channel H1 histograms.
    for (int i = 0; i < kNChannels; ++i) {
        auto     ch   = static_cast<NeutronChannel>(i);
        G4String tag  = ChannelTag(ch);
        G4String name = ChannelName(ch);

        am->CreateH1("h1_energy_" + tag,
                     "Neutron energy " + name + ";E_{n} [MeV];Counts",
                     500, 0., 50.);

        am->CreateH1("h1_theta_" + tag,
                     "Neutron polar angle " + name + ";#theta [deg];Counts",
                     180, 0., 180.);

        am->CreateH1("h1_phi_" + tag,
                     "Neutron azimuthal angle " + name + ";#phi [deg];Counts",
                     360, 0., 360.);
    }

    // H2: all-channel and per-channel E-theta maps.
    am->CreateH2("h2_E_theta_all",
                 "E vs #theta (all ch);E_{n} [MeV];#theta [deg]",
                 250, 0., 50., 90, 0., 180.);

    for (int i = 0; i < kNChannels; ++i) {
        auto     ch   = static_cast<NeutronChannel>(i);
        G4String tag  = ChannelTag(ch);
        G4String name = ChannelName(ch);

        am->CreateH2("h2_E_theta_" + tag,
                     "E vs #theta " + name + ";E_{n} [MeV];#theta [deg]",
                     250, 0., 50., 90, 0., 180.);
    }

    // Ntuple.
    am->CreateNtuple("particles", "Scored neutrons and alphas");
    am->CreateNtupleDColumn("E_MeV");        // 0
    am->CreateNtupleDColumn("theta_deg");    // 1
    am->CreateNtupleDColumn("phi_deg");      // 2
    am->CreateNtupleDColumn("originX_mm");   // 3
    am->CreateNtupleDColumn("originY_mm");   // 4
    am->CreateNtupleDColumn("originZ_mm");   // 5
    am->CreateNtupleSColumn("process");      // 6
    am->CreateNtupleIColumn("channel_id");   // 7
    am->CreateNtupleSColumn("channel_name"); // 8
    am->CreateNtupleIColumn("is_neutron");   // 9
    am->FinishNtuple();
}

//------------------------------------------------------------------------------
void RunAction::BeginOfRunAction(const G4Run* run)
{
    G4AccumulableManager::Instance()->Reset();

    const G4String outFile = BuildOutputFileName(run);

    auto* am = G4AnalysisManager::Instance();
    am->OpenFile(outFile);

    // Important: do not protect this with IsMaster(). Worker analysis managers
    // also need the histogram/ntuple definitions because workers fill them.
    BookHistogramsAndNtuple();

    if (IsMaster()) {
        G4cout << "\n=== Run " << run->GetRunID()
               << " -> " << outFile << " ===\n";
    }
}

//------------------------------------------------------------------------------
void RunAction::EndOfRunAction(const G4Run* run)
{
    if (run->GetNumberOfEvent() == 0) {
        return;
    }

    G4AccumulableManager::Instance()->Merge();

    const G4long   nP     = fTotalProtons.GetValue();
    const G4long   nN     = fTotalNeutrons.GetValue();
    const G4long   nA     = fTotalAlphas.GetValue();
    const G4double yieldN = (nP > 0) ? static_cast<G4double>(nN) / nP : 0.0;

    auto* am = G4AnalysisManager::Instance();

    if (IsMaster()) {
        am->FillH1(3, 0.5, yieldN);

        // Separate primary proton-induced neutron yield from secondary alpha-induced yield.
        // This does not change physics or counters; it only changes the printed report.
        const G4long nB11pn = fChannelNeutrons[static_cast<int>(NeutronChannel::B11_pn)].GetValue();
        const G4long nB10pn = fChannelNeutrons[static_cast<int>(NeutronChannel::B10_pn)].GetValue();
        const G4long nN14pn = fChannelNeutrons[static_cast<int>(NeutronChannel::N14_pn)].GetValue();

        const G4long nB11an = fChannelNeutrons[static_cast<int>(NeutronChannel::B11_alpha_n)].GetValue();
        const G4long nB10an = fChannelNeutrons[static_cast<int>(NeutronChannel::B10_alpha_n)].GetValue();
        const G4long nN14an = fChannelNeutrons[static_cast<int>(NeutronChannel::N14_alpha_n)].GetValue();

        const G4long nPN = nB11pn + nB10pn + nN14pn;
        const G4long nAN = nB11an + nB10an + nN14an;

        auto printNeutronChannel = [&](NeutronChannel ch, G4long nCh)
        {
            const G4double fracAll = (nN > 0) ? 100.0 * nCh / nN : 0.0;
            const G4double yieldCh = (nP > 0) ? static_cast<G4double>(nCh) / nP : 0.0;

            G4cout << "║  " << ChannelColour(ch)
                   << std::left << std::setw(18) << ChannelName(ch) << kReset
                   << std::right
                   << std::setw(10) << nCh
                   << "    " << std::scientific << std::setprecision(3) << yieldCh
                   << "   " << std::fixed << std::setprecision(1)
                   << std::setw(5) << fracAll << "%"
                   << "       ║\n";
        };

        auto printGroupLine = [&](const char* label, G4long nGroup)
        {
            const G4double yGroup = (nP > 0) ? static_cast<G4double>(nGroup) / nP : 0.0;
            const G4double fGroup = (nN > 0) ? 100.0 * nGroup / nN : 0.0;

            G4cout << "║  " << std::left << std::setw(18) << label
                   << std::right
                   << std::setw(10) << nGroup
                   << "    " << std::scientific << std::setprecision(3) << yGroup
                   << "   " << std::fixed << std::setprecision(1)
                   << std::setw(5) << fGroup << "%"
                   << "       ║\n";
        };

        G4cout << "\n╔══════════════════════════════════════════════════════════════╗\n"
               << "║  Run " << std::setw(4) << run->GetRunID()
               << " — Multi-channel p + BN_CATCHER summary              ║\n"
               << "╠══════════════════════════════════════════════════════════════╣\n"
               << std::left
               << "║  Primary protons   : " << std::setw(10) << nP
               << "                               ║\n"
               << "║  Total neutrons    : " << std::setw(10) << nN
               << "  Y(n/p) = " << std::scientific << std::setprecision(3)
               << yieldN << "         ║\n"
               << "║  Alpha secondaries : " << std::setw(10) << nA
               << "                               ║\n"
               << "╠══════════════════════════════════════════════════════════════╣\n"
               << "║  Channel         N(neutrons)    Y(n/p)       frac(all)      ║\n"
               << "╠══════════════════════════════════════════════════════════════╣\n"
               << "║  Primary proton-induced neutrons                            ║\n"
               << "╠══════════════════════════════════════════════════════════════╣\n";

        printNeutronChannel(NeutronChannel::B11_pn, nB11pn);
        printNeutronChannel(NeutronChannel::B10_pn, nB10pn);
        printNeutronChannel(NeutronChannel::N14_pn, nN14pn);
        printGroupLine("Total p,n", nPN);

        G4cout << "╠══════════════════════════════════════════════════════════════╣\n"
               << "║  Secondary alpha-induced neutrons                           ║\n"
               << "╠══════════════════════════════════════════════════════════════╣\n";

        printNeutronChannel(NeutronChannel::B11_alpha_n, nB11an);
        printNeutronChannel(NeutronChannel::B10_alpha_n, nB10an);
        printNeutronChannel(NeutronChannel::N14_alpha_n, nN14an);
        printGroupLine("Total alpha,n", nAN);

        G4cout << "╠══════════════════════════════════════════════════════════════╣\n";

        const G4long nAlphaFromP2A =
            fChannelAlphas[static_cast<int>(NeutronChannel::B11_p2alpha)].GetValue();

        G4cout << "║  " << ChannelColour(NeutronChannel::B11_p2alpha)
               << std::left << std::setw(18) << ChannelName(NeutronChannel::B11_p2alpha)
               << kReset << "  [no neutron]   alphas=" << std::setw(8)
               << nAlphaFromP2A << "             ║\n";

        // Safety check: if any future neutron-producing channel is added but not
        // included in the grouped table above, print it in an extra section.
        bool extraHeaderPrinted = false;
        for (int i = 0; i < kNChannels; ++i) {
            auto ch = static_cast<NeutronChannel>(i);
            if (ch == NeutronChannel::B11_pn ||
                ch == NeutronChannel::B10_pn ||
                ch == NeutronChannel::N14_pn ||
                ch == NeutronChannel::B11_alpha_n ||
                ch == NeutronChannel::B10_alpha_n ||
                ch == NeutronChannel::N14_alpha_n ||
                ch == NeutronChannel::B11_p2alpha) {
                continue;
            }

            const G4long nCh = fChannelNeutrons[i].GetValue();
            if (nCh <= 0) continue;

            if (!extraHeaderPrinted) {
                G4cout << "╠══════════════════════════════════════════════════════════════╣\n"
                       << "║  Other neutron-producing channels                           ║\n"
                       << "╠══════════════════════════════════════════════════════════════╣\n";
                extraHeaderPrinted = true;
            }
            printNeutronChannel(ch, nCh);
        }

        G4cout << "╚══════════════════════════════════════════════════════════════╝\n\n";
    }

    am->Write();
    am->CloseFile();
}

