#pragma once

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "G4String.hh"
#include "NeutronChannel.hh"

#include <array>

class G4Run;

class RunAction : public G4UserRunAction {
public:
    RunAction();
    ~RunAction() override;

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction  (const G4Run*) override;

    // Optional direct filename override. If this is not set, RunAction will
    // check the environment variable PBN_OUTPUT_FILE. If that is also absent,
    // it will build a filename from PBN_E_MEV or from PrimaryGeneratorAction.
    void SetOutputFileName(const G4String& name) { fOutputFileName = name; }

    void AddNeutrons(G4long n) { fTotalNeutrons += n; }
    void AddProtons (G4long n) { fTotalProtons  += n; }
    void AddAlphas  (G4long n) { fTotalAlphas   += n; }

    void AddChannelNeutrons(NeutronChannel ch, G4long n)
    {
        const int idx = static_cast<int>(ch);
        if (idx >= 0 && idx < static_cast<int>(fChannelNeutrons.size())) {
            fChannelNeutrons[idx] += n;
        }
    }

    void AddChannelAlphas(NeutronChannel ch, G4long n)
    {
        const int idx = static_cast<int>(ch);
        if (idx >= 0 && idx < static_cast<int>(fChannelAlphas.size())) {
            fChannelAlphas[idx] += n;
        }
    }

private:
    void BookHistogramsAndNtuple();
    G4String BuildOutputFileName(const G4Run* run) const;

    G4String fOutputFileName;

    G4Accumulable<G4long> fTotalNeutrons{0};
    G4Accumulable<G4long> fTotalProtons {0};
    G4Accumulable<G4long> fTotalAlphas  {0};

    // Size 8 includes OTHER at index 7. Only 0..kNChannels-1 are printed/booked
    // as normal physics channels.
    std::array<G4Accumulable<G4long>, 8> fChannelNeutrons;
    std::array<G4Accumulable<G4long>, 8> fChannelAlphas;
};
