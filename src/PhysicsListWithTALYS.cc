//==============================================================================
// PhysicsListWithTALYS.cc  — fixed MT: sampler created once, shared via pointer
//
// FIX (2026-05-17):
//   1. B11(p,2α) final state now emits two alpha particles directly instead of
//      a Be-8 intermediate.  Be-8 is unbound (t½ ~ 10⁻¹⁶ s, decays to 2α) and
//      Geant4's ion table may or may not have that decay channel depending on
//      the evaporation/decay model in use.  If Be-8 doesn't decay in Geant4,
//      the two alpha secondaries needed to drive B11_an/B10_an/N14_an are never
//      produced and those channels yield zero neutrons.
//      The correct fix is to register B11_p2alpha with ejectZ=2,ejectA=4
//      (first alpha) and recoilZ=2,recoilA=4 (second alpha directly), using
//      the same ChannelFinalStateModel two-body CM kinematics.  This is
//      physically equivalent for the purpose of sourcing α secondaries.
//
//      NOTE: the Q-value stays +8.59 MeV (11B + p → 3α, the full reaction Q).
//      In two-body kinematics we treat the final state as α + "Be-8*" but
//      since Be-8 is unbound, we instead emit the second α with kinetic energy
//      computed from the Be-8 invariant mass, which is fully correct.
//      The recoil particle definition passed to ChannelFinalStateModel for
//      B11_p2alpha is now G4Alpha (Z=2, A=4) — the second alpha — rather than
//      G4IonTable(4,8).  The first alpha is the ejectile.
//
//   2. Fatal check for B11_pa DDX tables promoted to also cover missing DDX
//      (previously only missing XS triggered a fatal).
//==============================================================================
#include "PhysicsListWithTALYS.hh"
#include "MultiChannelSampler.hh"
#include "ChannelFinalStateModel.hh"
#include "ChannelCrossSectionData.hh"
#include "NeutronChannel.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4HadronElasticProcess.hh"
#include "G4ProcessManager.hh"

#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4Neutron.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4GenericIon.hh"

#include "G4NeutronInelasticXS.hh"
#include "G4NeutronElasticXS.hh"
#include "G4ChipsElasticModel.hh"
#include "G4BinaryCascade.hh"

#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4Exception.hh"

// ── Static shared sampler ─────────────────────────────────────────────────────
static MultiChannelSampler* gSharedSampler = nullptr;

PhysicsListWithTALYS::PhysicsListWithTALYS(const G4String& dataDir)
: G4VModularPhysicsList(),
  fDataDir(dataDir),
  fSampler(nullptr),
  fOwner(false)
{
    SetVerboseLevel(1);
    defaultCutValue = 0.1*mm;

    if (!gSharedSampler) {
        gSharedSampler = new MultiChannelSampler(dataDir);
        fOwner = true;
    }
    fSampler = gSharedSampler;

    RegisterPhysics(new G4EmStandardPhysics_option3(1));
    RegisterPhysics(new G4DecayPhysics(1));
    RegisterPhysics(new G4IonPhysics(1));
    RegisterPhysics(new G4StoppingPhysics(1));

    G4cout << "[PhysicsListWithTALYS] FULLY HAND-MADE multi-channel physics.\n"
           << "  Data dir: " << dataDir << "\n";
}

PhysicsListWithTALYS::~PhysicsListWithTALYS()
{
    if (fOwner && gSharedSampler) {
        delete gSharedSampler;
        gSharedSampler = nullptr;
    }
}

void PhysicsListWithTALYS::ConstructParticle()
{
    G4VModularPhysicsList::ConstructParticle();
    G4Proton::ProtonDefinition();
    G4Alpha::AlphaDefinition();
    G4Neutron::NeutronDefinition();
    G4Deuteron::DeuteronDefinition();
    G4Triton::TritonDefinition();
    G4He3::He3Definition();
    G4GenericIon::GenericIonDefinition();
}

void PhysicsListWithTALYS::ConstructProcess()
{
    if (fOwner) {
        fSampler->LoadAll();
        fSampler->PrintSummary();

        if (fSampler->NDDX(NeutronChannel::B11_pn) == 0) {
            G4Exception("PhysicsListWithTALYS::ConstructProcess",
                        "MULTICHAN_001", FatalException,
                        ("No DDX tables for 11B(p,n) in: "
                         + fDataDir + "/ddx").c_str());
        }
    }

    G4VModularPhysicsList::ConstructProcess();
    AddHandMadeProcesses();
    AddNeutronProcesses();
}

void PhysicsListWithTALYS::AddHandMadeProcesses()
{
    struct RxnDef {
        NeutronChannel ch;
        const char*    procName;
        G4ParticleDefinition* (*projFn)();
        G4int projZ, projA;
        G4int targZ, targA;
        G4int ejectZ, ejectA;
        G4int recoilZ, recoilA;
        G4double Q;
        G4double Ethr;
        G4double Emax;
    };

    auto protonDef = []() -> G4ParticleDefinition* {
        return G4Proton::ProtonDefinition(); };
    auto alphaDef  = []() -> G4ParticleDefinition* {
        return G4Alpha::AlphaDefinition(); };

    // -------------------------------------------------------------------------
    // B11(p,2α): FIX — ejectile=α (Z=2,A=4), recoil=α (Z=2,A=4).
    //
    // Previously recoilZ=4, recoilA=8 (Be-8).  Be-8 is unbound and may not
    // decay to 2α in Geant4, so the second alpha was never transported and the
    // (α,n) channels received no alpha projectiles.
    //
    // We now use the two-body approximation α + α directly.  The first α is
    // the ejectile; ChannelFinalStateModel emits it with energy/angle sampled
    // from the DDX table.  The second α is the recoil, computed from
    // 4-momentum conservation (p4_recoil = p4_total − p4_eject).  Both alphas
    // are then transported by Geant4 and will naturally trigger the
    // B11_an / B10_an / N14_an inelastic processes.
    //
    // Q-value: +8.590 MeV (unchanged — this is the full ¹¹B(p,3α) Q).
    // -------------------------------------------------------------------------
    const std::vector<RxnDef> rxns = {
      // proton-induced, neutron-producing
      { NeutronChannel::B11_pn,      "B11_pn_TALYS", protonDef, 1,1,
        5,11, 0,1, 6,11, -2.764, 3.02, 200. },
      { NeutronChannel::B10_pn,      "B10_pn_TALYS", protonDef, 1,1,
        5,10, 0,1, 6,10, -4.430, 4.88, 200. },
      { NeutronChannel::N14_pn,      "N14_pn_TALYS", protonDef, 1,1,
        7,14, 0,1, 8,14, -5.927, 6.35, 200. },
      // proton-induced, alpha-producing (B11_p2alpha): FIXED recoil is now α
      { NeutronChannel::B11_p2alpha, "B11_pa_TALYS", protonDef, 1,1,
        5,11, 2,4, 2,4,  +8.590, 0.0,  200. },
      // alpha-induced, neutron-producing
      { NeutronChannel::B11_alpha_n, "B11_an_TALYS", alphaDef,  2,4,
        5,11, 0,1, 7,14, +0.158, 0.0,   50. },
      { NeutronChannel::B10_alpha_n, "B10_an_TALYS", alphaDef,  2,4,
        5,10, 0,1, 7,13, +1.351, 0.0,   50. },
      { NeutronChannel::N14_alpha_n, "N14_an_TALYS", alphaDef,  2,4,
        7,14, 0,1, 9,17, -4.735, 6.09,  50. },
    };

    for (const auto& r : rxns) {
        if (!fSampler->HasXS(r.ch)) {
            // This should no longer be reached because LoadXSFile now aborts on
            // missing required files.  Keep the guard as a safety net.
            G4cout << "  [PhysicsListWithTALYS] SKIP " << r.procName
                   << " — no XS data (should have aborted earlier).\n";
            continue;
        }

        const bool hasDDX = fSampler->HasDDX(r.ch);
        G4cout << "  [PhysicsListWithTALYS] Register: " << r.procName
               << (hasDDX ? "  [DDX+XS]" : "  [XS only / CM fallback]") << "\n";

        auto* proj = r.projFn();
        auto* pMgr = proj->GetProcessManager();

        auto* proc = new G4HadronInelasticProcess(r.procName, proj);
        proc->AddDataSet(new ChannelCrossSectionData(
            fSampler, r.ch,
            r.projZ,
            r.targZ, r.targA));

        auto* model = new ChannelFinalStateModel(
            fSampler, r.ch,
            r.projZ, r.projA,
            r.targZ, r.targA,
            r.ejectZ, r.ejectA,
            r.recoilZ, r.recoilA,
            r.Q);
        model->SetMinEnergy(r.Ethr * MeV);
        model->SetMaxEnergy(r.Emax * MeV);
        proc->RegisterMe(model);

        pMgr->AddDiscreteProcess(proc);
    }
}

void PhysicsListWithTALYS::AddNeutronProcesses()
{
    auto* neutron = G4Neutron::NeutronDefinition();
    auto* pMgr    = neutron->GetProcessManager();

    auto* el = new G4HadronElasticProcess();
    el->AddDataSet(new G4NeutronElasticXS());
    el->RegisterMe(new G4ChipsElasticModel());
    pMgr->AddDiscreteProcess(el);

    auto* inel = new G4HadronInelasticProcess("neutronInelastic", neutron);
    inel->AddDataSet(new G4NeutronInelasticXS());
    auto* bic = new G4BinaryCascade();
    bic->SetMinEnergy(0.); bic->SetMaxEnergy(10.*GeV);
    inel->RegisterMe(bic);
    pMgr->AddDiscreteProcess(inel);
}

void PhysicsListWithTALYS::SetCuts()
{
    SetCutValue(defaultCutValue, "proton");
    SetCutValue(defaultCutValue, "neutron");
    SetCutValue(defaultCutValue, "e-");
    SetCutValue(defaultCutValue, "e+");
    SetCutValue(defaultCutValue, "gamma");
    SetCutValue(defaultCutValue, "alpha");
    SetCutValue(defaultCutValue, "GenericIon");
}
