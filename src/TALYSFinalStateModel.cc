#include "TALYSFinalStateModel.hh"
#include "TALYSNeutronSampler.hh"

#include "G4DynamicParticle.hh"
#include "G4IonTable.hh"
#include "G4LorentzVector.hh"
#include "G4Neutron.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Exception.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cmath>
#include <string>

namespace {
constexpr G4double kThresholdMeV = 3.017064;

G4double SafeSqrt(G4double x)
{
    return std::sqrt(std::max(0.0, x));
}
}

TALYSFinalStateModel::TALYSFinalStateModel(TALYSNeutronSampler* sampler)
: G4HadronicInteraction("TALYSFinalStateModel_B11pn"),
  fSampler(sampler)
{
    SetMinEnergy(kThresholdMeV * MeV);
    SetMaxEnergy(200.0 * MeV);
}

G4bool TALYSFinalStateModel::IsApplicable(const G4HadProjectile& projectile,
                                          G4Nucleus& targetNucleus)
{
    const auto* def = projectile.GetDefinition();
    return def && def->GetParticleName() == "proton"
           && targetNucleus.GetZ_asInt() == 5
           && targetNucleus.GetA_asInt() == 11;
}

G4ThreeVector TALYSFinalStateModel::DirectionFromMuPhi(const G4ThreeVector& axis,
                                                       G4double mu,
                                                       G4double phi) const
{
    mu = std::max(-1.0, std::min(1.0, mu));
    const G4double sinTheta = std::sqrt(std::max(0.0, 1.0 - mu*mu));
    const G4ThreeVector local(sinTheta*std::cos(phi), sinTheta*std::sin(phi), mu);

    G4ThreeVector z = axis.unit();
    if (z.mag2() == 0.0) return local.unit();

    const G4ThreeVector ref = (std::abs(z.z()) < 0.999) ? G4ThreeVector(0.,0.,1.)
                                                         : G4ThreeVector(0.,1.,0.);
    const G4ThreeVector x = (ref.cross(z)).unit();
    const G4ThreeVector y = (z.cross(x)).unit();

    return (local.x()*x + local.y()*y + local.z()*z).unit();
}

G4HadFinalState* TALYSFinalStateModel::ApplyYourself(const G4HadProjectile& projectile,
                                                     G4Nucleus&)
{
    fLocalState.Clear();
    fLocalState.SetStatusChange(stopAndKill);

    const G4double Ep_MeV = projectile.GetKineticEnergy() / MeV;

    auto* neutronDef = G4Neutron::Neutron();
    auto* ionTable   = G4IonTable::GetIonTable();
    auto* B11Def     = ionTable->GetIon(5, 11, 0.0);
    auto* C11Ground  = ionTable->GetIon(6, 11, 0.0);

    if (!neutronDef || !B11Def || !C11Ground) {
        return &fLocalState;
    }

    const G4double mn   = neutronDef->GetPDGMass();
    const G4double mB11 = B11Def->GetPDGMass();
    const G4double mC11 = C11Ground->GetPDGMass();

    const G4LorentzVector p4total = projectile.Get4Momentum()
                                  + G4LorentzVector(0., 0., 0., mB11);
    const G4ThreeVector protonDir = projectile.Get4Momentum().vect().unit();

    auto addFinalState = [&](const G4LorentzVector& p4n,
                             const G4LorentzVector& p4c,
                             G4double excitation) -> G4bool {
        const G4double nKE = p4n.e() - mn;
        if (nKE <= 0.0 || p4n.vect().mag2() <= 0.0) return false;

        excitation = std::max(0.0, excitation);
        G4ParticleDefinition* C11 = ionTable->GetIon(6, 11, excitation);
        if (!C11) C11 = C11Ground;

        const G4double cMass = C11->GetPDGMass();
        const G4double cKE = p4c.e() - cMass;
        if (cKE < -1.e-6*MeV) return false;

        auto* neutron = new G4DynamicParticle(neutronDef, p4n.vect().unit(), nKE);
        fLocalState.AddSecondary(neutron);

        if (p4c.vect().mag2() > 0.0 && cKE > 0.0) {
            auto* recoil = new G4DynamicParticle(C11, p4c.vect().unit(), std::max(0.0, cKE));
            fLocalState.AddSecondary(recoil);
        }
        return true;
    };

    // First try to use the sampled TALYS/TENDL neutron energy-angle point.
    // Energy conservation is enforced by assigning the remaining invariant mass
    // to the C-11 recoil as excitation energy. If the sampled point is outside
    // the kinematically allowed region, resample. No analytic fallback is used.
    for (G4int trial = 0; trial < 64; ++trial) {
        TALYSNeutronSample sample = fSampler ? fSampler->Sample(Ep_MeV) : TALYSNeutronSample{};
        if (!sample.valid) break;

        const G4double En = std::max(1.e-9*MeV, sample.En_MeV * MeV);
        const G4double pn = SafeSqrt(En * (En + 2.0*mn));
        const G4ThreeVector nDir = DirectionFromMuPhi(protonDir,
                                                      sample.mu,
                                                      sample.phi);
        const G4LorentzVector p4n(nDir * pn, En + mn);
        const G4LorentzVector p4c = p4total - p4n;

        const G4double mRec2 = p4c.m2();
        if (mRec2 < mC11*mC11) continue;

        const G4double mRec = SafeSqrt(mRec2);
        const G4double excitation = mRec - mC11;

        // Reject extremely large apparent excitation. This usually means the
        // sampled energy-angle point is incompatible with reaction kinematics.
        if (excitation > 50.0*MeV) continue;

        if (addFinalState(p4n, p4c, excitation)) {
            return &fLocalState;
        }
    }

    // Strict mode: do not generate an approximate two-body fallback neutron.
    // If no sampled TALYS/TENDL point can be used, stop the run and fix the
    // DDX table coverage/format instead of silently mixing model sources.
    G4Exception("TALYSFinalStateModel::ApplyYourself",
                "TALYS_FS_001",
                FatalException,
                ("No valid TALYS/TENDL DDX neutron sample for Ep = "
                 + std::to_string(Ep_MeV) + " MeV").c_str());

    return &fLocalState;
}
