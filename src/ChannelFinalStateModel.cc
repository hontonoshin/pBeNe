//==============================================================================
// ChannelFinalStateModel.cc
//==============================================================================
#include "ChannelFinalStateModel.hh"
#include "MultiChannelSampler.hh"

#include "G4DynamicParticle.hh"
#include "G4IonTable.hh"
#include "G4Neutron.hh"
#include "G4Alpha.hh"
#include "G4LorentzVector.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Exception.hh"
#include "Randomize.hh"

#include <cmath>
#include <algorithm>

static G4double SafeSqrt(G4double x)
{
    return std::sqrt(std::max(0.0, x));
}

//------------------------------------------------------------------------------
ChannelFinalStateModel::ChannelFinalStateModel(
    MultiChannelSampler* sampler,
    NeutronChannel       channel,
    G4int projZ,  G4int projA,
    G4int targZ,  G4int targA,
    G4int ejectZ, G4int ejectA,
    G4int recoilZ,G4int recoilA,
    G4double Qvalue_MeV)
: G4HadronicInteraction("ChannelFSM_" + ChannelTag(channel)),
  fSampler(sampler),
  fChannel(channel),
  fProjZ(projZ),
  fProjA(projA),
  fTargZ(targZ),
  fTargA(targA),
  fEjectZ(ejectZ),
  fEjectA(ejectA),
  fRecoilZ(recoilZ),
  fRecoilA(recoilA),
  fQvalue(Qvalue_MeV * MeV)
{}

//------------------------------------------------------------------------------
G4bool ChannelFinalStateModel::IsApplicable(
    const G4HadProjectile& proj,
    G4Nucleus& target)
{
    return proj.GetDefinition() != nullptr
        && target.GetZ_asInt() == fTargZ
        && target.GetA_asInt() == fTargA;
}

//------------------------------------------------------------------------------
G4ThreeVector ChannelFinalStateModel::DirectionFromMuPhi(
    const G4ThreeVector& axis,
    G4double mu,
    G4double phi) const
{
    mu = std::max(-1.0, std::min(1.0, mu));

    const G4double sinTh = SafeSqrt(1.0 - mu * mu);
    const G4ThreeVector local(
        sinTh * std::cos(phi),
        sinTh * std::sin(phi),
        mu
    );

    G4ThreeVector z = axis.unit();

    if (z.mag2() == 0.0) {
        return local.unit();
    }

    const G4ThreeVector ref =
        (std::abs(z.z()) < 0.999)
        ? G4ThreeVector(0, 0, 1)
        : G4ThreeVector(0, 1, 0);

    const G4ThreeVector x = (ref.cross(z)).unit();
    const G4ThreeVector y = (z.cross(x)).unit();

    return (local.x() * x + local.y() * y + local.z() * z).unit();
}

//------------------------------------------------------------------------------
G4ParticleDefinition* ChannelFinalStateModel::GetParticleDef(
    G4int Z,
    G4int A,
    G4double excitation) const
{
    if (Z == 0 && A == 1) {
        return G4Neutron::Neutron();
    }

    if (Z == 2 && A == 4) {
        return G4Alpha::AlphaDefinition();
    }

    return G4IonTable::GetIonTable()->GetIon(Z, A, excitation);
}

//------------------------------------------------------------------------------
G4bool ChannelFinalStateModel::AddTwoBodyFinalStateCM(
    const G4LorentzVector& p4tot,
    const G4ThreeVector& projDir,
    G4ParticleDefinition* ejectDef,
    G4ParticleDefinition* recoilDef)
{
    const G4double m1 = ejectDef->GetPDGMass();
    const G4double m2 = recoilDef->GetPDGMass();

    const G4double s = p4tot.m2();

    if (s <= 0.0) {
        return false;
    }

    const G4double sqrtS = std::sqrt(s);

    if (sqrtS <= m1 + m2) {
        return false;
    }

    const G4double lambda =
        (s - (m1 + m2) * (m1 + m2)) *
        (s - (m1 - m2) * (m1 - m2));

    if (lambda <= 0.0) {
        return false;
    }

    const G4double pStar = std::sqrt(lambda) / (2.0 * sqrtS);

    const G4double e1cm = std::sqrt(m1 * m1 + pStar * pStar);
    const G4double e2cm = std::sqrt(m2 * m2 + pStar * pStar);

    const G4double mu  = 2.0 * G4UniformRand() - 1.0;
    const G4double phi = twopi * G4UniformRand();

    const G4ThreeVector dirCM = DirectionFromMuPhi(projDir, mu, phi);

    G4LorentzVector p4eCM(dirCM * pStar, e1cm);
    G4LorentzVector p4rCM(-dirCM * pStar, e2cm);

    const G4ThreeVector beta = p4tot.boostVector();

    p4eCM.boost(beta);
    p4rCM.boost(beta);

    const G4double ke1 = std::max(0.0, p4eCM.e() - m1);
    const G4double ke2 = std::max(0.0, p4rCM.e() - m2);

    if (p4eCM.vect().mag2() <= 0.0 || p4rCM.vect().mag2() <= 0.0) {
        return false;
    }

    fLocalState.AddSecondary(
        new G4DynamicParticle(ejectDef, p4eCM.vect().unit(), ke1)
    );

    fLocalState.AddSecondary(
        new G4DynamicParticle(recoilDef, p4rCM.vect().unit(), ke2)
    );

    return true;
}

//------------------------------------------------------------------------------
G4HadFinalState* ChannelFinalStateModel::ApplyYourself(
    const G4HadProjectile& projectile,
    G4Nucleus&)
{
    fLocalState.Clear();
    fLocalState.SetStatusChange(stopAndKill);

    const G4double Einc_MeV = projectile.GetKineticEnergy() / MeV;

    auto* ionTable = G4IonTable::GetIonTable();

    G4ParticleDefinition* ejectDef =
        GetParticleDef(fEjectZ, fEjectA, 0.0);

    G4ParticleDefinition* targDef =
        ionTable->GetIon(fTargZ, fTargA, 0.0);

    G4ParticleDefinition* recoilDef =
        ionTable->GetIon(fRecoilZ, fRecoilA, 0.0);

    if (!ejectDef || !targDef || !recoilDef) {
        fLocalState.SetStatusChange(isAlive);
        return &fLocalState;
    }

    const G4double mEject  = ejectDef->GetPDGMass();
    const G4double mTarg   = targDef->GetPDGMass();
    const G4double mRecoil = recoilDef->GetPDGMass();

    const G4LorentzVector p4tot =
        projectile.Get4Momentum() + G4LorentzVector(0, 0, 0, mTarg);

    const G4ThreeVector projDir =
        projectile.Get4Momentum().vect().unit();

    // -------------------------------------------------------------------------
    // Try DDX-sampled lab ejectile + exact recoil reconstruction.
    // Reject samples that are kinematically impossible.
    // -------------------------------------------------------------------------
    for (G4int trial = 0; trial < 128; ++trial) {
        TALYSNeutronSample s = fSampler->Sample(fChannel, Einc_MeV);

        if (!s.valid) {
            break;
        }

        const G4double Eout =
            std::max(1.0e-12 * MeV, s.En_MeV * MeV);

        const G4double pout =
            SafeSqrt(Eout * (Eout + 2.0 * mEject));

        const G4ThreeVector dir =
            DirectionFromMuPhi(projDir, s.mu, s.phi);

        const G4LorentzVector p4eject(
            dir * pout,
            Eout + mEject
        );

        const G4LorentzVector p4recoil =
            p4tot - p4eject;

        const G4double m2recoil = p4recoil.m2();

        if (m2recoil <= 0.0) {
            continue;
        }

        const G4double mRecEffective =
            std::sqrt(m2recoil);

        const G4double excitation =
            mRecEffective - mRecoil;

        // Too much excitation means this sampled lab ejectile is inconsistent
        // with two-body kinematics for this incident energy.
        if (excitation < -0.1 * MeV || excitation > 30.0 * MeV) {
            continue;
        }

        const G4double exc = std::max(0.0, excitation);

        G4ParticleDefinition* recoilFinal =
            (exc > 0.01 * MeV)
            ? ionTable->GetIon(fRecoilZ, fRecoilA, exc)
            : recoilDef;

        if (!recoilFinal) {
            recoilFinal = recoilDef;
        }

        const G4double mRecoilFinal =
            recoilFinal->GetPDGMass();

        const G4double recoilKE =
            p4recoil.e() - mRecoilFinal;

        if (recoilKE < -0.05 * MeV) {
            continue;
        }

        // Add ejectile.
        fLocalState.AddSecondary(
            new G4DynamicParticle(ejectDef, dir, Eout)
        );

        // Always add the recoil residual nucleus.
        // This is essential for Geant4 energy conservation checks.
        G4ThreeVector recoilDir = -dir;

        if (p4recoil.vect().mag2() > 0.0) {
            recoilDir = p4recoil.vect().unit();
        }

        fLocalState.AddSecondary(
            new G4DynamicParticle(
                recoilFinal,
                recoilDir,
                std::max(0.0, recoilKE)
            )
        );

        return &fLocalState;
    }

    // -------------------------------------------------------------------------
    // Fallback: exact two-body CM kinematics.
    // This is important near threshold where DDX tables are missing or too sparse.
    // -------------------------------------------------------------------------
    const G4bool ok =
        AddTwoBodyFinalStateCM(
            p4tot,
            projDir,
            ejectDef,
            recoilDef
        );

    if (!ok) {
        // If kinematically forbidden, do not destroy the primary.
        // This should rarely be reached if the cross section is zero below threshold.
        fLocalState.SetStatusChange(isAlive);
    }

    return &fLocalState;
}
