//==============================================================================
// DetectorConstruction.cc
//
// World   : G4_Galactic vacuum, ±5 cm box
// Catcher : BN_CATCHER (95% h-BN + 5% B₂O₃, 99% ¹¹B enrichment)
//           3×3 mm lateral × 1.5 mm thick, centred at origin (z = 0)
// Scoring : thin vacuum shell at r = 50 mm
//
// Protons are generated just upstream (z = −0.76 mm) pointing in +z.
// No pitcher, no spacer — the beam hits the catcher face directly.
//
// BN_CATCHER mass fractions (spec §2.3, locked 2026-05-12):
//   B-10  0.003947   B-11  0.429641
//   N-14  0.530052   N-15  0.002074
//   O-nat 0.034286
//==============================================================================
#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction()
{}

DetectorConstruction::~DetectorConstruction() {}

void DetectorConstruction::DefineMaterials()
{
    G4NistManager* nist = G4NistManager::Instance();
    fWorldMat = nist->FindOrBuildMaterial("G4_Galactic");

    G4Isotope* isoB10 = new G4Isotope("B10", 5, 10, 10.012937*g/mole);
    G4Isotope* isoB11 = new G4Isotope("B11", 5, 11, 11.009305*g/mole);
    G4Element* elB = new G4Element("BN_Boron", "B_BN", 2);
    elB->AddIsotope(isoB10,  1.0*perCent);
    elB->AddIsotope(isoB11, 99.0*perCent);

    G4Isotope* isoN14 = new G4Isotope("N14", 7, 14, 14.003074*g/mole);
    G4Isotope* isoN15 = new G4Isotope("N15", 7, 15, 15.000109*g/mole);
    G4Element* elN = new G4Element("BN_Nitrogen", "N_BN", 2);
    elN->AddIsotope(isoN14, 99.636*perCent);
    elN->AddIsotope(isoN15,  0.364*perCent);

    G4Element* elO = nist->FindOrBuildElement("O");

    fCatcherMat = new G4Material("BN_CATCHER", 2.12*g/cm3, 3);
    fCatcherMat->AddElement(elB, 0.433588);   // 0.003947+0.429641
    fCatcherMat->AddElement(elN, 0.532126);   // 0.530052+0.002074
    fCatcherMat->AddElement(elO, 0.034286);
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    DefineMaterials();

    // World: vacuum box 10x10x10 cm
    G4Box* worldSolid = new G4Box("World", 50.*mm, 50.*mm, 50.*mm);
    G4LogicalVolume* worldLV =
        new G4LogicalVolume(worldSolid, fWorldMat, "World");
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    G4VPhysicalVolume* worldPV =
        new G4PVPlacement(nullptr, {}, worldLV, "World", nullptr, false, 0, true);

    // Catcher: 3x3 mm x 1.5 mm BN_CATCHER, centred at z=0
    G4Box* catcherSolid = new G4Box("Catcher", 1.5*mm, 1.5*mm, 0.75*mm);
    fCatcherLV = new G4LogicalVolume(catcherSolid, fCatcherMat, "Catcher");
    auto* vis = new G4VisAttributes(G4Colour(0.2, 0.6, 1.0, 0.8));
    vis->SetForceSolid(true);
    fCatcherLV->SetVisAttributes(vis);
    new G4PVPlacement(nullptr, {}, fCatcherLV, "Catcher", worldLV, false, 0, true);

    // Virtual scoring sphere r=50 mm (thin shell)
    G4Sphere* sphSolid = new G4Sphere("ScoringSphere",
                                       49.9*mm, 50.*mm,
                                       0., 360.*deg, 0., 180.*deg);
    fScoringSphLV = new G4LogicalVolume(sphSolid, fWorldMat, "ScoringSphere");
    auto* svis = new G4VisAttributes(G4Colour(0., 1., 0., 0.12));
    svis->SetForceWireframe(true);
    fScoringSphLV->SetVisAttributes(svis);
    new G4PVPlacement(nullptr, {}, fScoringSphLV, "ScoringSphere",
                      worldLV, false, 0, true);

    G4cout <<
        "\n╔══════════════════════════════════════════╗\n"
        "║              Geometry                      ║\n"
        "╠════════════════════════════════════════════╣\n"
        "║  Catcher  : BN_CATCHER 3x3 mm x 1.5 mm     ║\n"
        "║             99% B-11, density=2.12 g/cm3   ║\n"
        "║             centred at origin z=0          ║\n"
        "║  World    : G4_Galactic vacuum, +/-5 cm    ║\n"
        "║  Scoring  : vacuum sphere r=50 mm          ║\n"
        "╚════════════════════════════════════════════╝\n\n";

    return worldPV;
}

void DetectorConstruction::ConstructSDandField() {}
