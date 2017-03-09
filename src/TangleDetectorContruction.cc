#include "TangleDetectorConstruction.hh"

#include <G4SubtractionSolid.hh>
#include "G4Material.hh"

#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"

#include "G4VisAttributes.hh"

#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"

TangleDetectorConstruction::TangleDetectorConstruction()
        : G4VUserDetectorConstruction(),
          cAbsorberPV_A(0),
          NAbsorberPV_A(0),
          EAbsorberPV_A(0),
          SAbsorberPV_A(0),
          WAbsorberPV_A(0),
          cAbsorberPV_B(0),
          NAbsorberPV_B(0),
          EAbsorberPV_B(0),
          SAbsorberPV_B(0),
          WAbsorberPV_B(0),

          cFoilPV_A(0),

          fCheckOverlaps(true)
{ }

TangleDetectorConstruction::~TangleDetectorConstruction()
{ }

G4VPhysicalVolume* TangleDetectorConstruction::Construct()
{

  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();

}

void TangleDetectorConstruction::DefineMaterials(){

  G4NistManager *nistManager = G4NistManager::Instance();
  G4bool isotopes = false;

  nistManager->FindOrBuildMaterial("G4_AIR");
  G4double density = 2.700*g/cm3, a = 26.98*g/mole;

  G4Element*  Al = nistManager->FindOrBuildElement("Al" , isotopes);
  G4Material* Aluminium = new G4Material("Al", density, 1);
  Aluminium->AddElement(Al, 1);



  G4Element*  O = nistManager->FindOrBuildElement("O" , isotopes);
  G4Element* Si = nistManager->FindOrBuildElement("Si", isotopes);
  G4Element* Lu = nistManager->FindOrBuildElement("Lu", isotopes);
  G4Element* Y = nistManager->FindOrBuildElement("Y", isotopes);

  G4Material* LYSO = new G4Material("LYSO", 7.4*g/cm3, 4);
  LYSO->AddElement(Lu, 6);
  LYSO->AddElement(Si, 10);
  LYSO->AddElement(O , 50);
  LYSO->AddElement(Y, 14);
}


G4VPhysicalVolume *TangleDetectorConstruction::DefineVolumes() {

  // Geometry parameters

  //World size
  G4double worldSizeXY = 30. * cm;
  G4double worldSizeZ = 30. * cm;

  //Crystal Parameters
  G4double cryst_dX = 4.*mm, cryst_dY = 4.*mm, cryst_dZ = 22 * mm;
  G4double foilThickness = 0.5*mm;
  G4double outBox_dX = cryst_dX+foilThickness, outBox_dY = cryst_dY+foilThickness, outBox_dZ = cryst_dZ;

  G4double zPos = 2.1*cm;

  // Get materials
  G4Material* defaultMaterial = G4Material::GetMaterial("G4_AIR");
  G4Material* cryst_mat = G4Material::GetMaterial("LYSO");
  G4Material* foil_mat = G4Material::GetMaterial("Al");


  // Creation of World volume where the apparatus will be placed.
  G4VSolid *worldS
          = new G4Box("World",           // its name
                      worldSizeXY / 2, worldSizeXY / 2, worldSizeZ / 2); // its size

  G4LogicalVolume *worldLV
          = new G4LogicalVolume(
                  worldS,           // its solid
                  defaultMaterial,  // its material
                  "World");         // its name

  G4VPhysicalVolume *worldPV
          = new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(),  // at (0,0,0)
                  worldLV,          // its logical volume
                  "World",          // its name
                  0,                // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                  fCheckOverlaps);  // checking overlaps

  //
  //Crystals wrapped in foil.  These will form the basis of the simulated experimental apparatus
  G4double gap = 0. * mm;        //a gap for wrapping
  //Foil wrapping
  G4Box* outerBox = new G4Box("Outerbox", outBox_dX/2, outBox_dY/2, outBox_dZ/2);
  G4Box* innerBox = new G4Box("Innerbox", cryst_dX/2, cryst_dY/2, cryst_dZ/2+5*mm);

  G4SubtractionSolid* foilS = new G4SubtractionSolid("Foilwrapping", outerBox, innerBox);
  //G4double foil_xpos = cryst_dX+2*foilThickness+gap, foil_ypos = cryst_dY+2*foilThickness+gap;

  G4double cryst_xpos = cryst_dX + gap, cryst_ypos = cryst_dY + gap;
  G4Box* crystalS = new G4Box("crystal", cryst_dX/2, cryst_dY/2, cryst_dZ/2);

  G4LogicalVolume* crystalLV
          = new G4LogicalVolume(
                  crystalS,     // its solid
                  cryst_mat,  // its material
                  "lysoCrystal");   // its name

  G4LogicalVolume* foilLV
          = new G4LogicalVolume(
                  foilS,     // its solid
                  foil_mat,  // its material
                  "foilwrap");   // its name


  //Placement of detector A elements
  cAbsorberPV_A
          = new G4PVPlacement(
          0,                // no rotation
          G4ThreeVector(0, 0, zPos),  // at (0,0,0)
          crystalLV,          // its logical volume
          "crystalPV_CA",    // its name
          worldLV,          // its mother  volume
          false,            // no boolean operation
          0,                // copy number
          fCheckOverlaps);  // checking overlaps
  cFoilPV_A
          = new G4PVPlacement(
          0,
          G4ThreeVector(0, 0, zPos),
          foilLV,
          "foilPV_CA",
          worldLV,
          false,
          0,
          fCheckOverlaps);

  NAbsorberPV_A
          = new G4PVPlacement(
          0,                // no rotation
          G4ThreeVector(0, cryst_ypos+foilThickness, zPos), // at (0,0,0)
          crystalLV,          // its logical volume
          "crystalPV_NA",    // its name
          worldLV,          // its mother  volume
          false,            // no boolean operation
          1,                // copy number
          fCheckOverlaps);  // checking overlaps
  nFoilPV_A
          = new G4PVPlacement(
          0,
          G4ThreeVector(0, cryst_ypos+foilThickness, zPos),
          foilLV,
          "foilPV_NA",
          worldLV,
          false,
          1,
          fCheckOverlaps);

  NEAbsorberPV_A
          = new G4PVPlacement(
          0,                // no rotation
          G4ThreeVector(cryst_xpos+foilThickness, cryst_ypos+foilThickness, zPos), // at (0,0,0)
          crystalLV,          // its logical volume
          "crystalPV_NEA",    // its name
          worldLV,          // its mother  volume
          false,            // no boolean operation
          2,                // copy number
          fCheckOverlaps);  // checking overlaps
  neFoilPV_A
          = new G4PVPlacement(
          0,
          G4ThreeVector(cryst_xpos+foilThickness, cryst_ypos+foilThickness, zPos),
          foilLV,
          "foilPV_NEA",
          worldLV,
          false,
          2,
          fCheckOverlaps);


  EAbsorberPV_A
          = new G4PVPlacement(
          0,                // no rotation
          G4ThreeVector(cryst_xpos+foilThickness, 0, zPos),  // at (0,0,0)
          crystalLV,          // its logical volume
          "crystalPV_EA",    // its name
          worldLV,          // its mother  volume
          false,            // no boolean operation
          3,                // copy number
          fCheckOverlaps);  // checking overlaps
  eFoilPV_A
          = new G4PVPlacement(
          0,
          G4ThreeVector(cryst_xpos+foilThickness, 0, zPos),
          foilLV,
          "foilPV_EA",
          worldLV,
          false,
          3,
          fCheckOverlaps);

  SEAbsorberPV_A
          = new G4PVPlacement(
          0,                // no rotation
          G4ThreeVector(cryst_xpos+foilThickness, -cryst_ypos-foilThickness, zPos), // at (0,0,0)
          crystalLV,          // its logical volume
          "crystalPV_SEA",    // its name
          worldLV,          // its mother  volume
          false,            // no boolean operation
          4,                // copy number
          fCheckOverlaps);  // checking overlaps
  seFoilPV_A
          = new G4PVPlacement(
          0,
          G4ThreeVector(cryst_xpos+foilThickness, -cryst_ypos-foilThickness, zPos),
          foilLV,
          "foilPV_SEA",
          worldLV,
          false,
          4,
          fCheckOverlaps);

  SAbsorberPV_A
          = new G4PVPlacement(
          0,                // no rotation
          G4ThreeVector(0, -cryst_ypos-foilThickness, zPos),  // at (0,0,0)
          crystalLV,          // its logical volume
          "crystalPV_SA",    // its name
          worldLV,          // its mother  volume
          false,            // no boolean operation
          5,                // copy number
          fCheckOverlaps);  // checking overlaps
  sFoilPV_A
          = new G4PVPlacement(
          0,
          G4ThreeVector(0, -cryst_ypos-foilThickness, zPos),
          foilLV,
          "foilPV_SA",
          worldLV,
          false,
          5,
          fCheckOverlaps);

  SWAbsorberPV_A
          = new G4PVPlacement(
          0,                // no rotation
          G4ThreeVector(-cryst_xpos-foilThickness, -cryst_ypos-foilThickness, zPos), // at (0,0,0)
          crystalLV,          // its logical volume
          "crystalPV_SWA",    // its name
          worldLV,          // its mother  volume
          false,            // no boolean operation
          6,                // copy number
          fCheckOverlaps);  // checking overlaps
  swFoilPV_A
          = new G4PVPlacement(
          0,
          G4ThreeVector(-cryst_xpos-foilThickness, -cryst_ypos-foilThickness, zPos),
          foilLV,
          "foilPV_SWA",
          worldLV,
          false,
          6,
          fCheckOverlaps);

  WAbsorberPV_A
          = new G4PVPlacement(
          0,                // no rotation
          G4ThreeVector(-cryst_xpos-foilThickness, 0, zPos),  // at (0,0,0)
          crystalLV,          // its logical volume
          "crystalPV_WA",    // its name
          worldLV,          // its mother  volume
          false,            // no boolean operation
          7,                // copy number
          fCheckOverlaps);  // checking overlaps
  wFoilPV_A
          = new G4PVPlacement(
          0,
          G4ThreeVector(-cryst_xpos-foilThickness, 0, zPos),
          foilLV,
          "foilPV_WA",
          worldLV,
          false,
          7,
          fCheckOverlaps);

  NWAbsorberPV_A
          = new G4PVPlacement(
          0,                // no rotation
          G4ThreeVector(-cryst_xpos-foilThickness, cryst_ypos+foilThickness, zPos), // at (0,0,0)
          crystalLV,          // its logical volume
          "crystalPV_NWA",    // its name
          worldLV,          // its mother  volume
          false,            // no boolean operation
          8,                // copy number
          fCheckOverlaps);  // checking overlaps
  nwFoilPV_A
          = new G4PVPlacement(
          0,
          G4ThreeVector(-cryst_xpos-foilThickness, cryst_ypos+foilThickness, zPos),
          foilLV,
          "foilPV_NWA",
          worldLV,
          false,
          8,
          fCheckOverlaps);


  //////////////////////////////////////////////////////////////////////
  //Placement of detector B elements
  cAbsorberPV_B
          = new G4PVPlacement(
          0,                // no rotation
          G4ThreeVector(0, 0, -zPos),  // at (0,0,0)
          crystalLV,          // its logical volume
          "crystalPV_CB",    // its name
          worldLV,          // its mother  volume
          false,            // no boolean operation
          9,                // copy number
          fCheckOverlaps);  // checking overlaps
  cFoilPV_B
          = new G4PVPlacement(
          0,
          G4ThreeVector(0, 0, -zPos),
          foilLV,
          "foilPV_CB",
          worldLV,
          false,
          9,
          fCheckOverlaps);

  NAbsorberPV_B
          = new G4PVPlacement(
          0,                // no rotation
          G4ThreeVector(0, cryst_ypos+foilThickness, -zPos), // at (0,0,0)
          crystalLV,          // its logical volume
          "crystalPV_NB",    // its name
          worldLV,          // its mother  volume
          false,            // no boolean operation
          10,                // copy number
          fCheckOverlaps);  // checking overlaps
  nFoilPV_B
          = new G4PVPlacement(
          0,
          G4ThreeVector(0, cryst_ypos+foilThickness, -zPos),
          foilLV,
          "foilPV_NB",
          worldLV,
          false,
          10,
          fCheckOverlaps);

  NEAbsorberPV_B
          = new G4PVPlacement(
          0,                // no rotation
          G4ThreeVector(cryst_xpos+foilThickness, cryst_ypos+foilThickness, -zPos), // at (0,0,0)
          crystalLV,          // its logical volume
          "crystalPV_NEB",    // its name
          worldLV,          // its mother  volume
          false,            // no boolean operation
          11,                // copy number
          fCheckOverlaps);  // checking overlaps
  neFoilPV_B
          = new G4PVPlacement(
          0,
          G4ThreeVector(cryst_xpos+foilThickness, cryst_ypos+foilThickness, -zPos),
          foilLV,
          "foilPV_NEB",
          worldLV,
          false,
          11,
          fCheckOverlaps);


  EAbsorberPV_B
          = new G4PVPlacement(
          0,                // no rotation
          G4ThreeVector(cryst_xpos+foilThickness, 0, -zPos),  // at (0,0,0)
          crystalLV,          // its logical volume
          "crystalPV_EB",    // its name
          worldLV,          // its mother  volume
          false,            // no boolean operation
          12,                // copy number
          fCheckOverlaps);  // checking overlaps
  eFoilPV_B
          = new G4PVPlacement(
          0,
          G4ThreeVector(cryst_xpos+foilThickness, 0, -zPos),
          foilLV,
          "foilPV_EB",
          worldLV,
          false,
          12,
          fCheckOverlaps);

  SEAbsorberPV_B
          = new G4PVPlacement(
          0,                // no rotation
          G4ThreeVector(cryst_xpos+foilThickness, -cryst_ypos-foilThickness, -zPos), // at (0,0,0)
          crystalLV,          // its logical volume
          "crystalPV_SEB",    // its name
          worldLV,          // its mother  volume
          false,            // no boolean operation
          13,                // copy number
          fCheckOverlaps);  // checking overlaps
  seFoilPV_B
          = new G4PVPlacement(
          0,
          G4ThreeVector(cryst_xpos+foilThickness, -cryst_ypos-foilThickness, -zPos),
          foilLV,
          "foilPV_SEB",
          worldLV,
          false,
          13,
          fCheckOverlaps);

  SAbsorberPV_B
          = new G4PVPlacement(
          0,                // no rotation
          G4ThreeVector(0, -cryst_ypos-foilThickness, -zPos),  // at (0,0,0)
          crystalLV,          // its logical volume
          "crystalPV_SB",    // its name
          worldLV,          // its mother  volume
          false,            // no boolean operation
          14,                // copy number
          fCheckOverlaps);  // checking overlaps
  sFoilPV_B
          = new G4PVPlacement(
          0,
          G4ThreeVector(0, -cryst_ypos-foilThickness, zPos),
          foilLV,
          "foilPV_SB",
          worldLV,
          false,
          14,
          fCheckOverlaps);

  SWAbsorberPV_B
          = new G4PVPlacement(
          0,                // no rotation
          G4ThreeVector(-cryst_xpos-foilThickness, -cryst_ypos-foilThickness, -zPos), // at (0,0,0)
          crystalLV,          // its logical volume
          "crystalPV_SWB",    // its name
          worldLV,          // its mother  volume
          false,            // no boolean operation
          15,                // copy number
          fCheckOverlaps);  // checking overlaps
  swFoilPV_B
          = new G4PVPlacement(
          0,
          G4ThreeVector(-cryst_xpos-foilThickness, -cryst_ypos-foilThickness, -zPos),
          foilLV,
          "foilPV_SWB",
          worldLV,
          false,
          15,
          fCheckOverlaps);

  WAbsorberPV_B
          = new G4PVPlacement(
          0,                // no rotation
          G4ThreeVector(-cryst_xpos-foilThickness, 0, -zPos),  // at (0,0,0)
          crystalLV,          // its logical volume
          "crystalPV_WB",    // its name
          worldLV,          // its mother  volume
          false,            // no boolean operation
          16,                // copy number
          fCheckOverlaps);  // checking overlaps
  wFoilPV_B
          = new G4PVPlacement(
          0,
          G4ThreeVector(-cryst_xpos-foilThickness, 0, -zPos),
          foilLV,
          "foilPV_WB",
          worldLV,
          false,
          16,
          fCheckOverlaps);

  NWAbsorberPV_B
          = new G4PVPlacement(
          0,                // no rotation
          G4ThreeVector(-cryst_xpos-foilThickness, cryst_ypos+foilThickness, -zPos), // at (0,0,0)
          crystalLV,          // its logical volume
          "crystalPV_NWB",    // its name
          worldLV,          // its mother  volume
          false,            // no boolean operation
          17,                // copy number
          fCheckOverlaps);  // checking overlaps
  nwFoilPV_B
          = new G4PVPlacement(
          0,
          G4ThreeVector(-cryst_xpos-foilThickness, cryst_ypos+foilThickness, -zPos),
          foilLV,
          "foilPV_NWB",
          worldLV,
          false,
          17,
          fCheckOverlaps);

/////////////////////////////////////////////////////////////////////////////
  worldLV->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes *simpleBoxVisAtt = new G4VisAttributes(G4Colour(0, 1.0, 1.0));
  simpleBoxVisAtt->SetVisibility(true);
  crystalLV->SetVisAttributes(simpleBoxVisAtt);


  // Always return the physical World
  return worldPV;


}









