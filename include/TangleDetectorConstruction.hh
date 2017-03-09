#ifndef TangleDetectorConstruction_hh
#define TangleDetectorConstruction_hh 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
//class G4LogicalVolume;

class TangleDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    TangleDetectorConstruction();
    virtual ~TangleDetectorConstruction();

public:
    virtual G4VPhysicalVolume* Construct();

    //get methods
    const G4VPhysicalVolume* GetCAbsorberPV_A() const;
    const G4VPhysicalVolume* GetNAbsorberPV_A() const;
    const G4VPhysicalVolume* GetNEAbsorberPV_A() const;
    const G4VPhysicalVolume* GetEAbsorberPV_A() const;
    const G4VPhysicalVolume* GetSEAbsorberPV_A() const;
    const G4VPhysicalVolume* GetSAbsorberPV_A() const;
    const G4VPhysicalVolume* GetSWAbsorberPV_A() const;
    const G4VPhysicalVolume* GetWAbsorberPV_A() const;
    const G4VPhysicalVolume* GetNWAbsorberPV_A() const;

    const G4VPhysicalVolume* GetCAbsorberPV_B() const;
    const G4VPhysicalVolume* GetNAbsorberPV_B() const;
    const G4VPhysicalVolume* GetNEAbsorberPV_B() const;
    const G4VPhysicalVolume* GetEAbsorberPV_B() const;
    const G4VPhysicalVolume* GetSEAbsorberPV_B() const;
    const G4VPhysicalVolume* GetSAbsorberPV_B() const;
    const G4VPhysicalVolume* GetSWAbsorberPV_B() const;
    const G4VPhysicalVolume* GetWAbsorberPV_B() const;
    const G4VPhysicalVolume* GetNWAbsorberPV_B() const;


private:

    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();

    G4VPhysicalVolume*   cAbsorberPV_A; // the absorber physical volume
    G4VPhysicalVolume*   NAbsorberPV_A; // the absorber physical volume
    G4VPhysicalVolume*   EAbsorberPV_A; // the absorber physical volume
    G4VPhysicalVolume*   SAbsorberPV_A; // the absorber physical volume
    G4VPhysicalVolume*   WAbsorberPV_A; // the absorber physical volume
    G4VPhysicalVolume*   NEAbsorberPV_A; // the absorber physical volume
    G4VPhysicalVolume*   SEAbsorberPV_A; // the absorber physical volume
    G4VPhysicalVolume*   SWAbsorberPV_A; // the absorber physical volume
    G4VPhysicalVolume*   NWAbsorberPV_A; // the absorber physical volume

    G4VPhysicalVolume*   cAbsorberPV_B; // the absorber physical volume
    G4VPhysicalVolume*   NAbsorberPV_B; // the absorber physical volume
    G4VPhysicalVolume*   EAbsorberPV_B; // the absorber physical volume
    G4VPhysicalVolume*   SAbsorberPV_B; // the absorber physical volume
    G4VPhysicalVolume*   WAbsorberPV_B; // the absorber physical volume
    G4VPhysicalVolume*   NEAbsorberPV_B; // the absorber physical volume
    G4VPhysicalVolume*   SEAbsorberPV_B; // the absorber physical volume
    G4VPhysicalVolume*   SWAbsorberPV_B; // the absorber physical volume
    G4VPhysicalVolume*   NWAbsorberPV_B; // the absorber physical volume

    //The wrappings for the crystals
    G4VPhysicalVolume*   cFoilPV_A;
    G4VPhysicalVolume*   nFoilPV_A;
    G4VPhysicalVolume*   eFoilPV_A;
    G4VPhysicalVolume*   sFoilPV_A;
    G4VPhysicalVolume*   wFoilPV_A;
    G4VPhysicalVolume*   neFoilPV_A;
    G4VPhysicalVolume*   seFoilPV_A;
    G4VPhysicalVolume*   swFoilPV_A;
    G4VPhysicalVolume*   nwFoilPV_A;
    G4VPhysicalVolume*   cFoilPV_B;
    G4VPhysicalVolume*   nFoilPV_B;
    G4VPhysicalVolume*   eFoilPV_B;
    G4VPhysicalVolume*   sFoilPV_B;
    G4VPhysicalVolume*   wFoilPV_B;
    G4VPhysicalVolume*   neFoilPV_B;
    G4VPhysicalVolume*   seFoilPV_B;
    G4VPhysicalVolume*   swFoilPV_B;
    G4VPhysicalVolume*   nwFoilPV_B;


    G4VPhysicalVolume*   fGapPV;      // the gap physical volume

    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps

};

// inline functions
//Functions to return the crystals in array A
inline const G4VPhysicalVolume* TangleDetectorConstruction::GetCAbsorberPV_A() const {
    return cAbsorberPV_A;
}
inline const G4VPhysicalVolume* TangleDetectorConstruction::GetNAbsorberPV_A() const {
    return NAbsorberPV_A;
}
inline const G4VPhysicalVolume* TangleDetectorConstruction::GetNEAbsorberPV_A() const {
    return NEAbsorberPV_A;
}
inline const G4VPhysicalVolume* TangleDetectorConstruction::GetEAbsorberPV_A() const {
    return EAbsorberPV_A;
}
inline const G4VPhysicalVolume* TangleDetectorConstruction::GetSEAbsorberPV_A() const {
    return SEAbsorberPV_A;
}
inline const G4VPhysicalVolume* TangleDetectorConstruction::GetSAbsorberPV_A() const {
    return SAbsorberPV_A;
}
inline const G4VPhysicalVolume* TangleDetectorConstruction::GetSWAbsorberPV_A() const {
    return SWAbsorberPV_A;
}
inline const G4VPhysicalVolume* TangleDetectorConstruction::GetWAbsorberPV_A() const {
    return WAbsorberPV_A;
}
inline const G4VPhysicalVolume* TangleDetectorConstruction::GetNWAbsorberPV_A() const {
    return NWAbsorberPV_A;
}


//Functions to return the crystals in array B
inline const G4VPhysicalVolume* TangleDetectorConstruction::GetCAbsorberPV_B() const {
    return cAbsorberPV_B;
}
inline const G4VPhysicalVolume* TangleDetectorConstruction::GetNAbsorberPV_B() const {
    return NAbsorberPV_B;
}
inline const G4VPhysicalVolume* TangleDetectorConstruction::GetNEAbsorberPV_B() const {
    return NEAbsorberPV_B;
}
inline const G4VPhysicalVolume* TangleDetectorConstruction::GetEAbsorberPV_B() const {
    return EAbsorberPV_B;
}
inline const G4VPhysicalVolume* TangleDetectorConstruction::GetSEAbsorberPV_B() const {
    return SEAbsorberPV_B;
}
inline const G4VPhysicalVolume* TangleDetectorConstruction::GetSAbsorberPV_B() const {
    return SAbsorberPV_B;
}
inline const G4VPhysicalVolume* TangleDetectorConstruction::GetSWAbsorberPV_B() const {
    return SWAbsorberPV_B;
}
inline const G4VPhysicalVolume* TangleDetectorConstruction::GetWAbsorberPV_B() const {
    return WAbsorberPV_B;
}
inline const G4VPhysicalVolume* TangleDetectorConstruction::GetNWAbsorberPV_B() const {
    return NWAbsorberPV_B;
}


#endif

