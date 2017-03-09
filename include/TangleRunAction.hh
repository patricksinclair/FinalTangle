#ifndef TangleRunAction_hh
#define TangleRunAction_hh

#include "G4UserRunAction.hh"

#include <vector>

class G4Run;

class TangleRunAction : public G4UserRunAction
{
public:
    struct Data {
        Data(): fTheta1(0.), fPhi1(0.), fTheta2(0.), fPhi2(0.) {}
        Data(G4double theta1, G4double phi1, G4double theta2, G4double phi2):
                fTheta1(theta1), fPhi1(phi1), fTheta2(theta2), fPhi2(phi2) {}
        // Default copy constructor
        // Default assigment operator
        G4double fTheta1, fPhi1, fTheta2, fPhi2;
    };

    struct crystalData{
        crystalData(): fPhiA(0.), fPhiB(0.) {}
        crystalData(G4double phiA, G4double phiB):
                fPhiA(phiA), fPhiB(phiB) {}

        G4double  fPhiA, fPhiB;
    };

    struct  crystalEnergies{
        crystalEnergies(): fEAC(0.), fEAN(0.), fEAE(0.), fEAS(0.), fEAW(0.),
                           fEBC(0.), fEBN(0.), fEBE(0.), fEBS(0.), fEBW(0.) {}
        crystalEnergies(G4double EAC, G4double EAN, G4double EAE, G4double EAS, G4double EAW,
                        G4double EBC, G4double EBN, G4double EBE, G4double EBS, G4double EBW):
                fEAC(EAC), fEAN(EAN), fEAE(EAE), fEAS(EAS), fEAW(EAW),
                fEBC(EBC), fEBN(EBN), fEBE(EBE), fEBS(EBS), fEBW(EBW){}

        G4double fEAC, fEAN, fEAE, fEAS, fEAW, fEBC, fEBN, fEBE, fEBS, fEBW;
    };

    TangleRunAction();
    virtual ~TangleRunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    void RecordData(const Data& data) {fData.push_back(data);}
    void RecordCrystalData(const crystalData& cryDat) {fcrystalData.push_back(cryDat);}
    void RecordCrystalEnergies(const crystalEnergies& cryEng) {fcrystalEnergies.push_back(cryEng);}

private:
    static TangleRunAction* fpMasterRunAction;
    std::vector<Data> fData;
    static std::vector<std::vector<Data>*> fPointers;

    std::vector<crystalData> fcrystalData;
    static std::vector<std::vector<crystalData>*> fcrystalPointers;

    std::vector<crystalEnergies> fcrystalEnergies;
    static std::vector<std::vector<crystalEnergies>*> fcrystalEnergiesPointers;
};

#endif

