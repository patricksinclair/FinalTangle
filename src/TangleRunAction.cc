#include "TangleRunAction.hh"

#include "G4Run.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"
#include <cassert>
#include <fstream>

TangleRunAction *TangleRunAction::fpMasterRunAction = 0;
std::vector<std::vector<TangleRunAction::Data> *> TangleRunAction::fPointers;
std::vector<std::vector<TangleRunAction::crystalData> *> TangleRunAction::fcrystalPointers;
std::vector<std::vector<TangleRunAction::crystalEnergies> *> TangleRunAction::fcrystalEnergiesPointers;

namespace {
    G4Mutex runActionMutex = G4MUTEX_INITIALIZER;
}

TangleRunAction::TangleRunAction() {
    if (G4Threading::IsMasterThread()) {
        fpMasterRunAction = this;
    }
    fData.reserve(10000000);
    fcrystalData.reserve(10000000);
    fcrystalEnergies.reserve(10000000);
    G4AutoLock lock(&runActionMutex);
    fpMasterRunAction->fPointers.push_back(&fData);
    fpMasterRunAction->fcrystalPointers.push_back(&fcrystalData);
    fpMasterRunAction->fcrystalEnergiesPointers.push_back(&fcrystalEnergies);
}

TangleRunAction::~TangleRunAction() {}

void TangleRunAction::BeginOfRunAction(const G4Run *) {
    fData.clear();
}

namespace {
    std::ofstream outFile("outFile.csv");
    std::ofstream crystScattersFile("crystScatters.csv");
    std::ofstream crystEnergiesFile("crystEnergies.csv");
}

void TangleRunAction::EndOfRunAction(const G4Run *run) {
    G4int nofEvents = run->GetNumberOfEvent();
    if (nofEvents == 0) return;

    G4String runType;
    if (G4Threading::IsMasterThread()) {
        runType = "Global Run";
        // Workers finished...
        outFile << "#theta1, phi1, theta2, phi2";
        for (const auto &pointer: fPointers) {
            for (const auto &data: *pointer) {


                outFile
                        << "\n" << data.fTheta1 << ", " << data.fPhi1
                        << ", " << data.fTheta2 << ", " << data.fPhi2;

            }
        }
        outFile << std::endl;

        crystScattersFile << "#phiA, phiB, dPhi\n";
        for(const auto &crypointer: fcrystalPointers) {
            for(const auto &cryData: *crypointer) {

                G4double dPhi = cryData.fPhiB - cryData.fPhiA;
                if(dPhi < 0.) dPhi += 360.;

                crystScattersFile
                        << cryData.fPhiA << " " << cryData.fPhiB << " " << dPhi << "\n";
            }
        }
        crystScattersFile << std::endl;

        for(const auto &cryEngpointer: fcrystalEnergiesPointers){
            for(const auto &cryEngData: *cryEngpointer){


                crystEnergiesFile
                        << cryEngData.fEAC << " " << cryEngData.fEAN << " " << cryEngData.fEAE << " "
                        << cryEngData.fEAS << " " << cryEngData.fEAW << "\n" << cryEngData.fEBC << " "
                        << cryEngData.fEBN << " " << cryEngData.fEBE << " " << cryEngData.fEBS << " "
                        << cryEngData.fEBW <<"\n";


            }
        }
        crystEnergiesFile << std::endl;
    } else {
        runType = "Local Run-";
    }

    G4cout
            << "\n----------------------End of " << runType << "------------------------"
            << "\n The run consists of " << nofEvents << " events."
            << "\n------------------------------------------------------------"
            << G4endl;
}
