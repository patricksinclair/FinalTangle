#include "AnnihilationPhotonsSteppingAction.hh"

#include "TangleRunAction.hh"

#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4MTRunManager.hh"
#include "G4EventManager.hh"
#include "G4TrackingManager.hh"
#include "G4SteppingManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Gamma.hh"
#include "G4Threading.hh"

#include "math.h"

AnnihilationPhotonsSteppingAction::AnnihilationPhotonsSteppingAction
        (TangleRunAction* runAction)
        : fpRunAction(runAction)
        , fAnnihilationPhotonFound1(false)
        , fComptonScatteringAnnihilationPhotonFound1(false)
        , fTrackID1(0)
        , fParentID1(0)
        , fTheta1(0.)
        , fPhi1(0.)

        ,secondPhotonFound(false)
        ,validTrack(true)
        ,good1(false), good2(false)
        ,hitCentreA(false), hitCentreB(false), hasScatteredA(false), hasScatteredB(false)
        ,hitNorthA(false), hitEastA(false), hitSouthA(false), hitWestA(false)
        ,hitNorthB(false), hitEastB(false), hitSouthB(false), hitWestB (false)
        ,stepCounter(0), stepCounter2(0)
        ,nIF_1(0), nIF_2(0)
        ,edepCA(0.), edepNA(0.), edepEA(0.), edepSA(0.), edepWA(0.)
        ,edepCB(0.), edepNB(0.), edepEB(0.), edepSB(0.), edepWB(0.)

{}



void AnnihilationPhotonsSteppingAction::BeginOfEventAction()
{
    fAnnihilationPhotonFound1 = false;
    fComptonScatteringAnnihilationPhotonFound1 = false;

    secondPhotonFound = false;
    validTrack = true;
    good1 = false; good2 = false;

    hitCentreA = false, hitCentreB = false;
    hasScatteredA = false, hasScatteredB = false;
    hitNorthA = false, hitEastA = false, hitSouthA = false, hitWestA = false;
    hitNorthB = false, hitEastB = false, hitSouthB = false, hitWestB = false;

    stepCounter = 0, stepCounter2 = 0;
    nIF_1 = 0, nIF_2 = 0;
    edepCA = 0., edepNA = 0., edepEA = 0., edepSA = 0., edepWA = 0.;
    edepCB = 0., edepNB = 0., edepEB = 0., edepSB = 0., edepWB = 0.;
    G4cout << "START OF EVENT" << G4endl;
}

void AnnihilationPhotonsSteppingAction::EndOfEventAction()
{
    G4cout << "END OF EVENT " << G4endl;
    G4cout << "   " << G4endl;

    if((good1 && good2) && (hitCentreA && hitCentreB)){

        crystEnergies.fEAC = edepCA, crystEnergies.fEAN = edepNA, crystEnergies.fEAE = edepEA;
        crystEnergies.fEAS = edepSA, crystEnergies.fEAW = edepWA, crystEnergies.fEBC = edepCB;
        crystEnergies.fEBN = edepNB, crystEnergies.fEBE = edepEB, crystEnergies.fEBS = edepSB;
        crystEnergies.fEBW = edepWB;

        if(hitNorthA) crystData.fPhiA = 0.;
        else if(hitEastA) crystData.fPhiA = 90.;
        else if(hitSouthA) crystData.fPhiA = 180.;
        else if(hitWestA) crystData.fPhiA = 270.;

        if(hitNorthB) crystData.fPhiB = 0.;
        if(hitEastB) crystData.fPhiB = 90.;
        if(hitSouthB) crystData.fPhiB = 180.;
        if(hitWestB) crystData.fPhiB = 270.;

        fpRunAction->RecordCrystalData(crystData);
        fpRunAction->RecordCrystalEnergies(crystEnergies);
    }
}

//#define AnnihilationPhotonsSteppingActionPrinting
#define AnnihilationPhotonsSteppingActionConsistencyCheck

namespace {
    void CalculateThetaPhi
            (const G4ThreeVector& v,
             const G4ThreeVector& z_axis,
                    // Output quantities
             G4ThreeVector& y_axis,
             G4ThreeVector& x_axis,
             G4double& cosTheta,
             G4double& theta,
             G4double& phi)
    {
        cosTheta = v*z_axis;
        theta = std::acos(cosTheta);
        // Make y' perpendicular to global x-axis.
        y_axis = (z_axis.cross(G4ThreeVector(1,0,0))).unit();
        x_axis = y_axis.cross(z_axis);
        const G4ThreeVector ontoXYPlane = v.cross(z_axis);
        // ontoXYPlane is a vector in the xy-plane, but perpendicular to the
        // projection of the scattered photon, so
        const G4double projection_x = -ontoXYPlane*y_axis;
        const G4double projection_y = ontoXYPlane*x_axis;
        phi = std::atan2(projection_y,projection_x);
    }
}

void AnnihilationPhotonsSteppingAction::UserSteppingAction(const G4Step* step)
{
    G4Track* track = step->GetTrack();
    G4StepPoint* postStepPoint = step->GetPostStepPoint();
    G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    G4double edep = (step->GetTotalEnergyDeposit())*pow(6.241509, 18);
    const G4VProcess* postProcessDefinedStep = postStepPoint->GetProcessDefinedStep();

    if (track->GetDefinition() != G4Gamma::Gamma()) return;

    G4cout << stepCounter << G4endl;

    //added by psinclair
    if(postProcessDefinedStep->GetProcessName() == "compt") stepCounter++;

    if(postProcessDefinedStep->GetProcessName() == "compt" &&
       fAnnihilationPhotonFound1 &&
       !secondPhotonFound) stepCounter2++;


    //changed from track->GetCurrentStepNumber() to stepCounter
    //if (stepCounter != 1) return;
    //G4cout << "This should be 1: n = " << track->GetCurrentStepNumber() << G4endl;
    //if(stepCounter > 1) G4cout <<"stepcounter failed" << G4endl;

    const G4VProcess* creatorProcess = track->GetCreatorProcess();
    if (creatorProcess == nullptr) return;
    if (creatorProcess->GetProcessName() != "annihil") return;

    G4StepPoint* preStepPoint = step->GetPreStepPoint();

    //stuff about postprocess defined step originally here. Have moved
    if (postProcessDefinedStep == nullptr) return;

    //if non-compton or non-transportation interaction happens, track is discounted (hopefully)
    //This introduced more dot product errors but at least now the no. of errors and
    //"good" readings are different
    if (postProcessDefinedStep->GetProcessName() != "compt" &&
        postProcessDefinedStep->GetProcessName() != "Transportation") {
        //G4cout << "ERROR: " << postProcessDefinedStep->GetProcessName()<< G4endl;
        validTrack = false;
        return;
    }


    if(hitCentreA){
        if(volume->GetName() == "crystalPV_CA"){
            edepCA += edep;
        } else if (volume->GetName() == "crystalPV_NA") {
            edepNA += edep;
            if(!hasScatteredA) {
                hitNorthA = true;
                hasScatteredA = true;
                G4cout << "SCATTERED NA " << G4endl;
            }
        } else if (volume->GetName() == "crystalPV_EA") {
            edepEA += edep;
            if(!hasScatteredA) {
                hitEastA = true;
                hasScatteredA = true;
                G4cout << "SCATTERED EA " << G4endl;
            }
        } else if (volume->GetName() == "crystalPV_SA") {
            edepSA += edep;
            if(!hasScatteredA) {
                hitSouthA = true;
                hasScatteredA = true;
                G4cout << "SCATTERED SA" << G4endl;
            }
        } else if (volume->GetName() == "crystalPV_WA") {
            edepWA += edep;
            if(!hasScatteredA) {
                hitWestA = true;
                hasScatteredA = true;
                G4cout << "SCATTERED WA" << G4endl;
            }
        }

    }else if(hitCentreB){
        if(volume->GetName() == "crystalPV_CB"){
            edepCB += edep;
        }else if (volume->GetName() == "crystalPV_NB") {
            edepNB += edep;
            if(!hasScatteredB) {
                hitNorthB = true;
                hasScatteredB = true;
                G4cout << "SCATTERED NB" << G4endl;
            }
        } else if (volume->GetName() == "crystalPV_EB") {
            edepEB += edep;
            if(!hasScatteredB) {
                hitEastB = true;
                hasScatteredB = true;
                G4cout << "SCATTERED EB" << G4endl;
            }
        } else if (volume->GetName() == "crystalPV_SB") {
            edepSB += edep;
            if(!hasScatteredB) {
                hitSouthB = true;
                hasScatteredB = true;
                G4cout << "SCATTERED SB" << G4endl;
            }
        } else if (volume->GetName() == "crystalPV_WB") {
            edepWB += edep;
            if(!hasScatteredB) {
                hitWestB = true;
                hasScatteredB = true;
                G4cout << "SCATTERED WB" << G4endl;
            }
        }
    }



    if(volume->GetName() == "crystalPV_CA" && !hitCentreA) {
        hitCentreA = true;
        edepCA += edep;
        G4cout << "HIT A" << G4endl;
    }
    if(volume->GetName() == "crystalPV_CB"  && !hitCentreB) {
        hitCentreB = true;
        edepCB += edep;
        G4cout << "HIT B" << G4endl;
    }









    // This is the frst step of an annihilation photon.
    ///////////////////////////////////////////////////////////////////////////
    if (!fAnnihilationPhotonFound1 && stepCounter == 1 && validTrack) {

        nIF_1++;

        G4cout << "Number of times in first if statement (should be 1): " << nIF_1 << G4endl;

        fAnnihilationPhotonFound1 = true;
        G4cout << "This shouldn't be 1 (track number): n = " << track->GetCurrentStepNumber() << G4endl;
        G4cout << "This should be 1 (stepCounter): n = " << stepCounter << G4endl;
        //transportation condition added psinclair
        if (postProcessDefinedStep->GetProcessName() != "compt" &&
            postProcessDefinedStep->GetProcessName() != "Transportation") return;

        fComptonScatteringAnnihilationPhotonFound1 = true;

        fTrackID1 = track->GetTrackID();
        fParentID1 = track->GetParentID();

        fPhoton1_z_axis = preStepPoint->GetMomentumDirection();

        G4ThreeVector photon1_y_axis;  // dummy, i.e., not used.
        G4ThreeVector photon1_x_axis;  // dummy, i.e., not used.
        CalculateThetaPhi
                (postStepPoint->GetMomentumDirection(),
                 fPhoton1_z_axis,
                 photon1_y_axis,
                 photon1_x_axis,
                 fCosTheta1,
                 fTheta1,
                 fPhi1);

#ifdef AnnihilationPhotonsSteppingActionPrinting
        G4cout
    << "\n  1st photon found: track ID: " << track->GetTrackID()
    << "\n  parent ID: " << track->GetParentID()
    << "\n  preStepPointPosition: " << preStepPoint->GetPosition()
    << "\n  postStepPointPosition: " << postStepPoint->GetPosition()
    << "\n  preStepPointMomentum: " << preStepPoint->GetMomentum()
    << "\n  postStepPointMomentum: " << postStepPoint->GetMomentum()
    << "\n  preStepPolarisation: " << preStepPoint->GetPolarization()
    << "\n  postStepPolarisation: " << postStepPoint->GetPolarization()
    << "\n  fTheta1: " << fTheta1
    << "\n  fPhi1: " << fPhi1
    << G4endl;
#endif  // AnnihilationPhotonsSteppingActionPrinting

        data.fTheta1 = fTheta1;
        data.fPhi1 = fPhi1;
        G4cout << "GOOD 1" << G4endl;
        good1 = true;
        ///////////////////////////////////////////////////////////////////////////
    } else if (fAnnihilationPhotonFound1 && stepCounter2 == 1 && !secondPhotonFound && validTrack) {
        ///////////////////////////////////////////////////////////////////////////
        nIF_2++;
        G4cout << "This shouldn't be 1 as well (else track number): n = " << track->GetCurrentStepNumber() << G4endl;
        G4cout << "This should be 1 as well (else stepCounter): n = " << stepCounter2 << G4endl;
        G4cout << "Number of times in second if statement (should be 1): " << nIF_2 << G4endl;

        // Second photon found
        //this boolean added by psinclair
        secondPhotonFound = true;
        // Unless BOTH annihilation photons undergo Compton scattering, do nothing.
        G4bool doNothing = false;
        if (!fComptonScatteringAnnihilationPhotonFound1) {
            doNothing = true;
#ifdef AnnihilationPhotonsSteppingActionPrinting
            G4cout << "First annihilation photon did not undergo Compton scattering." << G4endl;
#endif  // AnnihilationPhotonsSteppingActionPrinting
        }
        //Transportation condition added psinclair
        if (postProcessDefinedStep->GetProcessName() != "compt") {
            G4cout << "not second compton" << G4endl;
            doNothing = true;
#ifdef AnnihilationPhotonsSteppingActionPrinting
            G4cout << "Second annihilation photon did not undergo Compton scattering." << G4endl;
#endif  // AnnihilationPhotonsSteppingActionPrinting
        }
        if (doNothing) {
            //Reset for further possible annihilations in this event.
            fAnnihilationPhotonFound1 = false;
            fComptonScatteringAnnihilationPhotonFound1 = false;
            return;
        }
        G4cout << "MADE IT THIS FAR" << G4endl;

        const G4ThreeVector photon2_z_axis = preStepPoint->GetMomentumDirection();
        G4ThreeVector photon2_y_axis;
        G4ThreeVector photon2_x_axis;
        G4double originalCosTheta2;
        G4double originalTheta2;
        G4double originalPhi2;
        CalculateThetaPhi
                (postStepPoint->GetMomentumDirection(),
                 photon2_z_axis,
                 photon2_y_axis,
                 photon2_x_axis,
                 originalCosTheta2,
                 originalTheta2,
                 originalPhi2);

#ifdef AnnihilationPhotonsSteppingActionConsistencyCheck
        if (track->GetParentID() != fParentID1)
        {
            G4cout
                    << "\n  Annihilation photons do not have the same parent ID"
                    << "\n  track/parent IDs: " << fTrackID1 << '/' << fParentID1
                    << ',' << track->GetTrackID() << '/' << track->GetParentID()
                    << G4endl;

            //Reset for further possible annihilations in this event.
            fAnnihilationPhotonFound1 = false;
            fComptonScatteringAnnihilationPhotonFound1 = false;
            return;
        }
        const G4double dotProduct = preStepPoint->GetMomentumDirection().unit()*fPhoton1_z_axis;
        if (dotProduct > -0.999999) {
            G4cout <<
                   "\n  Annihilation photons not in opposite directions: dot product" << dotProduct
                   << G4endl;

            //Reset for further possible annihilations in this event.
            fAnnihilationPhotonFound1 = false;
            fComptonScatteringAnnihilationPhotonFound1 = false;
            return;
        }
#endif // AnnihilationPhotonsSteppingActionConsistencyCheck

#ifdef AnnihilationPhotonsSteppingActionPrinting
        G4cout
    << "\n  2nd photon found: track ID: " << track->GetTrackID()
    << "\n  parent ID: " << track->GetParentID()
    << "\n  preStepPointPosition: " << preStepPoint->GetPosition()
    << "\n  postStepPointPosition: " << postStepPoint->GetPosition()
    << "\n  preStepPointMomentum: " << preStepPoint->GetMomentum()
    << "\n  postStepPointMomentum: " << postStepPoint->GetMomentum()
    << "\n  preStepPolarisation: " << preStepPoint->GetPolarization()
    << "\n  postStepPolarisation: " << postStepPoint->GetPolarization()
    << "\n  originalTheta2: " << originalTheta2
    << "\n  originalPhi2: " << originalPhi2
    << G4endl;
#endif  // AnnihilationPhotonsSteppingActionPrinting

        // Calculate theta and phi of the Compton scatter of the second photon.
        // Draw the azimuthal angle from the entangled distribution:
        // A + B * cos(2*deltaPhi), or rather C + D * cos(2*deltaPhi), where
        // C = A / (A + |B|) and D = B / (A + |B|), so that maximum is 1.
        const G4double sin2Theta1 = 1.-fCosTheta1*fCosTheta1;
        const G4double sin2Theta2 = 1.-originalCosTheta2*originalCosTheta2;
        // Pryce and Ward
        const G4double A =
                ((std::pow(1.-fCosTheta1,3.))+2.)*(std::pow(1.-originalCosTheta2,3.)+2.)/
                ((std::pow(2.-fCosTheta1,3.)*std::pow(2.-originalCosTheta2,3.)));
        const G4double B = -(sin2Theta1*sin2Theta2)/
                           ((std::pow(2.-fCosTheta1,2.)*std::pow(2.-originalCosTheta2,2.)));
        const G4double C = A / (A + std::abs(B));
        const G4double D = B / (A + std::abs(B));
//    G4cout << "A,B,C,D: " << A << ',' << B  << ',' << C << ',' << D << G4endl;
//    // Snyder et al
//    const G4double& k0 = preStepPoint->GetKineticEnergy();
//    const G4double k1 = k0/(2.-fCosTheta1);
//    const G4double k2 = k0/(2.-originalCosTheta2);
//    const G4double gamma1 = k1/k0+k0/k1;
//    const G4double gamma2 = k2/k0+k0/k2;
//    const G4double A1 = gamma1*gamma2-gamma1*sin2Theta2-gamma2*sin2Theta1;
//    const G4double B1 = 2.*sin2Theta1*sin2Theta2;
//    // That's A1 + B1*sin2(deltaPhi) = A1 + B1*(0.5*(1.-cos(2.*deltaPhi).
//    const G4double ASnyder = A1 + 0.5*B1;
//    const G4double BSnyder = -0.5*B1;
//    const G4double CSnyder = ASnyder / (ASnyder + std::abs(BSnyder));
//    const G4double DSnyder = BSnyder / (ASnyder + std::abs(BSnyder));
//    G4cout << "A,B,C,D(Snyder): " << ASnyder << ',' << BSnyder << ',' << CSnyder << ',' << DSnyder << G4endl;

        // Sample delta phi
        G4double deltaPhi;
        const G4int maxCount = 999999;
        G4int iCount = 0;
        for (; iCount < maxCount; ++iCount) {
            deltaPhi = twopi * G4UniformRand();
            if (G4UniformRand() < C + D * cos(2.*deltaPhi)) break;
        }
        if (iCount >= maxCount ) {
            G4cout << "Random delta phi not found in " << maxCount << " tries - carry on anyway." << G4endl;
        }

        // Thus, the desired second photon azimuth
        G4double desiredPhi2 = deltaPhi - fPhi1;
        // Minus sign in above statement because of opposite coordinate system orientations.
        if (desiredPhi2 > pi) {
            desiredPhi2 -= twopi;
        }
        if (desiredPhi2 < -pi) {
            desiredPhi2 += twopi;
        }

        // Scattering angle is unchanged.
        const G4double desiredTheta2 = originalTheta2;

        G4ThreeVector newMomentumDirectionPrime;
        // In frame of second photon (denoted by "prime")
        newMomentumDirectionPrime.setRThetaPhi(1.,desiredTheta2,desiredPhi2);
        // Transform to global system
        // Some aliases
        const G4ThreeVector& v = newMomentumDirectionPrime;
        const G4ThreeVector& xp = photon2_x_axis;
        const G4ThreeVector& yp = photon2_y_axis;
        const G4ThreeVector& zp = photon2_z_axis;
        // In global system
        G4ThreeVector newMomentumDirection;
        newMomentumDirection.setX(v.x()*xp.x()+v.y()*yp.x()+v.z()*zp.x());
        newMomentumDirection.setY(v.x()*xp.y()+v.y()*yp.y()+v.z()*zp.y());
        newMomentumDirection.setZ(v.x()*xp.z()+v.y()*yp.z()+v.z()*zp.z());

#if defined AnnihilationPhotonsSteppingActionPrinting || defined AnnihilationPhotonsSteppingActionConsistencyCheck
        G4double newCosTheta2;
        G4double newTheta2;
        G4double newPhi2;
        CalculateThetaPhi
                (newMomentumDirection,
                 photon2_z_axis,
                 photon2_y_axis,
                 photon2_x_axis,
                 newCosTheta2,
                 newTheta2,
                 newPhi2);
#endif  // defined AnnihilationPhotonsSteppingActionPrinting || defined AnnihilationPhotonsSteppingActionConsistencyCheck

#ifdef AnnihilationPhotonsSteppingActionConsistencyCheck
        if (std::abs(newPhi2 - desiredPhi2) > 0.00001 || std::abs(newTheta2 - desiredTheta2) > 0.00001) {
            G4cout
                    << "\n  Inconsistent calculation of phi"
                    << "\n  originalTheta2: " << originalTheta2
                    << "\n  desiredTheta2: " << desiredTheta2
                    << "\n  newTheta2: " << newTheta2
                    << "\n  originalPhi2: " << originalPhi2
                    << "\n  desiredPhi2: " << desiredPhi2
                    << "\n  newPhi2: " << newPhi2
                    << G4endl;
        }
#endif // AnnihilationPhotonsSteppingActionConsistencyCheck

#ifdef AnnihilationPhotonsSteppingActionPrinting
        G4cout
    << "\n  originalTheta2: " << originalTheta2
    << "\n  desiredTheta2: " << desiredTheta2
    << "\n  newTheta2: " << newTheta2
    << "\n  originalPhi2: " << originalPhi2
    << "\n  desiredPhi2: " << desiredPhi2
    << "\n  newPhi2: " << newPhi2
    << G4endl;
#endif  // AnnihilationPhotonsSteppingActionPrinting

        track->SetMomentumDirection(newMomentumDirection);
        // And don't forget to do the same for the Compton electron.  <<<<<<<<< NEXT JOB

        data.fTheta2 = newTheta2;
        data.fPhi2 = newPhi2;
        fpRunAction->RecordData(data);
        G4cout << "GOOD 2" << G4endl;
        good2 = true;

        //Reset for further possible annihilations in this event.
        fAnnihilationPhotonFound1 = false;
        fComptonScatteringAnnihilationPhotonFound1 = false;

        //these were added by psinclair
        /* secondPhotonFound = false;
         stepCounter = 0;
         stepCounter2 = 0;
         */
        //G4cout << "END OF THE IFS" << G4endl;
    }
    ///////////////////////////////////////////////////////////////////////////

    //G4cout << "END OF THE CLASS" << G4endl;
    return;
}
