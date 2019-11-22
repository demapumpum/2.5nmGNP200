//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4PhysicalConstants.hh"
#include "G4Material.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0),
  solid(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
   
  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;
     
  // World
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_WATER");
  

      
  G4double world_prmax  = 500*um;      
  G4Orb* solidWorld =    
    new G4Orb("World",                      //its name
              world_prmax);		            //its size
                
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,         //its solid
                        world_mat,          //its material
                        "World");           //its name
						


  G4VPhysicalVolume* physWorld =                
    new G4PVPlacement(0,                    //no rotation
                    G4ThreeVector(),        //at position
                    logicWorld,             //its logical volume
                    "World",                //its name
                    0,                      //its mother  volume
                    false,                  //no boolean operation
                    0,                      //copy number
                    checkOverlaps);         //overlaps checking

  //Nanoparticle size parameters
  G4double shape1_prmin =  2.5*nm, shape1_prmax = 25*nm;

  // Detector
  //
  G4Material* det_mat = nist->FindOrBuildMaterial("G4_WATER");
      
  G4double shape2_prmin =  shape1_prmax, shape2_prmax = world_prmax;
  G4double shape2_pSPhi =  0.*deg, shape2_pDPhi = 360.*deg;
  G4double shape2_pSTheta = 0.*deg, shape2_pDTheta = 180.*deg;
  G4Sphere* solidDet =
    new G4Sphere("Detector", 
    shape2_prmin, shape2_prmax, shape2_pSPhi, shape2_pDPhi,
    shape2_pSTheta, shape2_pDTheta);
                
  G4LogicalVolume* logicDet =                         
    new G4LogicalVolume(solidDet,         //its solid
                        det_mat,          //its material
                        "Detector");      //its name


              
   new G4PVPlacement(0,                   //no rotation
                    G4ThreeVector(),      //at position
                    logicDet,             //its logical volume
                    "Detector",           //its name
                    logicWorld,           //its mother  volume
                    false,                //no boolean operation
                    0,                    //copy number
                    checkOverlaps);       //overlaps checking
  // Nanoparticle
  //
  G4Material* hollow = nist->FindOrBuildMaterial("G4_WATER");

  G4Material* Outernano_mat = nist->FindOrBuildMaterial("G4_Au");  
  
  if (shape1_prmin > 0){ 
  G4Orb* solidOuterNano =    
    new G4Orb("Nanoparticle", 
    shape1_prmax);
                      
  G4LogicalVolume* logicOuterNano =                         
    new G4LogicalVolume(solidOuterNano,      //its solid
                        Outernano_mat,       //its material
                        "Nanoparticle");     //its name
						

              
   new G4PVPlacement(0,                     //no rotation
                    G4ThreeVector(),         //at position
                    logicOuterNano,          //its logical volume
                    "Nanoparticle",          //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking 

  G4Orb* solidInnerNano =    
    new G4Orb("Nanoparticle", 
    shape1_prmin);
                      
  G4LogicalVolume* logicInnerNano =                         
    new G4LogicalVolume(solidInnerNano,       //its solid
                        hollow,               //its material
                        "Nanoparticle");      //its name
						

              
    solid = new G4PVPlacement(0,               //no rotation
                    G4ThreeVector(),          //at position
                    logicInnerNano,           //its logical volume
                    "InnerNanoparticle",      //its name
                    logicOuterNano,           //its mother  volume
                    false,                    //no boolean operation
                    0,                        //copy number
                    checkOverlaps);           //overlaps checking
  }else if (shape1_prmin == 0){

  G4Orb* solidOuterNano =    
    new G4Orb("Nanoparticle", 
    shape1_prmax);
                      
  G4LogicalVolume* logicOuterNano =                         
    new G4LogicalVolume(solidOuterNano,       //its solid
                        Outernano_mat,        //its material
                        "Nanoparticle");      //its name
						

              
   solid = new G4PVPlacement(0,               //no rotation
                    G4ThreeVector(),          //at position
                    logicOuterNano,           //its logical volume
                    "OuterNanoparticle",      //its name
                    logicWorld,               //its mother  volume
                    false,                    //no boolean operation
                    0,                        //copy number
                    checkOverlaps);           //overlaps checking 
  }
  

  // Set Detector as scoring volume

  fScoringVolume = logicDet;

  //always return the physical World
  //

  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
